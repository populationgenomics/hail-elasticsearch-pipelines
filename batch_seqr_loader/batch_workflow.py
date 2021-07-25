#!/usr/bin/env python3

"""
Driver for loading data into SEQR for the CPG. See the README for more information.

    - 2021/04/16 Michael Franklin and Vlad Savelyev
"""

import json
import logging
import os
import shutil
import subprocess
import tempfile
import hashlib
from collections import defaultdict
from os.path import join, dirname, abspath, splitext
from typing import Optional, List, Tuple, Set, Iterable
import pandas as pd
import click
import hailtop.batch as hb
from analysis_runner import dataproc
from hailtop.batch.job import Job
from find_inputs import find_inputs
from utils import file_exists

GATK_VERSION = '4.2.0.0'
GATK_CONTAINER = (
    f'australia-southeast1-docker.pkg.dev/cpg-common/images/gatk:{GATK_VERSION}'
)
PICARD_CONTAINER = f'us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8'
BAZAM_CONTAINER = f'australia-southeast1-docker.pkg.dev/cpg-common/images/bazam:v2'
SOMALIER_CONTAINER = 'brentp/somalier:latest'
PEDDY_CONTAINER = 'quay.io/biocontainers/peddy:0.4.8--pyh5e36f6f_0'

NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS = 50
NUMBER_OF_DATAPROC_WORKERS = 50
NUMBER_OF_GENOMICS_DB_INTERVALS = 10

REF_BUCKET = 'gs://cpg-reference/hg38/v1'
REF_FASTA = join(REF_BUCKET, 'Homo_sapiens_assembly38.fasta')
DBSNP_VCF = join(REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf')
UNPADDED_INTERVALS = join(REF_BUCKET, 'hg38.even.handcurated.20k.intervals')
SOMALIER_SITES = join(REF_BUCKET, 'sites.hg38.vcf.gz')

DATAPROC_PACKAGES = [
    'seqr-loader',
    'click',
    'google',
    'slackclient',
    'fsspec',
    'sklearn',
    'gcloud',
]

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.option(
    '--gvcf',
    'gvcfs',
    multiple=True,
    help='Glob path to input GVCF files.C orresponding tbi index files are expected '
    'exist',
)
@click.option(
    '--bam',
    '--cram',
    'crams',
    multiple=True,
    help='Glob path to input BAM or CRAM files. Corresponding bai/crai indices are '
    'expected to existVariants will be called with GATK HaplotypeCaller',
)
@click.option(
    '--bam-to-realign',
    '--cram-to-realign',
    '--fastq-to-realign',
    '--data-to-realign',
    'data_to_realign',
    multiple=True,
    help='Glob path with input BAM/CRAM/FASTQ files that need re-alignment. '
    'For BAMs and CRAMs, corresponding bai/crai indices are expected to exist. '
    'Files will be (re-)aligned with BWA, and variants will be called with '
    'GATK HaplotypeCaller',
)
@click.option(
    '--seqr-dataset',
    'dataset_name',
    required=True,
    help='Name of the seqr dataset. If a genomics DB already exists for this dataset, '
    'it will be extended with the new samples. When loading into the ES, '
    'the name will be suffixed with the dataset version (set by --version)',
)
@click.option('--version', 'dataset_version', type=str, required=True)
@click.option(
    '-n',
    '--namespace',
    'output_namespace',
    type=click.Choice(['main', 'test', 'tmp']),
    help='The bucket namespace to write the results to',
)
@click.option(
    '--ped-file', 'ped_fpath', help='Pedigree file to verify the relatedness and sex.'
)
@click.option('--keep-scratch', 'keep_scratch', is_flag=True)
@click.option(
    '--overwrite/--reuse',
    'overwrite',
    is_flag=True,
    help='if an intermediate or a final file exists, skip running the code '
    'that generates it.',
)
@click.option('--dry-run', 'dry_run', is_flag=True)
@click.option('--disable-validation', 'disable_validation', is_flag=True)
@click.option(
    '--remap-tsv',
    'remap_path',
    help='Path to a TSV file with two columns: s and seqr_id.',
)
@click.option(
    '--subset-tsv',
    'subset_path',
    help='Path to a TSV file with one column of sample IDs: s.',
)
@click.option(
    '--make-checkpoints',
    'make_checkpoints',
    is_flag=True,
    help='Create checkpoints for intermediate Hail data',
)
@click.option(
    '--skip-ped-checks',
    'skip_ped_checks',
    is_flag=True,
    help='Skip checking provided sex and pedigree against the inferred one',
)
@click.option('--vep-block-size', 'vep_block_size')
def main(
    gvcfs: List[str],
    crams: List[str],
    data_to_realign: List[str],
    dataset_name: str,
    dataset_version: str,
    output_namespace: str,
    ped_fpath: str,
    keep_scratch: bool,
    overwrite: bool,
    dry_run: bool,
    disable_validation: bool,  # pylint: disable=unused-argument
    remap_path: str,
    subset_path: str,
    make_checkpoints: bool,
    skip_ped_checks: bool,
    namespace: Optional[str] = None,
    vep_block_size: Optional[int] = None,  # pylint: disable=unused-argument
):  # pylint: disable=missing-function-docstring
    if not (gvcfs or crams or data_to_realign):
        raise click.BadParameter(
            'Specify at least one of the input parameters '
            '(can be multiple and/or repeated): '
            '--gvcf, --cram, --bam, --data-to-realign'
        )

    project = 'seqr'
    billing_project = os.getenv('HAIL_BILLING_PROJECT') or 'seqr'

    if output_namespace in ['test', 'tmp']:
        tmp_bucket_suffix = 'test-tmp'
    else:
        tmp_bucket_suffix = 'main-tmp'
    work_bucket = f'gs://cpg-{project}-{tmp_bucket_suffix}/seqr-loader-tmp/{dataset_name}/{dataset_version}'
    hail_bucket = os.environ.get('HAIL_BUCKET')
    if not hail_bucket or keep_scratch:
        # Scratch files are large, so we want to use the temporary bucket for them
        hail_bucket = f'{work_bucket}/hail'
    logger.info(
        f'Starting hail Batch with the project {billing_project}, '
        f'bucket {hail_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=hail_bucket.replace('gs://', ''),
    )
    b = hb.Batch(
        f'Seqr loading pipeline for dataset "{dataset_name}"'
        f', in the namespace "{output_namespace}"',
        backend=backend,
    )

    if output_namespace in ['test', 'main']:
        output_suffix = output_namespace
        web_bucket_suffix = f'{output_namespace}-web'
    else:
        output_suffix = 'test-tmp'
        web_bucket_suffix = 'test-tmp'
    genomicsdb_bucket = f'gs://cpg-{project}-{output_suffix}/datasets/{dataset_name}/{dataset_version}/genomicsdbs'
    fingerprints_bucket = f'gs://cpg-{project}-{output_suffix}/datasets/{dataset_name}/{dataset_version}/fingerprints'
    web_bucket = f'gs://cpg-{project}-{web_bucket_suffix}/datasets/{dataset_name}/{dataset_version}'

    local_tmp_dir = tempfile.mkdtemp()

    samples_df, ped_fpath = find_inputs(
        gvcfs,
        crams,
        data_to_realign,
        local_tmp_dir=local_tmp_dir,
        work_bucket=work_bucket,
        ped_fpath=ped_fpath,
    )

    reference = b.read_input_group(
        base=REF_FASTA,
        fai=REF_FASTA + '.fai',
        dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
        + '.dict',
    )
    bwa_reference = b.read_input_group(
        base=REF_FASTA,
        fai=REF_FASTA + '.fai',
        dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
        + '.dict',
        sa=REF_FASTA + '.sa',
        amb=REF_FASTA + '.amb',
        bwt=REF_FASTA + '.bwt',
        ann=REF_FASTA + '.ann',
        pac=REF_FASTA + '.pac',
    )

    # Aligning available FASTQs
    align_fastq_jobs, samples_df = _make_realign_jobs(
        b=b,
        samples_df=samples_df,
        file_type='fastq_to_realign',
        reference=bwa_reference,
        work_bucket=work_bucket,
        overwrite=overwrite,
    )

    # Indexing input CRAMs, BAMs, GVCFs when index files are missing
    index_jobs, samples_df = _make_index_jobs(b, samples_df)

    # Realigning CRAMs
    realign_cram_jobs, samples_df = _make_realign_jobs(
        b=b,
        samples_df=samples_df,
        file_type='cram_to_realign',
        reference=bwa_reference,
        work_bucket=work_bucket,
        overwrite=overwrite,
        depends_on=index_jobs,
    )

    ped_check_j = None
    if not skip_ped_checks:
        ped_check_j, ped_fpath = _pedigree_checks(
            b=b,
            samples_df=samples_df,
            reference=reference,
            sites=b.read_input(SOMALIER_SITES),
            ped_file=b.read_input(ped_fpath),
            overwrite=overwrite,
            fingerprints_bucket=fingerprints_bucket,
            web_bucket=web_bucket,
            web_url=f'https://{namespace}-web.populationgenomics.org.au/{project}',
            depends_on=align_fastq_jobs + realign_cram_jobs,
        )

    genotype_jobs, samples_df = _make_haplotypecaller_jobs(
        b=b,
        samples_df=samples_df,
        reference=reference,
        work_bucket=work_bucket,
        overwrite=overwrite,
        depends_on=align_fastq_jobs
        + realign_cram_jobs
        + ([ped_check_j] if ped_check_j else []),
    )

    joint_genotype_job, gathered_vcf_path = _make_joint_genotype_jobs(
        b=b,
        genomicsdb_bucket=genomicsdb_bucket,
        samples_df=samples_df,
        reference=reference,
        dbsnp=DBSNP_VCF,
        work_bucket=work_bucket,
        local_tmp_dir=local_tmp_dir,
        overwrite=overwrite,
        depends_on=genotype_jobs,
    )

    annotated_mt_path = join(dirname(gathered_vcf_path), 'annotated.mt')
    if overwrite or not file_exists(annotated_mt_path):
        annotate_job = dataproc.hail_dataproc_job(
            b,
            f'batch_seqr_loader/scripts/make_annotated_mt.py '
            f'--source-path {gathered_vcf_path} '
            f'--dest-mt-path {annotated_mt_path} '
            f'--bucket {join(work_bucket, "seqr_load")} '
            + (f'--disable-validation ' if disable_validation else '')
            + (f'--make-checkpoints ' if make_checkpoints else '')
            + (f'--remap-tsv ' if remap_path else '')
            + (f'--subset-tsv ' if subset_path else '')
            + (f'--vep-block-size ' if vep_block_size else ''),
            max_age='8h',
            packages=DATAPROC_PACKAGES,
            num_secondary_workers=NUMBER_OF_DATAPROC_WORKERS,
            job_name='make_annotated_mt.py',
            vep='GRCh38',
            depends_on=[joint_genotype_job],
        )
    else:
        annotate_job = b.new_job('make_annotated_mt.py [reuse]')

    dataproc.hail_dataproc_job(
        b,
        f'batch_seqr_loader/scripts/load_to_es.py '
        f'--mt-path {annotated_mt_path} '
        f'--es-index {dataset_name}-{dataset_version} '
        f'--es-index-min-num-shards 1 '
        f'--genome-version GRCh38 '
        f'{"--prod" if namespace == "main" else ""}',
        max_age='8h',
        packages=DATAPROC_PACKAGES,
        num_secondary_workers=10,
        job_name='load_to_es.py',
        depends_on=[annotate_job],
        scopes=['cloud-platform'],
    )

    b.run(dry_run=dry_run, delete_scratch_on_exit=not keep_scratch)
    shutil.rmtree(local_tmp_dir)


def _pedigree_checks(
    b: hb.Batch,
    samples_df: pd.DataFrame,
    reference: hb.ResourceGroup,
    sites: hb.ResourceFile,
    ped_file: hb.ResourceFile,
    overwrite: bool,  # pylint: disable=unused-argument
    fingerprints_bucket: str,
    web_bucket: str,
    web_url: str,
    depends_on: Optional[List[Job]] = None,
) -> Tuple[Job, str]:
    """
    Add somalier and peddy based jobs that infer relatedness and sex, compare that
    to the provided PED file, and attempt to recover it. If unable to recover, cancel
    the further workflow jobs.

    Returns a job, and a bucket path to a fixed PED file if able to recover.
    """

    extract_jobs = []
    fp_file_by_sample = dict()
    for sn, input_path, input_index in zip(
        samples_df['s'], samples_df['file'], samples_df['index']
    ):
        fp_file_by_sample[sn] = join(fingerprints_bucket, f'{sn}.somalier')
        if can_reuse(fp_file_by_sample[sn], overwrite):
            extract_jobs.append(b.new_job(f'Somalier extract, {sn} [reuse]'))
        else:
            j = b.new_job(f'Somalier extract, {sn}')
            j.image(SOMALIER_CONTAINER)
            j.cpu(2)
            j.memory('highmem')  # ~ 4G/core ~ 8G
            if input_path.endswith('.bam'):
                j.storage(f'200G')
            elif input_path.endswith('.cram'):
                j.storage(f'50G')
            else:
                j.storage(f'10G')
            if depends_on:
                j.depends_on(*depends_on)

            input_file = b.read_input_group(
                base=input_path,
                index=input_index,
            )

            j.command(
                f"""set -ex
                
                somalier extract -d extracted/ --sites {sites} -f {reference.base} \\
                {input_file['base']}
                
                mv extracted/*.somalier {j.output_file}
                """
            )
            b.write_output(j.output_file, fp_file_by_sample[sn])
            extract_jobs.append(j)

    relate_j = b.new_job(f'Somalier relate')
    relate_j.image(SOMALIER_CONTAINER)
    relate_j.cpu(1)
    relate_j.memory('standard')  # ~ 4G/core ~ 4G
    # Size of one somalier file is 212K, so we add another G only if the number of
    # samples is >4k
    relate_j.storage(f'{1 + len(extract_jobs) // 4000 * 1}G')
    relate_j.depends_on(*extract_jobs)
    fp_files = [b.read_input(fp) for sn, fp in fp_file_by_sample.items()]
    relate_j.command(
        f"""set -e

        cat {ped_file} | grep -v Family.ID > samples.ped 

        somalier relate \\
        {' '.join(fp_files)} \\
        --ped samples.ped \\
        -o related \\
        --infer

        ls
        mv related.html {relate_j.output_html}
        mv related.pairs.tsv {relate_j.output_pairs}
        mv related.samples.tsv {relate_j.output_samples}
      """
    )

    # Copy somalier outputs to buckets
    sample_hash = hash_sample_names(samples_df['s'])
    prefix = join(fingerprints_bucket, sample_hash, 'somalier')
    somalier_samples_path = f'{prefix}.samples.tsv'
    somalier_pairs_path = f'{prefix}.pairs.tsv'
    b.write_output(relate_j.output_samples, somalier_samples_path)
    b.write_output(relate_j.output_pairs, somalier_pairs_path)
    # Copy somalier HTML to the web bucket
    rel_path = join('loader', sample_hash[:10], 'somalier.html')
    somalier_html_path = join(web_bucket, rel_path)
    somalier_html_url = f'{web_url}/{rel_path}'
    b.write_output(relate_j.output_html, somalier_html_path)

    check_j = b.new_job(f'Check relatedness and sex')
    check_j.image(PEDDY_CONTAINER)
    check_j.cpu(1)
    check_j.memory('standard')  # ~ 4G/core ~ 4G
    with open(join(dirname(abspath(__file__)), 'check_pedigree.py')) as f:
        script = f.read()
    check_j.command(
        f"""set -e
cat <<EOT >> check_pedigree.py
{script}
EOT
python check_pedigree.py \
--somalier-samples {relate_j.output_samples} \
--somalier-pairs {relate_j.output_pairs} \
{('--somalier-html ' + somalier_html_url) if somalier_html_url else ''}
    """
    )

    check_j.depends_on(relate_j)
    return check_j, somalier_samples_path


def _make_index_jobs(
    b: hb.Batch,
    samples_df: pd.DataFrame,
) -> Tuple[List[Job], pd.DataFrame]:
    jobs = []
    cram_df = samples_df[samples_df['type'].isin(['cram', 'cram_to_realign'])]
    for sn, fpath, index_fpath in zip(cram_df['s'], cram_df['file'], cram_df['index']):
        if index_fpath is None:
            job_name = f'Index alignment file, {sn}'
            file = b.read_input(fpath)
            j = b.new_job(job_name)
            jobs.append(j)
            j.image(BAZAM_CONTAINER)
            j.storage(('50G' if fpath.endswith('.cram') else '150G'))
            j.command(f'samtools index {file} {j.output_crai}')
            index_fpath = splitext(fpath)[0] + (
                '.crai' if fpath.endswith('.cram') else '.bai'
            )
            b.write_output(j.output_crai, index_fpath)
            samples_df.loc[sn, 'index'] = index_fpath
            jobs.append(j)
    gvcf_df = samples_df[samples_df['type'].isin(['gvcf'])]

    for sn, fpath, index_fpath in zip(cram_df['s'], gvcf_df['file'], gvcf_df['index']):
        if index_fpath is None:
            job_name = f'Index GVCF, {sn}'
            file = b.read_input(fpath)
            j = b.new_job(job_name)
            jobs.append(j)
            j.image('quay.io/biocontainers/tabix')
            j.storage('10G')
            j.command(f'tabix -p vcf {file} && ln {file}.tbi {j.output_tbi}')
            index_fpath = fpath + '.tbi'
            b.write_output(j.output_tbi, index_fpath)
            samples_df.loc[sn, 'index'] = index_fpath
            jobs.append(j)
    return jobs, samples_df


def _make_realign_jobs(
    b: hb.Batch,
    samples_df: pd.DataFrame,
    file_type: str,  # 'fastq_to_realign', 'cram_to_realign'
    reference: hb.ResourceGroup,
    work_bucket: str,
    overwrite: bool,
    depends_on: Optional[List[Job]] = None,
) -> Tuple[List[Job], pd.DataFrame]:
    """
    Takes all samples with a 'file' of 'type'='fastq_to_realign'|'cram_to_realign'
    in `samples_df`, runs BWA to realign reads again, and sets a new 'file'
    of 'type'='cram'.

    When the input is CRAM/BAM, uses Bazam to stream reads to BWA.
    """
    jobs = []
    realign_df = samples_df[samples_df['type'] == file_type]
    for sn, file1, file2, index in zip(
        realign_df['s'], realign_df['file'], realign_df['file2'], realign_df['index']
    ):
        job_name = f'BWA align, {sn}'
        output_cram_path = join(work_bucket, f'{sn}.cram')
        if can_reuse(output_cram_path, overwrite):
            jobs.append(b.new_job(f'{job_name} [reuse]'))
        else:
            assert (
                (file1.endswith('.cram') or file1.endswith('.bam')) and index or file2
            )
            if file1.endswith('.cram') or file1.endswith('.bam'):
                use_bazam = True
                file1 = b.read_input_group(base=file1, index=index)
            else:
                use_bazam = False
                file1 = b.read_input(file1)
                file2 = b.read_input(file2)

            j = b.new_job(job_name)
            jobs.append(j)
            j.image(BAZAM_CONTAINER)
            total_cpu = 16
            if use_bazam:
                bazam_cpu = 6 // 2
                bwa_cpu = 20 // 2
                bamsormadup_cpu = 6 // 2
            else:
                bazam_cpu = 0
                bwa_cpu = 24
                bamsormadup_cpu = 8
            j.cpu(total_cpu)
            j.memory('standard')
            j.storage('300G')
            j.declare_resource_group(
                output_cram={
                    'cram': '{root}.cram',
                    'crai': '{root}.crai',
                }
            )
            if depends_on:
                j.depends_on(*depends_on)

            if use_bazam:
                extract_fq_cmd = (
                    f'bazam -Xmx16g -Dsamjdk.reference_fasta={reference.base}'
                    f' -n{bazam_cpu} -bam {file1.base}'
                )
            else:
                extract_fq_cmd = ''

            rg_line = f'@RG\\tID:{sn}\\tSM:~{sn}'

            # BWA command options:
            # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
            # -p     smart pairing (ignoring in2.fq)
            # -v3    minimum score to output [30]
            # -t16   threads
            # -Y     use soft clipping for supplementary alignments
            # -R     read group header line such as '@RG\tID:foo\tSM:bar'
            # -M     mark shorter split hits as secondary
            j.command(
                f"""
set -o pipefail
set -ex

(while true; do df -h; pwd; du -sh $(dirname {j.output_cram.cram})/*; sleep 600; done) &

{extract_fq_cmd} | \\
bwa mem -K 100000000 {'-p' if use_bazam else ''} -v3 -t{bwa_cpu} -Y \\
  -R '{rg_line}' {reference.base} \\
  {'/dev/stdin' if use_bazam else file1} {'-' if use_bazam else file2} | \\
bamsormadup inputformat=sam threads={bamsormadup_cpu} SO=coordinate \\
  M={j.duplicate_metrics} outputformat=sam | \\
samtools view -T {reference.base} -O cram -o {j.output_cram.cram}

samtools index -@{total_cpu} {j.output_cram.cram} {j.output_cram.crai}

df -h; pwd; du -sh $(dirname {j.output_cram.cram})/*
            """
            )
            b.write_output(j.output_cram, splitext(output_cram_path)[0])
            b.write_output(
                j.duplicate_metrics,
                join(work_bucket, 'bwa', f'{sn}-duplicate-metrics.csv'),
            )

        samples_df.loc[sn, 'file'] = output_cram_path
        samples_df.loc[sn, 'index'] = splitext(output_cram_path)[0] + '.crai'
        samples_df.loc[sn, 'type'] = 'cram'

    return jobs, samples_df


def _make_haplotypecaller_jobs(
    b: hb.Batch,
    samples_df: pd.DataFrame,
    reference: hb.ResourceGroup,
    work_bucket: str,  # pylint: disable=unused-argument
    overwrite: bool,  # pylint: disable=unused-argument
    depends_on: Optional[List[Job]] = None,
) -> Tuple[List[Job], pd.DataFrame]:
    """
    Takes all samples with a 'file' of 'type'='bam' in `samples_df`,
    and runs HaplotypeCaller on them, and sets a new 'file' of 'type'='gvcf'

    HaplotypeCaller is run in an interval-based sharded way, with per-interval
    HaplotypeCaller jobs defined in a nested loop.
    """
    intervals_j = None

    merge_gvcf_jobs = []
    bams_df = samples_df[samples_df['type'] == 'cram']
    for sn, bam_fpath in zip(bams_df['s'], bams_df['file']):
        output_gvcf_path = join(work_bucket, f'{sn}.g.vcf.gz')
        if can_reuse(output_gvcf_path, overwrite):
            merge_gvcf_jobs.append(b.new_job(f'HaplotypeCaller, {sn} [reuse]'))
        else:
            haplotype_caller_jobs = []
            if intervals_j is None:
                intervals_j = _add_split_intervals_job(
                    b,
                    UNPADDED_INTERVALS,
                    NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
                    REF_FASTA,
                )
                if depends_on:
                    intervals_j.depends_on(*depends_on)
            for idx in range(NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS):
                haplotype_caller_jobs.append(
                    _add_haplotype_caller_job(
                        b,
                        bam_fpath,
                        interval=intervals_j.intervals[f'interval_{idx}'],
                        reference=reference,
                        sample_name=sn,
                        interval_idx=idx,
                        number_of_intervals=NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
                        depends_on=depends_on,
                    )
                )
            merge_gvcf_jobs.append(
                _add_merge_gvcfs_job(
                    b,
                    gvcfs=[j.output_gvcf for j in haplotype_caller_jobs],
                    output_gvcf_path=output_gvcf_path,
                    sample_name=sn,
                )
            )
        samples_df.loc[sn, 'file'] = output_gvcf_path
        samples_df.loc[sn, 'index'] = output_gvcf_path + '.tbi'
        samples_df.loc[sn, 'type'] = 'gvcf'

    return merge_gvcf_jobs, samples_df


def _make_joint_genotype_jobs(
    b: hb.Batch,
    genomicsdb_bucket: str,
    samples_df: pd.DataFrame,
    reference: hb.ResourceGroup,
    dbsnp: str,
    work_bucket: str,
    local_tmp_dir: str,
    overwrite: bool,
    depends_on: Optional[List[Job]] = None,
) -> Tuple[Job, str]:
    """
    Assumes all samples have a 'file' of 'type'='gvcf' in `samples_df`.
    Adds samples to the GenomicsDB and runs joint genotyping on them.
    Outputs a multi-sample VCF under `output_vcf_path`.
    """
    intervals_j = _add_split_intervals_job(
        b,
        UNPADDED_INTERVALS,
        NUMBER_OF_GENOMICS_DB_INTERVALS,
        REF_FASTA,
    )
    if depends_on:
        intervals_j.depends_on(*depends_on)

    genotype_vcf_jobs = []
    genotyped_vcfs = []

    samples_to_be_in_db_per_interval = dict()
    for idx in range(NUMBER_OF_GENOMICS_DB_INTERVALS):
        genomicsdb_gcs_path = join(
            genomicsdb_bucket,
            f'interval_{idx}_outof_{NUMBER_OF_GENOMICS_DB_INTERVALS}.tar',
        )

        import_gvcfs_job, samples_to_be_in_db = _add_import_gvcfs_job(
            b=b,
            genomicsdb_gcs_path=genomicsdb_gcs_path,
            samples_df=samples_df,
            work_bucket=work_bucket,
            local_tmp_dir=local_tmp_dir,
            interval=intervals_j.intervals[f'interval_{idx}'],
            interval_idx=idx,
            number_of_intervals=NUMBER_OF_GENOMICS_DB_INTERVALS,
        )

        samples_to_be_in_db_per_interval[idx] = samples_to_be_in_db
        samples_hash = hash_sample_names(samples_to_be_in_db)
        output_vcf_path = join(
            work_bucket, 'jointly-called', samples_hash, f'interval_{idx}.vcf.gz'
        )
        genotype_vcf_job = _add_gatk_genotype_gvcf_job(
            b,
            genomicsdb=b.read_input(genomicsdb_gcs_path),
            interval=intervals_j.intervals[f'interval_{idx}'],
            reference=reference,
            dbsnp=dbsnp,
            overwrite=overwrite,
            interval_idx=idx,
            number_of_samples=len(samples_to_be_in_db),
            number_of_intervals=NUMBER_OF_GENOMICS_DB_INTERVALS,
            output_vcf_path=output_vcf_path,
        )
        if import_gvcfs_job:
            genotype_vcf_job.depends_on(import_gvcfs_job)

        genotype_vcf_jobs.append(genotype_vcf_job)
        genotyped_vcfs.append(genotype_vcf_job.output_vcf)

    # We expect each DB to contain the same set of samples
    intervals_per_sample_sets = defaultdict(set)
    for interval_idx, samples_to_be_in_db in samples_to_be_in_db_per_interval.items():
        intervals_per_sample_sets[frozenset(samples_to_be_in_db)].add(interval_idx)
    if len(intervals_per_sample_sets) > 1:
        logger.critical(
            f'GenomicsDB for each interval is expected to contain the same '
            f'set of samples. Got: {intervals_per_sample_sets}'
        )

    samples_to_be_in_db = list(samples_to_be_in_db_per_interval.values())[0]
    samples_hash = hash_sample_names(samples_to_be_in_db)
    output_vcf_path = join(
        work_bucket, 'jointly-called', samples_hash, f'gathered.vcf.gz'
    )
    final_gathered_vcf_job = _add_final_gather_vcf_step(
        b,
        input_vcfs=genotyped_vcfs,
        overwrite=overwrite,
        output_vcf_path=output_vcf_path,
    )

    return final_gathered_vcf_job, output_vcf_path


def _add_split_intervals_job(
    b: hb.Batch,
    interval_list: str,
    scatter_count: int,
    ref_fasta: str,
) -> Job:
    """
    Split genome into intervals to parallelise GnarlyGenotyper.

    Returns: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job(f'Make {scatter_count} intervals')
    j.image(GATK_CONTAINER)
    java_mem = 3
    j.memory('standard')  # ~ 4G/core ~ 4G
    j.storage('16G')
    j.declare_resource_group(
        intervals={
            f'interval_{idx}': f'{{root}}/{str(idx).zfill(4)}-scattered.interval_list'
            for idx in range(scatter_count)
        }
    )

    j.command(
        f"""set -e

    # Modes other than INTERVAL_SUBDIVISION will produce an unpredicted number 
    # of intervals. But we have to expect exactly the {scatter_count} number of 
    # output files because our workflow is not dynamic.
    gatk --java-options -Xms{java_mem}g SplitIntervals \\
      -L {interval_list} \\
      -O {j.intervals} \\
      -scatter {scatter_count} \\
      -R {ref_fasta} \\
      -mode INTERVAL_SUBDIVISION
      """
    )
    # Could save intervals to a bucket here to avoid rerunning the job
    return j


def _add_haplotype_caller_job(
    b: hb.Batch,
    bam_fpath: str,
    interval: hb.ResourceFile,
    reference: hb.ResourceGroup,
    sample_name: str,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    depends_on: Optional[List[Job]] = None,
) -> Job:
    """
    Run HaplotypeCaller on an input BAM or CRAM, and output GVCF
    """
    job_name = 'HaplotypeCaller'
    if interval_idx is not None:
        job_name += f', {sample_name} {interval_idx}/{number_of_intervals}'

    j = b.new_job(job_name)
    j.image(GATK_CONTAINER)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage('8G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )
    if depends_on:
        j.depends_on(*depends_on)

    j.command(
        f"""set -e
    (while true; do df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']})/*; free -m; sleep 300; done) &

    gatk --java-options "-Xms{java_mem}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
      HaplotypeCaller \\
      -R {reference.base} \\
      -I {bam_fpath} \\
      -L {interval} \\
      -O {j.output_gvcf['g.vcf.gz']} \\
      -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \\
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \\
      -ERC GVCF \\

    df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']})/*; free -m
    """
    )
    return j


def _add_merge_gvcfs_job(
    b: hb.Batch,
    gvcfs: List[hb.ResourceGroup],
    output_gvcf_path: Optional[str],
    sample_name: str,
) -> Job:
    """
    Combine by-interval GVCFs into a single sample GVCF file
    """

    job_name = f'Merge {len(gvcfs)} GVCFs, {sample_name}'
    j = b.new_job(job_name)
    j.image(PICARD_CONTAINER)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage(f'{len(gvcfs) * 1.5 + 2}G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )

    input_cmd = ' '.join(f'INPUT={g["g.vcf.gz"]}' for g in gvcfs)

    j.command(
        f"""set -e

    (while true; do df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']})/*; free -m; sleep 300; done) &

    java -Xms{java_mem}g -jar /usr/picard/picard.jar \
      MergeVcfs {input_cmd} OUTPUT={j.output_gvcf['g.vcf.gz']}

    df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']})/*; free -m
      """
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, output_gvcf_path.replace('.g.vcf.gz', ''))
    return j


def _add_import_gvcfs_job(
    b: hb.Batch,
    genomicsdb_gcs_path: str,
    samples_df: pd.DataFrame,
    work_bucket: str,
    local_tmp_dir: str,
    interval: hb.ResourceFile,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    depends_on: Optional[List[Job]] = None,
) -> Tuple[Optional[Job], Set[str]]:
    """
    Add GVCFs to a genomics database (or create a new instance if it doesn't exist)
    Returns a Job, or None if no new samples to add
    """
    new_samples_df = samples_df[samples_df['type'] == 'gvcf']

    if file_exists(genomicsdb_gcs_path):
        # Check if sample exists in the DB already
        genomicsdb_metadata = join(local_tmp_dir, f'callset-{interval_idx}.json')
        # This command will download the DB metadata file locally.
        # The `-O` argument to `tar` means "write the file being extracted to the stdout",
        # and the file to be extracted is specified as a positional argument to `tar`.
        cmd = (
            f'gsutil cat {genomicsdb_gcs_path} | '
            f'tar -O --extract workspace/callset.json > {genomicsdb_metadata}'
        )
        logger.info(cmd)
        subprocess.run(cmd, check=False, shell=True)

        with open(genomicsdb_metadata) as f:
            db_metadata = json.load(f)
        samples_in_db = set(s['sample_name'] for s in db_metadata['callsets'])
        new_samples = set(new_samples_df.s)
        samples_to_add = new_samples - samples_in_db
        samples_to_skip = new_samples & samples_in_db
        if samples_to_skip:
            logger.warning(
                f'Samples {samples_to_skip} already exist in the DB '
                f'{genomicsdb_gcs_path}, skipping adding them'
            )
        if samples_to_add:
            logger.info(
                f'Will add samples {samples_to_add} into the DB '
                f'{genomicsdb_gcs_path}'
            )
        else:
            logger.warning(f'Nothing will be added into the DB {genomicsdb_gcs_path}')

        samples_to_add_df = new_samples_df[new_samples_df.s.isin(samples_to_add)]
        samples_will_be_in_db = samples_in_db | samples_to_add

        # Update existing DB
        genomicsdb_param = '--genomicsdb-update-workspace-path workspace'
        genomicsdb = b.read_input(genomicsdb_gcs_path)
        untar_genomicsdb_cmd = f'tar -xf {genomicsdb}'
        job_name = 'Adding to GenomicsDB'
    else:
        # Initiate new DB
        genomicsdb_param = '--genomicsdb-workspace-path workspace'
        untar_genomicsdb_cmd = ''
        job_name = 'Creating GenomicsDB'

        samples_to_add = set(new_samples_df.s)
        samples_to_skip = set()
        samples_to_add_df = new_samples_df
        samples_will_be_in_db = samples_to_add

    if not samples_to_add:
        return None, samples_will_be_in_db

    sample_map_fpath = join(work_bucket, 'work', 'sample_name.csv')
    samples_to_add_df[['s', 'file']].to_csv(
        sample_map_fpath, sep='\t', header=False, index=False
    )
    sample_name_map = b.read_input(sample_map_fpath)

    if interval_idx is not None:
        job_name += f' {interval_idx}/{number_of_intervals}'

    j = b.new_job(job_name)
    j.image(GATK_CONTAINER)
    j.cpu(16)
    java_mem = 14
    j.memory('lowmem')  # ~ 1G/core ~ 14.4G
    # 1G + 1G per sample divided by the number of intervals
    j.storage(f'{1 + len(samples_will_be_in_db) * 1 // number_of_intervals}G')
    if depends_on:
        j.depends_on(*depends_on)

    j.declare_resource_group(output={'tar': '{root}.tar'})
    j.command(
        f"""set -e

    # We've seen some GenomicsDB performance regressions related to intervals, 
    # so we're going to pretend we only have a single interval
    # using the --merge-input-intervals arg. There's no data in between since 
    # we didn't run HaplotypeCaller over those loci so we're not wasting any compute

    # The memory setting here is very important and must be several GiB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    
    # The batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    
    {untar_genomicsdb_cmd}

    (while true; do df -h; pwd; du -sh $(dirname {j.output['tar']})/*; free -m; sleep 300; done) &

    echo "Adding samples: {', '.join(samples_to_add)}"
    {f'echo "Skipping adding samples that are already in the DB: '
     f'{", ".join(samples_to_skip)}"' if samples_to_skip else ''}

    gatk --java-options -Xms{java_mem}g \
      GenomicsDBImport \
      {genomicsdb_param} \
      --batch-size 50 \
      -L {interval} \
      --sample-name-map {sample_name_map} \
      --reader-threads 14 \
      --merge-input-intervals \
      --consolidate

    do df -h; pwd; du -sh $(dirname {j.output['tar']})/*; free -m

    tar -cf {j.output['tar']} workspace

    do df -h; pwd; du -sh $(dirname {j.output['tar']})/*; free -m
    """
    )
    b.write_output(j.output, genomicsdb_gcs_path.replace('.tar', ''))
    return j, samples_will_be_in_db


def _add_gatk_genotype_gvcf_job(
    b: hb.Batch,
    genomicsdb: hb.ResourceFile,
    interval: hb.ResourceFile,
    reference: hb.ResourceGroup,
    dbsnp: str,
    overwrite: bool,
    number_of_samples: int,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Run joint-calling on all samples in a genomics database
    """
    job_name = 'Joint-genotype'
    if interval_idx is not None:
        job_name += f' {interval_idx}/{number_of_intervals}'

    if output_vcf_path and file_exists(output_vcf_path) and not overwrite:
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(GATK_CONTAINER)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 8G
    # 4G (fasta+fai+dict) + 1G per sample divided by the number of intervals
    j.storage(f'{4 + number_of_samples * 1 // number_of_intervals}G')
    j.declare_resource_group(
        output_vcf={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
        }
    )

    j.command(
        f"""set -e
        
    (while true; do df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})/*; free -m; sleep 300; done) &

    tar -xf {genomicsdb}

    df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})/*; free -m

    gatk --java-options -Xms{java_mem}g \\
      GenotypeGVCFs \\
      -R {reference.base} \\
      -O {j.output_vcf['vcf.gz']} \\
      -D {dbsnp} \\
      --only-output-calls-starting-in-intervals \\
      -V gendb://workspace \\
      -L {interval} \\
      --merge-input-intervals

    df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})/*; free -m
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_final_gather_vcf_step(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceGroup],
    overwrite: bool,
    output_vcf_path: str = None,
) -> Job:
    """
    Combines per-interval scattered VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    job_name = f'Gather {len(input_vcfs)} VCFs'

    if output_vcf_path and file_exists(output_vcf_path) and not overwrite:
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(GATK_CONTAINER)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage(f'{1 + len(input_vcfs) * 1}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    (while true; do df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})/*; free -m; sleep 300; done) &

    # --ignore-safety-checks makes a big performance difference so we include it in 
    # our invocation. This argument disables expensive checks that the file headers 
    # contain the same set of genotyped samples and that files are in order 
    # by position of first record.
    gatk --java-options -Xms{java_mem}g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      {input_cmdl} \\
      --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}
    
    df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})/*; free -m
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j


def can_reuse(fpath: str, overwrite: bool) -> bool:
    """
    Checks if the file `fpath` exists and we are not overwriting
    """
    if not file_exists(fpath):
        return False
    elif overwrite:
        logger.info(f'File {fpath} exists and will be overwritten')
        return False
    else:
        logger.info(f'Reusing existing {fpath}. Use --overwrite to overwrite')
        return True


def hash_sample_names(sample_names: Iterable[str]) -> str:
    """
    Return a unique hash string from a from a set of strings
    :param sample_names: set of strings
    :return: a string hash
    """
    for sn in sample_names:
        assert ' ' not in sn, sn
    return hashlib.sha224(' '.join(sorted(sample_names)).encode()).hexdigest()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
