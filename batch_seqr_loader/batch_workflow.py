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
import sys
import tempfile
import hashlib
from collections import defaultdict
from os.path import join, dirname, abspath, splitext
from typing import Optional, List, Tuple, Set, Iterable, Dict
import pandas as pd
import click
import hailtop.batch as hb
from analysis_runner import dataproc
from hailtop.batch.job import Job

from find_inputs import find_inputs, pull_inputs_from_sm_server
from vqsr import make_vqsr_jobs
import utils
from utils import file_exists


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.option(
    '--use-sm-server',
    'use_sm_server',
    is_flag=True,
    help='Query the sample metadata server for samples that need reprocessing, and'
    'update back statuses on progress',
)
@click.option(
    '--sm-server-db-name',
    'sm_server_db_name',
    default='seqr',
    help='Override the server metadata project/DB name',
)
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
    use_sm_server: bool,
    sm_server_db_name: str,
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

    if use_sm_server:
        if any([gvcfs, crams, data_to_realign]):
            logger.critical(
                'Can\'t mix direct inputs (--cram, --gvcf, etc) and --use-sm-server'
            )
            sys.exit(1)
        samples_df, ped_fpath = pull_inputs_from_sm_server(sm_server_db_name, ped_fpath)
    else:
        samples_df, ped_fpath = find_inputs(
            gvcfs,
            crams,
            data_to_realign,
            local_tmp_dir=local_tmp_dir,
            work_bucket=work_bucket,
            ped_fpath=ped_fpath,
        )

    reference = b.read_input_group(
        base=utils.REF_FASTA,
        fai=utils.REF_FASTA + '.fai',
        dict=utils.REF_FASTA.replace('.fasta', '')
        .replace('.fna', '')
        .replace('.fa', '')
        + '.dict',
    )
    bwa_reference = b.read_input_group(
        base=utils.REF_FASTA,
        fai=utils.REF_FASTA + '.fai',
        dict=utils.REF_FASTA.replace('.fasta', '')
        .replace('.fna', '')
        .replace('.fa', '')
        + '.dict',
        sa=utils.REF_FASTA + '.sa',
        amb=utils.REF_FASTA + '.amb',
        bwt=utils.REF_FASTA + '.bwt',
        ann=utils.REF_FASTA + '.ann',
        pac=utils.REF_FASTA + '.pac',
    )
    noalt_regions = b.read_input(join(utils.REF_BUCKET, 'noalt.bed'))

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
            sites=b.read_input(utils.SOMALIER_SITES),
            ped_file=b.read_input(ped_fpath),
            overwrite=overwrite,
            fingerprints_bucket=fingerprints_bucket,
            web_bucket=web_bucket,
            web_url=f'https://{output_namespace}-web.populationgenomics.org.au/{project}',
            depends_on=align_fastq_jobs + realign_cram_jobs,
        )

    produce_gvcf_jobs, samples_df = _make_produce_gvcfs_jobs(
        b=b,
        samples_df=samples_df,
        reference=reference,
        work_bucket=work_bucket,
        overwrite=overwrite,
        noalt_regions=noalt_regions,
        depends_on=align_fastq_jobs
        + realign_cram_jobs
        + ([ped_check_j] if ped_check_j else []),
    )

    joint_genotype_job, gathered_vcf_path = _make_joint_genotype_jobs(
        b=b,
        genomicsdb_bucket=genomicsdb_bucket,
        samples_df=samples_df,
        reference=reference,
        dbsnp=utils.DBSNP_VCF,
        work_bucket=work_bucket,
        local_tmp_dir=local_tmp_dir,
        overwrite=overwrite,
        depends_on=produce_gvcf_jobs,
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
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=utils.NUMBER_OF_DATAPROC_WORKERS,
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
        f'{"--prod" if output_namespace == "main" else ""}',
        max_age='8h',
        packages=utils.DATAPROC_PACKAGES,
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
            j.image(utils.SOMALIER_IMAGE)
            j.memory('standard')
            if input_path.endswith('.bam'):
                j.cpu(4)
                j.storage(f'200G')
            elif input_path.endswith('.cram'):
                j.cpu(4)
                j.storage(f'50G')
            else:
                j.cpu(2)
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
    relate_j.image(utils.SOMALIER_IMAGE)
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
    check_j.image(utils.PEDDY_IMAGE)
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
            j.image(utils.BAZAM_IMAGE)
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
                file2 or (file1.endswith('.cram') or file1.endswith('.bam')) and index
            )
            j = b.new_job(job_name)
            jobs.append(j)
            j.image(utils.BAZAM_IMAGE)
            total_cpu = 32
            use_bazam = file1.endswith('.cram') or file1.endswith('.bam')
            if use_bazam:
                bazam_cpu = 4
                bwa_cpu = 24
                bamsormadup_cpu = 4
            else:
                bazam_cpu = 0
                bwa_cpu = 32
                bamsormadup_cpu = 8
            j.cpu(total_cpu)
            j.memory('standard')
            j.storage('500G')
            j.declare_resource_group(
                output_cram={
                    'cram': '{root}.cram',
                    'crai': '{root}.crai',
                }
            )
            if depends_on:
                j.depends_on(*depends_on)

            if use_bazam:
                file1 = b.read_input_group(base=file1, index=index)
                r1_param = (
                    f'<(bazam -Xmx16g -Dsamjdk.reference_fasta={reference.base}'
                    f' -n{bazam_cpu} -bam {file1.base})'
                )
                r2_param = '-'
            else:
                files1 = [b.read_input(f1) for f1 in file1.split(',')]
                files2 = [b.read_input(f1) for f1 in file2.split(',')]
                r1_param = f'<(cat {" ".join(files1)})'
                r2_param = f'<(cat {" ".join(files2)})'

            rg_line = f'@RG\\tID:{sn}\\tSM:{sn}'

            # BWA command options:
            # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
            # -p     smart pairing (ignoring in2.fq)
            # -t16   threads
            # -Y     use soft clipping for supplementary alignments
            # -R     read group header line such as '@RG\tID:foo\tSM:bar'
            j.command(
                f"""
set -o pipefail
set -ex

(while true; do df -h; pwd; du -sh $(dirname {j.output_cram.cram}); sleep 600; done) &

bwa mem -K 100000000 {'-p' if use_bazam else ''} -t{bwa_cpu} -Y \\
  -R '{rg_line}' {reference.base} {r1_param} {r2_param} | \\
bamsormadup inputformat=sam threads={bamsormadup_cpu} SO=coordinate \\
  M={j.duplicate_metrics} outputformat=sam \\
  tmpfile=$(dirname {j.output_cram.cram})/bamsormadup-tmp | \\
samtools view -T {reference.base} -O cram -o {j.output_cram.cram}

samtools index -@{total_cpu} {j.output_cram.cram} {j.output_cram.crai}

df -h; pwd; du -sh $(dirname {j.output_cram.cram})
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


def _make_produce_gvcfs_jobs(
    b: hb.Batch,
    samples_df: pd.DataFrame,
    reference: hb.ResourceGroup,
    noalt_regions: hb.ResourceFile,
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
    merge_gvcf_job_by_sn = dict()
    cram_df = samples_df[samples_df['type'] == 'cram']
    for sn, bam_fpath in zip(cram_df['s'], cram_df['file']):
        called_gvcf_path = join(work_bucket, 'raw', f'{sn}.g.vcf.gz')
        if can_reuse(called_gvcf_path, overwrite):
            j = b.new_job(f'HaplotypeCaller, {sn} [reuse]')
        else:
            if intervals_j is None:
                intervals_j = _add_split_intervals_job(
                    b,
                    utils.UNPADDED_INTERVALS,
                    utils.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
                    utils.REF_FASTA,
                )
                if depends_on:
                    intervals_j.depends_on(*depends_on)
            haplotype_caller_jobs = []
            for idx in range(utils.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS):
                haplotype_caller_jobs.append(
                    _add_haplotype_caller_job(
                        b,
                        bam_fpath,
                        interval=intervals_j.intervals[f'interval_{idx}'],
                        reference=reference,
                        sample_name=sn,
                        interval_idx=idx,
                        number_of_intervals=utils.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
                        depends_on=depends_on,
                    )
                )
            j = _add_merge_gvcfs_job(
                b,
                gvcfs=[j.output_gvcf for j in haplotype_caller_jobs],
                output_gvcf_path=called_gvcf_path,
                sample_name=sn,
            )
        merge_gvcf_job_by_sn[sn] = j
        samples_df.loc[sn, 'file'] = called_gvcf_path
        samples_df.loc[sn, 'index'] = called_gvcf_path + '.tbi'
        samples_df.loc[sn, 'type'] = 'gvcf'

    gvcf_df = samples_df[samples_df['type'] == 'gvcf']
    jobs = []
    for sn, gvcf_fpath, tbi_path in zip(gvcf_df['s'], gvcf_df['file'], gvcf_df['index']):
        called_gvcf_path = join(work_bucket, f'{sn}.g.vcf.gz')
        reblock_gvcf_job = _add_reblock_gvcf_job(
            b,
            b.read_input_group(**{'g.vcf.gz': gvcf_fpath, 'g.vcf.gz.tbi': tbi_path}),
            overwrite,
        )
        if merge_gvcf_job_by_sn.get(sn):
            reblock_gvcf_job.depends_on(merge_gvcf_job_by_sn[sn])
        subset_to_noalt_job = _add_subset_noalt_step(
            b,
            reblock_gvcf_job.output_gvcf,
            noalt_regions,
            overwrite,
        )
        jobs.append(subset_to_noalt_job)
        samples_df.loc[sn, 'file'] = called_gvcf_path
        samples_df.loc[sn, 'index'] = called_gvcf_path + '.tbi'
        samples_df.loc[sn, 'type'] = 'gvcf'

    return jobs, samples_df


def _add_reblock_gvcf_job(
    b: hb.Batch,
    input_gvcf: hb.ResourceGroup,
    overwrite: bool,
    output_gvcf_path: Optional[str] = None,
) -> Job:
    """
    Runs ReblockGVCF to annotate with allele-specific VCF INFO fields
    required for recalibration
    """
    job_name = 'ReblockGVCF'
    if output_gvcf_path and file_exists(output_gvcf_path) and not overwrite:
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'30G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )

    j.command(
        f"""
    gatk --java-options "-Xms{mem_gb - 1}g" \\
        ReblockGVCF \\
        -V {input_gvcf['g.vcf.gz']} \\
        -do-qual-approx \\
        -O {j.output_gvcf['g.vcf.gz']} \\
        --create-output-variant-index true"""
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, output_gvcf_path.replace('.g.vcf.gz', ''))
    return j


def _add_subset_noalt_step(
    b: hb.Batch,
    input_gvcf: hb.ResourceGroup,
    noalt_regions: str,
    overwrite: bool,
    output_gvcf_path: Optional[str] = None,
) -> Job:
    """
    1. Subset GVCF to main chromosomes to avoid downstream errors
    2. Removes the DS INFO field that is added to some HGDP GVCFs to avoid errors
       from Hail about mismatched INFO annotations
    """
    job_name = 'SubsetToNoalt'
    if output_gvcf_path and file_exists(output_gvcf_path) and not overwrite:
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.BCFTOOLS_IMAGE)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'30G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )
    j.command(
        f"""set -e

    bcftools view \\
        {input_gvcf['g.vcf.gz']} \\
        -T {noalt_regions} \\
        | bcftools annotate -x INFO/DS \\
        -o {j.output_gvcf['g.vcf.gz']} \\
        -Oz

    bcftools index --tbi {j.output_gvcf['g.vcf.gz']}
        """
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, output_gvcf_path.replace('.g.vcf.gz', ''))
    return j


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
    is_small_callset = len(samples_df) < 1000
    # 1. For small callsets, we don't apply the ExcessHet filtering.
    # 2. For small callsets, we gather the VCF shards and collect QC metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics
    # collected from them.
    is_huge_callset = len(samples_df) >= 100000
    # For huge callsets, we allocate more memory for the SNPs Create Model step

    intervals_j = _add_split_intervals_job(
        b,
        utils.UNPADDED_INTERVALS,
        utils.NUMBER_OF_GENOMICS_DB_INTERVALS,
        utils.REF_FASTA,
    )
    if depends_on:
        intervals_j.depends_on(*depends_on)

    samples_to_be_in_db_per_interval = dict()
    genomics_gcs_path_per_interval = dict()
    import_gvcfs_job_per_interval = dict()
    for idx in range(utils.NUMBER_OF_GENOMICS_DB_INTERVALS):
        genomics_gcs_path_per_interval[idx] = join(
            genomicsdb_bucket,
            f'interval_{idx}_outof_{utils.NUMBER_OF_GENOMICS_DB_INTERVALS}.tar',
        )

        import_gvcfs_job, samples_to_be_in_db = _add_import_gvcfs_job(
            b=b,
            genomicsdb_gcs_path=genomics_gcs_path_per_interval[idx],
            samples_df=samples_df,
            work_bucket=work_bucket,
            local_tmp_dir=local_tmp_dir,
            interval=intervals_j.intervals[f'interval_{idx}'],
            interval_idx=idx,
            number_of_intervals=utils.NUMBER_OF_GENOMICS_DB_INTERVALS,
        )
        import_gvcfs_job_per_interval[idx] = import_gvcfs_job
        samples_to_be_in_db_per_interval[idx] = samples_to_be_in_db

    # Post-import-genomicsdb check:
    # we expect each DB to contain the same set of samples
    intervals_per_sample_sets = defaultdict(set)
    for interval_idx, samples_to_be_in_db in samples_to_be_in_db_per_interval.items():
        intervals_per_sample_sets[frozenset(samples_to_be_in_db)].add(interval_idx)
    if len(intervals_per_sample_sets) > 1:
        logger.critical(
            f'GenomicsDB for each interval is expected to contain the same '
            f'set of samples. Got: {intervals_per_sample_sets}'
        )

    make_site_only_jobs = []
    scattered_vcf_by_interval: Dict[int, hb.ResourceGroup] = dict()

    for idx in range(utils.NUMBER_OF_GENOMICS_DB_INTERVALS):
        samples_hash = hash_sample_names(samples_to_be_in_db_per_interval[idx])
        site_only_vcf_path = join(
            work_bucket, 'jointly-called', samples_hash, f'interval_{idx}.vcf.gz'
        )
        if site_only_vcf_path and file_exists(site_only_vcf_path) and not overwrite:
            make_site_only_jobs.append(b.new_job('Joint genotyping [reuse]'))
            scattered_vcf_by_interval[idx] = b.read_input_group(
                **{
                    'vcf.gz': site_only_vcf_path,
                    'vcf.gz.tbi': site_only_vcf_path + '.tbi',
                }
            )
        else:
            genotype_vcf_job = _add_genotype_gvcfs_job(
                b,
                genomicsdb=b.read_input(genomics_gcs_path_per_interval[idx]),
                interval=intervals_j.intervals[f'interval_{idx}'],
                reference=reference,
                dbsnp=dbsnp,
                overwrite=overwrite,
                interval_idx=idx,
                number_of_samples=len(samples_to_be_in_db_per_interval[idx]),
                number_of_intervals=utils.NUMBER_OF_GENOMICS_DB_INTERVALS,
            )
            if import_gvcfs_job_per_interval.get(idx):
                genotype_vcf_job.depends_on(import_gvcfs_job_per_interval.get(idx))

            vcf = genotype_vcf_job.output_vcf
            if not is_small_callset:
                exccess_filter_job = _add_exccess_het_filter(
                    b,
                    input_vcf=vcf,
                    overwrite=overwrite,
                    interval=intervals_j.intervals[f'interval_{idx}'],
                )
                vcf = exccess_filter_job.output_vcf
            make_site_only_job = _add_make_sites_only_job(
                b,
                input_vcf=vcf,
                overwrite=overwrite,
            )
            make_site_only_jobs.append(make_site_only_job)
            scattered_vcf_by_interval[idx] = make_site_only_job.output_vcf

    samples_to_be_in_db = list(samples_to_be_in_db_per_interval.values())[0]
    samples_hash = hash_sample_names(samples_to_be_in_db)
    gathered_vcf_path = join(
        work_bucket, 'jointly-called', samples_hash, f'gathered.vcf.gz'
    )
    final_gathered_vcf_job = _add_final_gather_vcf_job(
        b,
        input_vcfs=list(scattered_vcf_by_interval.values()),
        overwrite=overwrite,
        output_vcf_path=gathered_vcf_path,
    )

    vqsr_bucket = join(work_bucket, 'vqsrred', samples_hash)
    vqsred_vcf_path = join(vqsr_bucket, f'output.vcf.gz')
    if overwrite or not utils.file_exists(vqsred_vcf_path):
        vqsr_job = make_vqsr_jobs(
            b,
            input_vcf_gathered=gathered_vcf_path,
            input_vcfs_scattered=list(scattered_vcf_by_interval.values()),
            is_small_callset=is_small_callset,
            is_huge_callset=is_huge_callset,
            work_bucket=vqsr_bucket,
            web_bucket=vqsr_bucket,
            depends_on=[final_gathered_vcf_job],
            vqsr_params_d={
                'snp_filter_level': 99.7,
                'indel_filter_level': 99.0,
            },
            intervals=intervals_j.intervals,
            scatter_count=utils.NUMBER_OF_GENOMICS_DB_INTERVALS,
            output_vcf_path=vqsred_vcf_path,
        )
    else:
        vqsr_job = b.new_job('VQSR [reuse]')

    return vqsr_job, vqsred_vcf_path


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
    j.image(utils.GATK_IMAGE)
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
    output_gvcf_path: Optional[str] = None,
    overwrite: bool = False,
) -> Job:
    """
    Run HaplotypeCaller on an input BAM or CRAM, and output GVCF
    """
    job_name = 'HaplotypeCaller'
    if interval_idx is not None:
        job_name += f', {sample_name} {interval_idx}/{number_of_intervals}'
    if output_gvcf_path and not overwrite and utils.file_exists(output_gvcf_path):
        return b.new_job(f'{job_name} [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
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
    (while true; do df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m; sleep 300; done) &

    gatk --java-options "-Xms{java_mem}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
      HaplotypeCaller \\
      -R {reference.base} \\
      -I {bam_fpath} \\
      -L {interval} \\
      -O {j.output_gvcf['g.vcf.gz']} \\
      -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \\
      -GQB 20 \
      -ERC GVCF \\

    df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m
    """
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, j.output_gvcf)
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
    j.image(utils.PICARD_IMAGE)
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

    (while true; do df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m; sleep 300; done) &

    java -Xms{java_mem}g -jar /usr/picard/picard.jar \
      MergeVcfs {input_cmd} OUTPUT={j.output_gvcf['g.vcf.gz']}

    df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m
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
    j.image(utils.GATK_IMAGE)
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

    (while true; do df -h; pwd; du -sh $(dirname {j.output['tar']}); free -m; sleep 300; done) &

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

    df -h; pwd; du -sh $(dirname {j.output['tar']}); free -m

    tar -cf {j.output['tar']} workspace

    df -h; pwd; du -sh $(dirname {j.output['tar']}); free -m
    """
    )
    b.write_output(j.output, genomicsdb_gcs_path.replace('.tar', ''))
    return j, samples_will_be_in_db


def _add_genotype_gvcfs_job(
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
    job_name = 'Joint genotyping: GenotypeGVCFs'
    if interval_idx is not None:
        job_name += f' {interval_idx}/{number_of_intervals}'

    if output_vcf_path and utils.file_exists(output_vcf_path) and not overwrite:
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
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
        
    (while true; do df -h; pwd; free -m; sleep 300; done) &

    tar -xf {genomicsdb}

    df -h; pwd; free -m

    gatk --java-options -Xms{java_mem}g \\
      GenotypeGVCFs \\
      -R {reference.base} \\
      -O {j.output_vcf['vcf.gz']} \\
      -D {dbsnp} \\
      --only-output-calls-starting-in-intervals \\
      -V gendb://workspace \\
      -L {interval} \\
      --merge-input-intervals

    df -h; pwd; free -m
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_gnarly_genotyper_job(
    b: hb.Batch,
    genomicsdb: hb.ResourceFile,
    reference: hb.ResourceGroup,
    dbsnp: hb.ResourceGroup,
    overwrite: bool,
    number_of_samples: int,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    interval: Optional[hb.ResourceGroup] = None,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Runs GATK GnarlyGenotyper on a combined_gvcf VCF bgzipped file.

    GnarlyGenotyper performs "quick and dirty" joint genotyping on large cohorts,
    pre-called with HaplotypeCaller, and post-processed with ReblockGVCF.

    HaplotypeCaller must be used with `-ERC GVCF` or `-ERC BP_RESOLUTION` to add
    genotype likelihoods.

    ReblockGVCF must be run to add all the annotations necessary for VQSR:
    QUALapprox, VarDP, RAW_MQandDP.

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    job_name = 'Joint genotyping: GnarlyGenotyper'
    if interval_idx is not None:
        job_name += f' {interval_idx}/{number_of_intervals}'

    if output_vcf_path and utils.file_exists(output_vcf_path) and not overwrite:
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    # GnarlyGenotyper crashes with NullPointerException when using standard GATK docker
    j.image(utils.GNARLY_IMAGE)
    j.cpu(2)
    j.memory(f'32G')
    # 4G (fasta+fai+dict) + 1G per sample divided by the number of intervals
    j.storage(f'{4 + number_of_samples * 1 // number_of_intervals}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -e

    (while true; do df -h; pwd; free -m; sleep 300; done) &

    tar -xf {genomicsdb}

    df -h; pwd; free -m

    gatk --java-options -Xms8g \\
      GnarlyGenotyper \\
      -R {reference.base} \\
      -O {j.output_vcf['vcf.gz']} \\
      -D {dbsnp.base} \\
      --only-output-calls-starting-in-intervals \\
      --keep-all-sites \\
      -V gendb://workspace \\
      {f'-L {interval} ' if interval else ''} \\
      --create-output-variant-index"""
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_exccess_het_filter(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    excess_het_threshold: float = 54.69,
    interval: Optional[hb.ResourceGroup] = None,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Filter a large cohort callset on Excess Heterozygosity.

    The filter applies only to large callsets (`not is_small_callset`)

    Requires all samples to be unrelated.

    ExcessHet estimates the probability of the called samples exhibiting excess
    heterozygosity with respect to the null hypothesis that the samples are unrelated.
    The higher the score, the higher the chance that the variant is a technical artifact
    or that there is consanguinuity among the samples. In contrast to Inbreeding
    Coefficient, there is no minimal number of samples for this annotation.

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    job_name = 'Joint genotyping: ExcessHet filter'
    if output_vcf_path and utils.file_exists(output_vcf_path) and not overwrite:
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'32G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail

    # Captring stderr to avoid Batch pod from crashing with OOM from millions of
    # warning messages from VariantFiltration, e.g.:
    # > JexlEngine - ![0,9]: 'ExcessHet > 54.69;' undefined variable ExcessHet
    gatk --java-options -Xms3g \\
      VariantFiltration \\
      --filter-expression 'ExcessHet > {excess_het_threshold}' \\
      --filter-name ExcessHet \\
      {f'-L {interval} ' if interval else ''} \\
      -O {j.output_vcf['vcf.gz']} \\
      -V {input_vcf['vcf.gz']} \\
      2> {j.stderr}
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_make_sites_only_job(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Create sites-only VCF with only site-level annotations.
    Speeds up the analysis in the AS-VQSR modeling step.

    Returns: a Job object with a single output j.sites_only_vcf of type ResourceGroup
    """
    job_name = 'Joint genotyping: MakeSitesOnlyVcf'
    if output_vcf_path and utils.file_exists(output_vcf_path) and not overwrite:
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'32G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      MakeSitesOnlyVcf \\
      -I {input_vcf['vcf.gz']} \\
      -O {j.output_vcf['vcf.gz']}
      """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_final_gather_vcf_job(
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
    j.image(utils.GATK_IMAGE)
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

    (while true; do df -h; pwd free -m; sleep 300; done) &

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
    
    df -h; pwd; free -m
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
