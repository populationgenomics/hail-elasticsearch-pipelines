#!/usr/bin/env python3

"""
Let this be the entrypoint / driver for loading data into SEQR for the CPG
See the README for more information. This is WIP.

    - 2021/04/16 Michael Franklin and Vlad Savelyev
"""
import json
import logging
import math
import os
import re
import shutil
import subprocess
import tempfile
import hashlib
from collections import defaultdict
from os.path import join, basename, dirname, abspath
from typing import Optional, List, Tuple, Dict, Set, Iterable
import pandas as pd
import click
import hailtop.batch as hb
from analysis_runner import dataproc
from hailtop.batch.job import Job
from google.cloud import storage

GATK_VERSION = '4.2.0.0'
GATK_CONTAINER = (
    f'australia-southeast1-docker.pkg.dev/cpg-common/images/gatk:{GATK_VERSION}'
)
PICARD_CONTAINER = f'us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8'
BAZAM_CONTAINER = f'australia-southeast1-docker.pkg.dev/fewgenomes/images/bazam:v2'
SOMALIER_CONTAINER = 'brentp/somalier:latest'
PEDDY_CONTAINER = 'quay.io/biocontainers/peddy:0.4.8--pyh5e36f6f_0'

# Fixed scatter count, because we are storing a genomics DB per each interval
NUMBER_OF_INTERVALS = 10

REF_BUCKET = 'gs://cpg-reference/hg38/v0'
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

logger = logging.getLogger('seqr-loader')
logger.setLevel('INFO')


@click.command()
@click.option(
    '--gvcf-bucket',
    'gvcf_buckets',
    multiple=True,
    help='Bucket path with input GVCF files along with tbi indices',
)
@click.option(
    '--bam-bucket',
    '--cram-bucket',
    'bam_buckets',
    multiple=True,
    help='Bucket path with input BAM or CRAM files',
)
@click.option(
    '--bam-to-realign-bucket',
    '--cram-to-realign-bucket',
    'bam_to_realign_buckets',
    multiple=True,
    help='Bucket path with input BAM or CRAM files that need re-alignment',
)
@click.option(
    '--data-bucket',
    'data_bucket',
    required=True,
    help='Base bucket path to store genomics DBs, fingerprints, and matrix tables',
)
@click.option(
    '--ped-file', 'ped_fpath', help='Pedigree file to verify the relatedness and sex'
)
@click.option(
    '--work-bucket',
    '--bucket',
    'work_bucket',
    required=True,
    help='Base bucket path to store temporary and intermediate files',
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
@click.option('--vep-block-size', 'vep_block_size')
def main(
    gvcf_buckets: List[str],
    bam_buckets: List[str],
    bam_to_realign_buckets: List[str],
    data_bucket: str,
    ped_fpath: str,
    work_bucket: str,
    keep_scratch: bool,
    overwrite: bool,
    dry_run: bool,
    disable_validation: bool,  # pylint: disable=unused-argument
    remap_path: str,
    subset_path: str,
    make_checkpoints: bool,
    vep_block_size: Optional[int] = None,  # pylint: disable=unused-argument
):
    """
    Entry point for a batch workflow
    """
    if not (gvcf_buckets or bam_buckets or bam_to_realign_buckets):
        raise click.BadParameter(
            'Specify at least one of the input parameters '
            '(can be multiple and/or repeated): '
            '--gvcf-bucket, --bam-bucket, --bam-to-realign-bucket'
        )

    billing_project = os.getenv('HAIL_BILLING_PROJECT') or 'seqr'

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
    b = hb.Batch('Somalier check', backend=backend)

    local_tmp_dir = tempfile.mkdtemp()

    samples_df = find_inputs(
        gvcf_buckets,
        bam_buckets,
        bam_to_realign_buckets,
        local_tmp_dir=local_tmp_dir,
        ped_fpath=ped_fpath,
    )

    reference = b.read_input_group(
        base=REF_FASTA,
        fai=REF_FASTA + '.fai',
        dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
        + '.dict',
    )

    ped_check_j = _pedigree_checks(
        b=b,
        samples_df=samples_df,
        reference=reference,
        sites=b.read_input(SOMALIER_SITES),
        ped_file=b.read_input(ped_fpath),
        work_bucket=work_bucket,
        overwrite=overwrite,
        fingerprints_bucket=join(data_bucket, 'fingerprints'),
    )

    # realign_bam_jobs = _make_realign_bam_jobs(
    #     b=b,
    #     samples_df=samples_df,
    #     reference=reference,
    #     work_bucket=work_bucket,
    #     overwrite=overwrite,
    # )

    intervals_j = _add_split_intervals_job(
        b,
        UNPADDED_INTERVALS,
        NUMBER_OF_INTERVALS,
        REF_FASTA,
    )
    intervals_j.depends_on(ped_check_j)

    genotype_jobs, samples_df = _make_genotype_jobs(
        b=b,
        intervals=intervals_j.intervals,
        samples_df=samples_df,
        reference=reference,
        work_bucket=work_bucket,
        overwrite=overwrite,
    )

    joint_genotype_job, gathered_vcf_path = _make_joint_genotype_jobs(
        b=b,
        intervals=intervals_j.intervals,
        genomicsdb_bucket=join(data_bucket, 'genomicsdb'),
        samples_df=samples_df,
        reference=reference,
        dbsnp=b.read_input_group(base=DBSNP_VCF, idx=DBSNP_VCF + '.idx'),
        work_bucket=work_bucket,
        local_tmp_dir=local_tmp_dir,
        overwrite=overwrite,
        depends_on=genotype_jobs,
    )

    annotated_mt_path = join(dirname(gathered_vcf_path), 'annotated.mt')
    dataproc.hail_dataproc_job(
        b,
        f'batch_seqr_loader/seqr_load.py '
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
        num_secondary_workers=2,
        job_name='seqr_load.py',
        vep='GRCh38',
        depends_on=[joint_genotype_job],
    )

    b.run(dry_run=dry_run, delete_scratch_on_exit=not keep_scratch)
    shutil.rmtree(local_tmp_dir)


# def _make_realign_bam_jobs(
#     b: hb.Batch,
#     samples_df: pd.DataFrame,
#     reference: hb.ResourceGroup,
#     work_bucket: str,
#     overwrite: bool,
# ) -> List[Job]:
#     return []


def _pedigree_checks(
    b: hb.Batch,
    samples_df: pd.DataFrame,
    reference: hb.ResourceGroup,
    sites: hb.ResourceFile,
    ped_file: hb.ResourceFile,
    work_bucket: str,  # pylint: disable=unused-argument
    overwrite: bool,  # pylint: disable=unused-argument
    fingerprints_bucket: str,
    depends_on: Optional[List[Job]] = None,
) -> Job:
    """
    Add somalier and peddy based jobs that infer relatedness
    and sex, and compare that to the provided ped file, and fail
    if anything don't match
    """

    extract_jobs = []
    fp_file_by_sample = dict()
    for sn, input_path, input_index in zip(
        samples_df['s'], samples_df['file'], samples_df['index']
    ):
        fp_file_by_sample[sn] = join(fingerprints_bucket, f'{sn}.somalier')
        if can_reuse(fp_file_by_sample[sn], overwrite):
            extract_jobs.append(b.new_job(f'Somalier Extract, {sn} [reuse]'))
        else:
            j = b.new_job(f'Somalier extract, {sn}')
            j.image(SOMALIER_CONTAINER)
            j.memory(f'8G')
            j.storage(f'10G')
            if depends_on:
                j.depends_on(*depends_on)

            input_file = b.read_input_group(
                base=input_path,
                index=input_index,
            )

            j.command(
                f"""set -e
        
                somalier extract \\
                -d extracted/ \\
                --sites {sites} \\
                -f {reference.base} \\
                {input_file['base']}
                
                mv extracted/{sn}.somalier {j.output_file}
                """
            )
            b.write_output(j.output_file, fp_file_by_sample[sn])
            extract_jobs.append(j)

    relate_j = b.new_job(f'Somalier relate')
    relate_j.image(SOMALIER_CONTAINER)
    relate_j.memory(f'8G')
    relate_j.storage(f'10G')
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

    sample_hash = hash_sample_names(samples_df['s'])
    prefix = join(fingerprints_bucket, sample_hash, 'somalier')
    somalier_html_path = f'{prefix}.html'
    somalier_samples_path = f'{prefix}.samples.tsv'
    somalier_pairs_path = f'{prefix}.pairs.tsv'
    b.write_output(relate_j.output_html, somalier_html_path)
    b.write_output(relate_j.output_samples, somalier_samples_path)
    b.write_output(relate_j.output_pairs, somalier_pairs_path)

    check_j = b.new_job(f'Check relatedness and sex')
    check_j.image(PEDDY_CONTAINER)
    check_j.memory(f'8G')
    check_j.storage(f'10G')
    with open(join(dirname(abspath(__file__)), 'check_pedigree.py')) as f:
        script = f.read()
    check_j.command(
        f"""set -e
cat <<EOT >> check_pedigree.py
{script}
EOT
python check_pedigree.py --somalier-prefix {prefix}
    """
    )

    check_j.depends_on(relate_j)
    return check_j


def _make_genotype_jobs(
    b: hb.Batch,
    intervals: hb.ResourceGroup,
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
    merge_gvcf_jobs = []
    bams_df = samples_df[samples_df['type'] == 'bam']
    for sn, bam_fpath in zip(bams_df['s'], bams_df['file']):
        output_gvcf_path = join(work_bucket, f'{sn}.g.vcf.gz')
        if can_reuse(output_gvcf_path, overwrite):
            merge_gvcf_jobs.append(b.new_job(f'HaplotypeCaller, {sn} [reuse]'))
        else:
            # if bam_fpath.endswith('.cram'):
            #     index_path = os.path.splitext(bam_fpath)[0] + '.crai'
            # else:
            #     index_path = os.path.splitext(bam_fpath)[0] + '.bai'
            # input_bam = b.read_input_group(
            #     bam=bam_fpath,
            #     bai=index_path,
            # )
            haplotype_caller_jobs = []
            for idx in range(NUMBER_OF_INTERVALS):
                haplotype_caller_jobs.append(
                    _add_haplotype_caller_job(
                        b,
                        bam_fpath,
                        interval=intervals[f'interval_{idx}'],
                        reference=reference,
                        sample_name=sn,
                        interval_idx=idx,
                        number_of_intervals=NUMBER_OF_INTERVALS,
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
    intervals: hb.ResourceGroup,
    genomicsdb_bucket: str,
    samples_df: pd.DataFrame,
    reference: hb.ResourceGroup,
    dbsnp: hb.ResourceGroup,
    work_bucket: str,
    local_tmp_dir: str,
    overwrite: bool,
    depends_on: List[Job] = None,
) -> Tuple[Job, str]:
    """
    Assumes all samples have a 'file' of 'type'='gvcf' in `samples_df`.
    Adds samples to the GenomicsDB and runs joint genotyping on them.
    Outputs a multi-sample VCF under `output_vcf_path`.
    """

    genotype_vcf_jobs = []
    genotyped_vcfs = []

    samples_to_be_in_db_per_interval = dict()

    for idx in range(NUMBER_OF_INTERVALS):
        genomicsdb_gcs_path = join(
            genomicsdb_bucket, f'interval_{idx}_outof_{NUMBER_OF_INTERVALS}.tar'
        )

        import_gvcfs_job, samples_to_be_in_db = _add_import_gvcfs_job(
            b=b,
            genomicsdb_gcs_path=genomicsdb_gcs_path,
            samples_df=samples_df,
            work_bucket=work_bucket,
            local_tmp_dir=local_tmp_dir,
            interval=intervals[f'interval_{idx}'],
            interval_idx=idx,
            number_of_intervals=NUMBER_OF_INTERVALS,
            depends_on=depends_on,
        )

        samples_to_be_in_db_per_interval[idx] = samples_to_be_in_db
        samples_hash = hash_sample_names(samples_to_be_in_db)
        output_vcf_path = join(
            work_bucket, 'jointly-called', samples_hash, f'interval_{idx}.vcf.gz'
        )
        genotype_vcf_job = _add_gatk_genotype_gvcf_job(
            b,
            genomicsdb=b.read_input(genomicsdb_gcs_path),
            interval=intervals[f'interval_{idx}'],
            reference=reference,
            dbsnp=dbsnp,
            overwrite=overwrite,
            interval_idx=idx,
            number_of_intervals=NUMBER_OF_INTERVALS,
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


def find_inputs(
    gvcf_buckets: List[str],
    bam_buckets: List[str],
    bam_to_realign_buckets: List[str],
    local_tmp_dir: str,
    ped_fpath: Optional[str] = None,
) -> pd.DataFrame:  # pylint disable=too-many-branches
    """
    Read the inputs assuming a standard CPG storage structure.
    :param gvcf_buckets: buckets to find GVCF files
    :param bam_buckets: buckets to find BAM files
        (will be passed to HaplotypeCaller to produce GVCFs)
    :param bam_to_realign_buckets: buckets to find BAM files
        (will be re-aligned with BWA before passing to HaplotypeCaller)
    :param ped_fpath: pedigree file
    :param local_tmp_dir: temporary local directory
    :return: a dataframe with the pedigree information and paths to gvcfs
    """
    input_buckets_by_type = dict(
        gvcf=gvcf_buckets,
        bam=bam_buckets,
        bam_to_realign=bam_to_realign_buckets,
    )
    found_files_by_type: Dict[str, List[str]] = {t: [] for t in input_buckets_by_type}
    found_indices_by_type: Dict[str, List[str]] = {t: [] for t in input_buckets_by_type}

    for input_type, buckets in input_buckets_by_type.items():
        patterns = ['*.g.vcf.gz'] if input_type == 'gvcf' else ['*.bam', '*.cram']
        for ib in buckets:
            for pattern in patterns:
                cmd = f"gsutil ls '{join(ib, pattern)}'"
                try:
                    found_files_by_type[input_type].extend(
                        line.strip()
                        for line in subprocess.check_output(cmd, shell=True)
                        .decode()
                        .split()
                    )
                except subprocess.CalledProcessError:
                    pass

    for input_type, fpaths in found_files_by_type.items():
        for fp in fpaths:
            index = fp + '.tbi'
            if fp.endswith('.g.vcf.gz'):
                if not file_exists(index):
                    logger.critical(f'Not found TBI index for file {fp}')
            elif fp.endswith('.bam'):
                index = fp + '.bai'
                if not file_exists(index):
                    index = re.sub('.bam$', '.bai', fp)
                    if not file_exists(index):
                        logger.critical(f'Not found BAI index for file {fp}')
            elif fp.endswith('.cram'):
                index = fp + '.crai'
                if not file_exists(index):
                    index = re.sub('.cram$', '.crai', fp)
                    if not file_exists(index):
                        logger.critical(f'Not found CRAI index for file {fp}')
            else:
                logger.critical(f'Unrecognised input file extention {fp}')
            found_indices_by_type[input_type].append(index)

    def _get_base_name(file_path):
        return re.sub('(.bam|.cram|.g.vcf.gz)$', '', os.path.basename(file_path))

    if ped_fpath:
        local_ped_fpath = join(local_tmp_dir, basename(ped_fpath))
        subprocess.run(
            f'gsutil cat {ped_fpath} | cut -f1-6 | grep -v ^Family.ID > {local_ped_fpath}',
            check=False,
            shell=True,
        )
        df = pd.read_csv(
            local_ped_fpath,
            delimiter='\t',
            names=[
                'Family.ID',
                'Individual.ID',
                'Paternal.ID',
                'Maternal.ID',
                'Gender',
                'Phenotype',
            ],
        )
        df = df.rename(columns={'Individual.ID': 's'})
        df = df.set_index('s', drop=False)
        ped_snames = list(df['s'])

        # Checking that file base names have a 1-to-1 match with the samples in PED
        # First, checking the match of PED sample names to input files
        all_input_snames = []
        for fpaths in found_files_by_type.values():
            all_input_snames.extend([_get_base_name(fp) for fp in fpaths])

        for ped_sname in ped_snames:
            matching_files = [
                input_sname
                for input_sname in all_input_snames
                if ped_sname == input_sname
            ]
            if len(matching_files) > 1:
                logging.warning(
                    f'Multiple input files found for the sample {ped_sname}:'
                    f'{matching_files}'
                )
            elif len(matching_files) == 0:
                logging.warning(f'No files found for the sample {ped_sname}')

        # Second, checking the match of input files to PED sample names, and filling a dict
        for input_type in found_files_by_type:
            for fp, index in zip(
                found_files_by_type[input_type], found_indices_by_type[input_type]
            ):
                input_sname = _get_base_name(fp)
                matching_sn = [sn for sn in ped_snames if sn == input_sname]
                if len(matching_sn) > 1:
                    logging.warning(
                        f'Multiple samples found for the input {input_sname}:'
                        f'{matching_sn}'
                    )
                elif len(matching_sn) == 0:
                    logging.warning(f'No samples found for the input {input_sname}')
                else:
                    df.loc[matching_sn[0], 'type'] = input_type
                    df.loc[matching_sn[0], 'file'] = fp
                    df.loc[matching_sn[0], 'index'] = index

        df = df[df.file.notnull()]

    else:
        data: Dict[str, List] = dict(s=[], file=[], type=[])
        for input_type in found_files_by_type:
            for fp, index in zip(
                found_files_by_type[input_type], found_indices_by_type[input_type]
            ):
                data['s'].append(_get_base_name(fp))
                data['file'].append(fp)
                data['index'].append(index)
                data['type'].append(input_type)

        df = pd.DataFrame(data=data).set_index('s', drop=False)

    return df


def _add_split_intervals_job(
    b: hb.Batch,
    interval_list: str,
    scatter_count: int,
    ref_fasta: str,
    disk_size: int = 30,
) -> Job:
    """
    Split genome into intervals to parallelise GnarlyGenotyper.

    Returns: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job(f'Make {scatter_count} intervals')
    j.image(GATK_CONTAINER)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'{disk_size}G')
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
    gatk --java-options -Xms{mem_gb - 1}g SplitIntervals \\
      -L {interval_list} \\
      -O {j.intervals} \\
      -scatter {scatter_count} \\
      -R {ref_fasta} \\
      -mode INTERVAL_SUBDIVISION
      """
    )
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

    # We need disk to localize the sharded input and output due to the scatter for
    # HaplotypeCaller. If we take the number we are scattering by and reduce by 20,
    # we will have enough disk space to account for the fact that the data is quite
    # uneven across the shards.
    scatter_divisor = max((number_of_intervals - 20), 1)
    input_and_output_size = 40
    reference_data_size = 20
    disk_size = math.ceil(input_and_output_size / scatter_divisor) + reference_data_size

    job_name = 'HaplotypeCaller'
    if interval_idx is not None:
        job_name += f', {sample_name} {interval_idx}/{number_of_intervals}'

    j = b.new_job('HaplotypeCaller')
    j.image(GATK_CONTAINER)
    j.cpu(2)
    mem_gb = 64
    j.memory(f'{mem_gb}G')
    j.storage(f'{disk_size}G')
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
    gatk --java-options "-Xms{mem_gb - 1}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      -R {reference.base} \
      -I {bam_fpath} \
      -L {interval} \
      -O {j.output_gvcf['g.vcf.gz']} \
      -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      -ERC GVCF \
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
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'25G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )

    input_cmd = ' '.join(f'INPUT={g["g.vcf.gz"]}' for g in gvcfs)

    j.command(
        f"""set -e
    java -Xms{mem_gb - 1}g -jar /usr/picard/picard.jar \
      MergeVcfs {input_cmd} OUTPUT={j.output_gvcf['g.vcf.gz']}
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
        subprocess.run(
            f'gsutil cat {genomicsdb_gcs_path} | tar -O --extract workspace/callset.json > {genomicsdb_metadata}',
            check=False,
            shell=True,
        )

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
    mem_gb = 16
    j.memory(f'{mem_gb}G')
    j.storage(f'100G')
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
    
    gatk --java-options -Xms{mem_gb - 1}g \
      GenomicsDBImport \
      {genomicsdb_param} \
      --batch-size 50 \
      -L {interval} \
      --sample-name-map {sample_name_map} \
      --reader-threads 5 \
      --merge-input-intervals \
      --consolidate

    tar -cf {j.output['tar']} workspace
    """
    )
    b.write_output(j.output, genomicsdb_gcs_path.replace('.tar', ''))
    return j, samples_will_be_in_db


def _add_gatk_genotype_gvcf_job(
    b: hb.Batch,
    genomicsdb: hb.ResourceFile,
    interval: hb.ResourceFile,
    reference: hb.ResourceGroup,
    dbsnp: hb.ResourceGroup,
    overwrite: bool,
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
    mem_gb = 16
    j.memory(f'{mem_gb}G')
    j.storage(f'200G')
    j.declare_resource_group(
        output_vcf={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
        }
    )

    j.command(
        f"""set -e
        
    tar -xf {genomicsdb}

    gatk --java-options -Xms{mem_gb - 1}g \\
      GenotypeGVCFs \\
      -R {reference.base} \\
      -O {j.output_vcf['vcf.gz']} \\
      -D {dbsnp.base} \\
      --only-output-calls-starting-in-intervals \\
      -V gendb://workspace \\
      -L {interval} \\
      --merge-input-intervals
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
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'100G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in 
    # our invocation. This argument disables expensive checks that the file headers 
    # contain the same set of genotyped samples and that files are in order 
    # by position of first record.
    gatk --java-options -Xms{mem_gb - 1}g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      {input_cmdl} \\
      --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}"""
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j


def file_exists(path: str) -> bool:
    """
    Check if the object exists, where the object can be:
        * local file
        * local directory
        * Google Storage object
        * Google Storage URL representing a *.mt or *.ht Hail data,
          in which case it will check for the existence of a
          *.mt/_SUCCESS or *.ht/_SUCCESS file.
    :param path: path to the file/directory/object/mt/ht
    :return: True if the object exists
    """
    if path.startswith('gs://'):
        bucket = path.replace('gs://', '').split('/')[0]
        path = path.replace('gs://', '').split('/', maxsplit=1)[1]
        path = path.rstrip('/')  # '.mt/' -> '.mt'
        if any(path.endswith(f'.{suf}') for suf in ['mt', 'ht']):
            path = os.path.join(path, '_SUCCESS')
        gs = storage.Client()
        return gs.get_bucket(bucket).get_blob(path)
    return os.path.exists(path)


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


def hash_sample_names(samples_names: Iterable[str]) -> str:
    """
    Return a unique hash string from a from a set of strings
    :param samples_names: set of strings
    :return: a string hash
    """
    return hashlib.sha224(' '.join(sorted(samples_names)).encode()).hexdigest()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
