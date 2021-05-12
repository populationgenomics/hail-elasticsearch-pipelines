#!/usr/bin/env python3

"""
Let this be the entrypoint / driver for loading data into SEQR for the CPG
See the README for more information. This is WIP.

    - 2021/04/16 Michael Franklin
"""

import logging
import os
import shutil
import subprocess
import tempfile
from os.path import join, basename
from typing import Optional, List
import pandas as pd
import click
import hailtop.batch as hb
from analysis_runner import dataproc
from hailtop.batch.job import Job

GATK_CONTAINER = 'broadinstitute/gatk:4.1.8.0'

BROAD_REF_BUCKET = 'gs://gcp-public-data--broad-references/hg38/v0'
REF_FASTA = join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.fasta')
DBSNP_VCF = join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf')
UNPADDED_INTERVALS = join(BROAD_REF_BUCKET, 'hg38.even.handcurated.20k.intervals')

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
    '--gvcf-bucket', 'gvcf_buckets', multiple=True, required=True,
)
@click.option(
    '--ped-file', 'ped_fpath', required=True,
)
@click.option(
    '--dataset', 'dataset_name', required=True,
)
@click.option(
    '--dest-mt-path', 'dest_mt_path', required=True,
)
@click.option(
    '--work-bucket', '--bucket', 'work_bucket', required=True,
)
@click.option('--keep-scratch', 'keep_scratch', is_flag=True)
@click.option('--dry-run', 'dry_run', is_flag=True)
@click.option(
    '--billing-project',
    'billing_project',
    type=str,
    default=os.getenv('HAIL_BILLING_PROJECT'),
)
@click.option('--genome-version', 'genome_version', default='GRCh38')
@click.option(
    '--disable-validation', 'disable_validation', is_flag=True
)
@click.option(
    '--sample-type', 'sample_type', type=click.Choice(['WGS', 'WES']), default='WGS'
)
@click.option(
    '--dataset-type', 'dataset_type', type=click.Choice(['VARIANTS', 'SV']), 
    default='VARIANTS'
)
@click.option(
    '--remap-tsv', 'remap_path',
    help='Path to a TSV file with two columns: s and seqr_id.'
)
@click.option(
    '--subset-tsv', 'subset_path',
    help='Path to a TSV file with one column of sample IDs: s.'
)
@click.option(
    '--vep-config', 'vep_config_json_path',
    help='Path of hail vep config .json file'
)
@click.option(
    '--vep-block-size', 'vep_block_size'
)
def main(
    gvcf_buckets: List[str],
    ped_fpath: str,
    dataset_name: str,
    dest_mt_path: str,
    work_bucket: str,
    keep_scratch: bool,
    dry_run: bool,
    billing_project: Optional[str],
    genome_version: str,
    disable_validation: bool,
    dataset_type: str,
    sample_type: str,
    remap_path: str = None,
    subset_path: str = None,
    vep_config_json_path: Optional[str] = None,
    vep_block_size: Optional[int] = None,
):
    """
    Entry point for a batch workflow
    """
    if not billing_project:
        raise click.BadParameter(
            '--billing-project has to be specified'
        )

    hail_bucket = os.environ.get('HAIL_BUCKET') or join(work_bucket, 'hail')
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=hail_bucket.replace('gs://', ''),
    )
    b = hb.Batch('Seqr loader', backend=backend)

    samples_df = find_inputs(gvcf_buckets, ped_fpath=ped_fpath)
    gvcfs = [
        b.read_input_group(**{'g.vcf.gz': gvcf, 'g.vcf.gz.tbi': gvcf + '.tbi'})
        for gvcf in list(samples_df.gvcf)
    ]

    reference = b.read_input_group(
        base=REF_FASTA,
        fai=REF_FASTA + '.fai',
        dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '') + '.dict',
    )
    dbsnp = b.read_input_group(base=DBSNP_VCF, idx=DBSNP_VCF + '.idx')

    scatter_count_scale_factor = 0.15
    scatter_count = int(round(scatter_count_scale_factor * len(gvcfs)))
    scatter_count = max(scatter_count, 2)
    split_intervals_job = add_split_intervals_job(
        b,
        UNPADDED_INTERVALS,
        scatter_count,
        reference,
    )
    intervals = split_intervals_job.intervals

    genotype_vcf_jobs = []
    sample_map_fpath = join(work_bucket, 'work', 'sample_name.csv')
    samples_df[['s', 'gvcf']].to_csv(sample_map_fpath,
                                     sep='\t', header=False, index=False)
    for idx in range(scatter_count):
        import_gvcfs_job = add_import_gvcfs_job(
            b=b,
            sample_name_map=b.read_input(sample_map_fpath),
            interval=intervals[f'interval_{idx}'],
        )
        genotype_vcf_job = add_gatk_genotype_gvcf_job(
            b,
            workspace_tar=import_gvcfs_job.genomicsdb_tar,
            interval=intervals[f'interval_{idx}'],
            reference=reference,
            dbsnp=dbsnp,
            disk_size=100,
        )
        genotype_vcf_jobs.append(genotype_vcf_job)

    final_gathered_vcf_job = add_final_gather_vcf_step(
        b,
        input_vcfs=[j.output_vcf for j in genotype_vcf_jobs],
        disk_size=200,
        output_vcf_path=join(work_bucket, f'{dataset_name}.vcf'),
    )

    j = dataproc.hail_dataproc_job(
        b,
        f'batch_seqr_loader/seqr_load.py '
        f'--source-path {final_gathered_vcf_job.output_vcf} '
        f'--dest-mt-path {dest_mt_path} --',
        max_age='8h',
        packages=DATAPROC_PACKAGES,
        num_secondary_workers=2,
        job_name='seqr_load.py',
        vep='GRCh38',
        depends_on=[final_gathered_vcf_job],
    )

    b.run(dry_run=dry_run, delete_scratch_on_exit=not keep_scratch)


def find_inputs(
    input_buckets: List[str],
    ped_fpath: str,
) -> pd.DataFrame:  # pylint disable=too-many-branches
    """
    Read the inputs assuming a standard CPG storage structure.
    :param input_buckets: buckets to find GVCFs and CSV metadata files.
    :param ped_fpath: pedigree file
    :return: a dataframe with the pedigree information and paths to gvcfs
    """
    gvcf_paths: List[str] = []
    for ib in input_buckets:
        cmd = f'gsutil ls \'{join(ib, "*.g.vcf.gz")}\''
        gvcf_paths.extend(
            line.strip()
            for line in subprocess.check_output(cmd, shell=True).decode().split()
        )

    local_tmp_dir = tempfile.mkdtemp()
    local_ped_fpath = join(local_tmp_dir, basename(ped_fpath))
    subprocess.run(
        f'gsutil cp {ped_fpath} {local_ped_fpath}', check=False, shell=True
    )
    df = pd.read_csv(local_ped_fpath, delimiter='\t')
    shutil.rmtree(local_tmp_dir)
    df['gvcf'] = ''
    df = df.set_index('Individual.ID', drop=False)
    df = df.rename(columns={'Individual.ID': 's'})

    sample_names = list(df['s'])

    # Checking 1-to-1 match of sample names to GVCFs
    for sn in sample_names:
        matching_gvcfs = [gp for gp in gvcf_paths if sn in gp]
        if len(matching_gvcfs) > 1:
            logging.warning(
                f'Multiple GVCFs found for the sample {sn}:' f'{matching_gvcfs}'
            )
        elif len(matching_gvcfs) == 0:
            logging.warning(f'No GVCFs found for the sample {sn}')

    # Checking 1-to-1 match of GVCFs to sample names, and filling a dict
    for gp in gvcf_paths:
        matching_sn = [sn for sn in sample_names if sn in gp]
        if len(matching_sn) > 1:
            logging.warning(
                f'Multiple samples found for the GVCF {gp}:' f'{matching_sn}'
            )
        elif len(matching_sn) == 0:
            logging.warning(f'No samples found for the GVCF {gp}')
        else:
            df.loc[matching_sn[0], ['gvcf']] = gp
    df = df[df.gvcf.notnull()]
    return df


def add_split_intervals_job(
    b: hb.Batch,
    interval_list: hb.ResourceFile,
    scatter_count: int,
    ref_fasta: hb.ResourceGroup,
    disk_size: int = 30,
) -> Job:
    """
    Split genome into intervals to parallelise GnarlyGenotyper.

    Returns: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job('SplitIntervals')
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

    gatk --java-options -Xms{mem_gb - 1}g SplitIntervals \\
      -L {interval_list} \\
      -O {j.intervals} \\
      -scatter {scatter_count} \\
      -R {ref_fasta.base} \\
      -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
      """
    )
    return j


def add_import_gvcfs_job(
    b: hb.Batch,
    sample_name_map: hb.ResourceFile,
    interval: hb.ResourceFile,
    disk_size: int = 30,
):
    j = b.new_job('ImportGVCFs')
    j.image(GATK_CONTAINER)
    mem_gb = 16
    j.memory(f'{mem_gb}G')
    j.storage(f'{disk_size}G')
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

    gatk --java-options -Xms{mem_gb - 1}g \
      GenomicsDBImport \
      --genomicsdb-workspace-path /genomicsdb_workspace \
      --batch-size 50 \
      -L {interval} \
      --sample-name-map {sample_name_map} \
      --reader-threads 5 \
      --merge-input-intervals \
      --consolidate

    tar -cf {j.genomicsdb_tar} /genomicsdb_workspace"""
    )
    return j


def add_gatk_genotype_gvcf_job(
    b: hb.Batch,
    workspace_tar: hb.ResourceGroup,
    interval: hb.ResourceFile,
    reference: hb.ResourceGroup,
    dbsnp: hb.ResourceGroup,
    disk_size: int = 100,
) -> Job:
    j = b.new_job('GenotypeGVCFs')
    j.image(GATK_CONTAINER)
    j.memory(f'32G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
        }
    )

    j.command(
        f"""set -e
        
    tar -xf {workspace_tar}
    WORKSPACE=$(basename {workspace_tar} .tar)

    gatk \\
      --java-options -Xms8g \\
      --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \\
      GenotypeGVCFs \\
      -R {reference.base} \\
      -O {j.output_vcf} \\
      -D {dbsnp.base} \\
      -G StandardAnnotation -G AS_StandardAnnotation \\
      --only-output-calls-starting-in-intervals \\
      -V gendb://$WORKSPACE \\
      -L {interval} \\
      --merge-input-intervals
    """
    )
    return j


def add_final_gather_vcf_step(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceGroup],
    disk_size: int,
    output_vcf_path: str = None,
) -> Job:
    """
    Combines recalibrated VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    j = b.new_job('FinalGatherVcf')
    j.image(GATK_CONTAINER)
    j.memory(f'8G')
    j.storage(f'{disk_size}G')
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
    gatk --java-options -Xms6g \\
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


if __name__ == '__main__':
    main()  # pylint: disable=E1120
