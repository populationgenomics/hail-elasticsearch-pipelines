#!/usr/bin/env python3

"""
Let this be the entrypoint / driver for loading data into SEQR for the CPG
See the README for more information. This is WIP.

    - 2021/04/16 Michael Franklin
"""

import logging
import os
from os.path import join
from typing import Optional, List

import click
import hailtop.batch as hb
from analysis_runner import dataproc

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
    '--source-path', 'source_paths', multiple=True, required=True,
)
@click.option(
    '--dest-path', 'dest_path', required=True,
)
@click.option(
    '--work-bucket', 'work_bucket', required=True,
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
    source_paths: List[str],
    dest_path: str,
    work_bucket: str,
    keep_scratch: bool,
    dry_run: bool,
    billing_project: str,
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
    Args:
        source_paths: Path or list of paths of VCFs to be loaded.
        dest_path: Path to write the matrix table.
        genome_version: Reference Genome Version (37 or 38)
        disable_validation: Disable checking whether the dataset matches the specified genome version and WGS vs. WES sample type.
        dataset_type: VARIANTS or SV.
        remap_path: Path to a tsv file with two columns: s and seqr_id.
        subset_path: Path to a tsv file with one column of sample IDs: s.
        vep_config_json_path: Path of hail vep config .json file
        vep_block_size: Block size to parallelize VEP annotation

    Returns:
    """
    if not dry_run:
        if not billing_project:
            raise click.BadParameter(
                '--billing-project has to be specified (unless --dry-run is set)'
            )

    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=join(work_bucket, 'hail').replace('gs://', ''),
    )
    b = hb.Batch('Seqr loader', backend=backend)
    sp_cmdl = ' '.join([f'--source-path {sp}' for sp in source_paths])
    j = dataproc.hail_dataproc_job(
        b,
        f'batch_seqr_loader/seqr_load.py {sp_cmdl} --dest-path {dest_path} --',
        max_age='8h',
        packages=DATAPROC_PACKAGES,
        num_secondary_workers=2,
        job_name='Load to Seqr',
        vep='GRCh38',
    )
    j.always_run()

    b.run(dry_run=dry_run, delete_scratch_on_exit=not keep_scratch)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
