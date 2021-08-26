#!/usr/bin/env python3

"""Converts an annotated seqr dataset in Hail MatrixTable format to the Arrow IPC format"""

import click
import os
import hailtop.batch as hb
from analysis_runner import dataproc, output_path

PARQUET_TO_ARROW_SHARD_COUNT = 10
MT_TO_PARQUET_PY = 'mt_to_parquet.py'
PARQUET_TO_ARROW_PY = 'parquet_to_arrow.py'


@click.command()
@click.option(
    '--input', help='Input path for annotated Hail MatrixTable', required=True
)
def main(input):
    """Script entry point."""

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.getenv('HAIL_BUCKET'),
    )

    batch = hb.Batch(name='seqr table conversion', backend=service_backend)

    parquet_path = output_path('annotated_ht.parquet', 'tmp')
    arrow_path = output_path('annotated_ht.arrow')

    # mt_to_parquet_job = dataproc.hail_dataproc_job(
    #    batch,
    #    f'{MT_TO_PARQUET_PY} --input="{input}" --output="{parquet_path}"',
    #    max_age='4h',
    #    num_secondary_workers=10,
    #    packages=['click'],
    #    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    #    job_name='MT to Parquet',
    # )

    mt_to_parquet_job = batch.new_job('no-op')
    mt_to_parquet_job.image(os.getenv('DRIVER_IMAGE'))
    mt_to_parquet_job.command('echo hi')

    parquet_to_arrow_src = open(PARQUET_TO_ARROW_PY).read()
    for shard_index in range(PARQUET_TO_ARROW_SHARD_COUNT):
        parquet_to_arrow_job = batch.new_job(
            f'Parquet to Arrow {shard_index + 1}/{PARQUET_TO_ARROW_SHARD_COUNT}'
        )

        parquet_to_arrow_job.image(os.getenv('DRIVER_IMAGE'))
        parquet_to_arrow_job.depends_on(mt_to_parquet_job)

        parquet_to_arrow_job.command(
            'gcloud -q auth activate-service-account --key-file=/gsa-key/key.json'
        )
        parquet_to_arrow_job.command(
            f'cat <<EOF > {PARQUET_TO_ARROW_PY}\n{parquet_to_arrow_src}\nEOF'
        )
        parquet_to_arrow_job.command(
            f'python3 {PARQUET_TO_ARROW_PY} '
            f'--input={parquet_path} '
            f'--output={arrow_path} '
            f'--shard_index={shard_index} '
            f'--shard_count={PARQUET_TO_ARROW_SHARD_COUNT}'
        )

    # Don't wait, which avoids resubmissions if this job gets preempted.
    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter