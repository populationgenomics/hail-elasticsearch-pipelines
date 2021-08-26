#!/usr/bin/env python3

"""Converts an annotated seqr dataset in Hail MatrixTable format to the Arrow IPC format"""

import click
import os
import hailtop.batch as hb
from analysis_runner import dataproc


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

    dataproc.hail_dataproc_job(
        batch,
        f'convert_to_parquet.py --input="{input}"',
        max_age='4h',
        num_secondary_workers=10,
        packages=['click'],
        init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
        job_name='convert to Parquet',
    )

    # Don't wait, which avoids resubmissions if this job gets preempted.
    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter