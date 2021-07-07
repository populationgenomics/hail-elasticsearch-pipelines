#!/usr/bin/env python3

import os
import hailtop.batch as hb
from analysis_runner import dataproc

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'), bucket=os.getenv('HAIL_BUCKET')
)

batch = hb.Batch(name='Parquet conversion', backend=service_backend)

dataproc.hail_dataproc_job(
    batch,
    'annotated_mt_to_parquet.py',
    max_age='12h',
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name='convert',
)

batch.run()
