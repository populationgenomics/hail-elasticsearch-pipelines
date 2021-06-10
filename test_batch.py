#!/usr/bin/env python3


import os
import hailtop.batch as hb
from analysis_runner import dataproc


billing_project = os.getenv('HAIL_BILLING_PROJECT')
hail_bucket = os.environ.get('HAIL_BUCKET')
backend = hb.ServiceBackend(
    billing_project=billing_project,
    bucket=hail_bucket.replace('gs://', ''),
)
b = hb.Batch('Test VEP', backend=backend)

hail_dataproc_job(
    b,
    'hail_with_vep.py ',
    max_age='1h',
    num_secondary_workers=1,
    job_name='hail_with_vep.py',
    vep='GRCh38'
)

b.run()


