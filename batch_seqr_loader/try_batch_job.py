#!/usr/bin/env python3

import hailtop.batch as hb
from os.path import join
import sys, os

hail_bucket = os.getenv('HAIL_BUCKET')
billing_project = os.getenv('HAIL_BILLING_PROJECT')
hail_token = os.getenv('HAIL_TOKEN')

backend = hb.ServiceBackend(
    billing_project=billing_project,
    bucket=hail_bucket.replace('gs://', ''),
    token=hail_token,
)
batch = hb.Batch(backend=backend, name='test-batch')

j = batch.new_job(name='test-job')
j.command(
    f'''
echo "test" > {j.ofile}
'''
)
j.cpu(2)
# j.image('australia-southeast1-docker.pkg.dev/cpg-common/images/aspera:v1')
batch.run(wait=False)
