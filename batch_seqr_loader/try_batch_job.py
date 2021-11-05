#!/usr/bin/env python3

import hailtop.batch as hb
from os.path import join
import sys, os

backend = hb.ServiceBackend()
batch = hb.Batch(backend=backend, name='test-batch')

j = batch.new_job(name='test-job-from-analysis-runner')
j.command(
    f'''
echo "test" > {j.ofile}
'''
)
j.cpu(2)
# j.image('australia-southeast1-docker.pkg.dev/cpg-common/images/aspera:v1')
batch.run(wait=False)
