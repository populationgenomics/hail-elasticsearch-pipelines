#!/usr/bin/env python3

"""
Submit test_simulate_sm_worklfow.py
"""

import os
from os.path import join, dirname, abspath
import hailtop.batch as hb


SM_CONTAINER = 'australia-southeast1-docker.pkg.dev/cpg-common/images/sm-api:latest'
BAZAM_CONTAINER = 'australia-southeast1-docker.pkg.dev/cpg-common/images/bazam:v2'
REF_BUCKET = 'gs://cpg-reference/hg38/v1'
TARGET_BUCKET = 'gs://cpg-seqr-test-tmp/hg38/v1'
REF_FASTA = join(REF_BUCKET, 'Homo_sapiens_assembly38.fasta')
TARGET_FASTA = join(TARGET_BUCKET, 'Homo_sapiens_assembly38.fasta')


def _make_test_simulate_sm_workflow_job(
    _b: hb.Batch,
):
    j = _b.new_job('Test simulate SM workflow')
    j.image(SM_CONTAINER)

    with open(join(dirname(abspath(__file__)), 'test_simulate_sm_worklfow.py')) as f:
        script = f.read()
    j.command(
        f"""set -e
cat <<EOT >> test_simulate_sm_worklfow.py
{script}
EOT
export SM_USE_SERVICE_ACCOUNT=true
export SM_DEV_DB_PROJECT=vladdev
export SM_ENVIRONMENT=PRODUCTION
python test_simulate_sm_worklfow.py
    """
    )


billing_project = os.getenv('HAIL_BILLING_PROJECT') or 'seqr'
hail_bucket = os.environ['HAIL_BUCKET']
print(
    f'Starting hail Batch with the project {billing_project}, ' f'bucket {hail_bucket}'
)
backend = hb.ServiceBackend(
    billing_project=billing_project,
    bucket=hail_bucket.replace('gs://', ''),
)
b = hb.Batch(backend=backend, name='test')
_make_test_simulate_sm_workflow_job(b)
b.run(open=True)
