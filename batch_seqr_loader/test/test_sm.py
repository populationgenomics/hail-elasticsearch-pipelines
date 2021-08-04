#!/usr/bin/env python3

"""
Simulate communication with the SM server during the workflow
"""

import os
import string
import random
from collections import defaultdict
from typing import Union, List

from sample_metadata import AnalysisUpdateModel
from sample_metadata.api import SampleApi, AnalysisApi
from sample_metadata.models.new_sample import NewSample
from sample_metadata.models.analysis_model import AnalysisModel


PROJ = os.environ.get('SM_DEV_DB_PROJECT', 'sm_dev')


def test_simulate_joint_calling_pipeline():
    """
    Simulates events of the joint-calling workflow
    """

    # test_run_id = 'AFYPXR'
    # sample_ids = 'CPG620, CPG638, CPG646'.split(', ')

    # Unique test run ID to avoid clashing with previous test run samples
    test_run_id = os.environ.get(
        'SM_DV_TEST_RUN_ID',
        ''.join(
            random.choice(string.ascii_uppercase + string.digits) for _ in range(6)
        ),
    )
    print(f'Test run ID: {test_run_id}')
    sapi = SampleApi()
    print('sapi = SampleApi()')
    s1 = NewSample(
        external_id=f'NA12878-from-fq-{test_run_id}',
        type='blood',
        meta={
            'reads': [
                [
                    'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878_L001_R1.fq',
                    'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878_L002_R1.fq',
                ],
                [
                    'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878_L001_R2.fq',
                    'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878_L002_R2.fq',
                ],
            ]
        },
    )
    print('s1 = NewSample(')
    sample_id = sapi.create_new_sample(PROJ, s1)
    print(f'Created {sample_id}')


if __name__ == '__main__':
    test_simulate_joint_calling_pipeline()
