#!/usr/bin/env python3

"""
Simulate communication with the SM server during the workflow
"""

import os
from datetime import datetime
from sample_metadata.api import SampleApi
from sample_metadata.models.new_sample import NewSample


PROJ = os.environ.get('SM_DEV_DB_PROJECT', 'sm_dev')


sapi = SampleApi()
s1 = NewSample(
    external_id='210303_A00692_0190_ML211637',
    type='blood',
    meta={
        'reads': [
            [
                'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L001_R1.fastq.gz',
                'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L002_R1.fastq.gz',
            ],
            [
                'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L001_R2.fastq.gz',
                'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L002_R2.fastq.gz',
            ],
        ]
    },
)
s2 = NewSample(
    external_id='20W000094-FAM000347',
    type='blood',
    meta={'reads': 'gs://cpg-acute-care-test-upload/20W000094-FAM000347.cram'},
)
sample_ids = []
for s in (s1, s2):
    print(f'{datetime.now()}: sapi.create_new_sample() for {s.external_id}')
    sample_ids.append(sapi.create_new_sample(PROJ, s))
print(f'{datetime.now()}: added samples {", ".join(sample_ids)}')
