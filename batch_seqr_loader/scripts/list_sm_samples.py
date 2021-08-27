#!/usr/bin/env python3
# pylint: skip-file

"""
List samples in the SM server DB
"""
from sample_metadata import SampleUpdateModel
from sample_metadata.api import SampleApi, AnalysisApi

sapi = SampleApi()
aapi = AnalysisApi()

project = 'seqr-test'

samples = sapi.get_samples(
    body_get_samples_by_criteria_api_v1_sample_post={
        'project_ids': [project],
        'active': True,
    }
)

for s in samples:
    print(s)
    if s['external_id'] == 'NA12878-fastq':
        sapi.update_sample(s['id'], sample_update_model=SampleUpdateModel(active=False))
