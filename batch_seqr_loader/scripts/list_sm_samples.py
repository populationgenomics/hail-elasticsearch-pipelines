#!/usr/bin/env python3
# pylint: skip-file

"""
List samples in the SM server DB
"""


import os
from sample_metadata.api import SampleApi, AnalysisApi

PROJ = os.environ.get('SM_DEV_DB_PROJECT', 'vladdev')

sapi = SampleApi()
aapi = AnalysisApi()

samples = sapi.get_samples()
projects = set()
for s in samples:
    projects.add(s['meta'].get('project'))
print(f'Projects: {projects}')
print('Samples:')
for s in samples:
    print(s)
