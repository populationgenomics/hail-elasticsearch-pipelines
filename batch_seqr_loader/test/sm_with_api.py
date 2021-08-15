#!/usr/bin/env python3
# pylint: skip-file

"""
Reads input from the sample_metadata server using the generated API
"""

import os
from typing import Dict, List
import pandas as pd

from sample_metadata.api import SampleApi, AnalysisApi
from sample_metadata.models.new_sample import NewSample
from sample_metadata.models.analysis_model import AnalysisModel


PROJ = os.environ.get('SM_DEV_DB_PROJECT', 'vladdev')


def main():
    check_sample()


def add_sample():
    sapi = SampleApi()
    aapi = AnalysisApi()

    new_sample = NewSample(
        external_id='MYSAMPLE2', type='blood', meta={'is_test': 'true'}
    )
    sample_id = sapi.create_new_sample(PROJ, new_sample)

    analysis = AnalysisModel(
        sample_ids=[sample_id],
        type='gvcf',
        output='gs://mysample.g.vcf',
        status='completed',
    )

    analysis_id = aapi.create_new_analysis(PROJ, analysis)

    print(
        f'Inserted sample with ID: {sample_id}, inserted analysis with ID: {analysis_id}'
    )


def check_sample():
    sapi = SampleApi()
    s = sapi.get_sample_by_external_id('MYSAMPLE2', project=PROJ)
    print(s)


def list_samples():
    sapi = SampleApi()
    samples = sapi.get_samples(
        body_get_samples_by_criteria_api_v1_sample_post={
            'project_ids': ['acute-care'],
            'active': True,
        }
    )
    print(samples)


def get_incomplete_analysis():
    aapi = AnalysisApi()
    incomplete_analyses = []
    for t in ['cram', 'gvcf', 'joint-calling']:
        for status in ['queued', 'in-progress']:
            a = aapi.get_analysis_by_type_and_status(t, status, PROJ)
            incomplete_analyses.append(a)
    return incomplete_analyses


def start_analysis():
    data: Dict[str, List] = {
        'Family.ID': [],
        's': [],
        'Paternal.ID': [],
        'Maternal.ID': [],
        'Sex': [],
        'Phenotype': [],
        'file': [],
        'file2': [],
        'index': [],
        'type': [],
    }

    incomplete_analysis = get_incomplete_analysis()
    if incomplete_analysis:
        print(f'Error: incomplete analysis exist: {incomplete_analysis}')
    else:
        sapi = SampleApi()
        samples = sapi.get_samples('dev')
        samples_with_cram = [s for s in samples if 'cram' in s['meta']]
        for s in samples_with_cram:
            file = s['meta']['cram']
            index = file + '.crai'
            data['s'].append(s['id'])
            data['Family.ID'].append(s['id'])
            data['Paternal.ID'].append('0')
            data['Maternal.ID'].append('0')
            data['Sex'].append('0')
            data['Phenotype'].append('0')
            data['file'].append(file)
            data['file2'].append(None)
            data['index'].append(index)
            data['type'].append('cram')
    df = pd.DataFrame(data=data).set_index('s', drop=False)
    print(df)


if __name__ == '__main__':
    list_samples()
    # start_analysis()
