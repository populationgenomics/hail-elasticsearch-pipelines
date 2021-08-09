#!/usr/bin/env python3

"""
Simulate communication with the SM server during the workflow
"""

import os
import string
import random
from collections import defaultdict
from datetime import datetime
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
    print(f'{datetime.now()} Test run ID: {test_run_id}')

    print(f'{datetime.now()} Populate samples for test run {test_run_id}')
    sample_ids = _jc_pipeline_add_samples(test_run_id)
    print()

    print(f'{datetime.now()} Add/update analyses, reads -> cram')
    _jc_pipeline_submit_analyses()
    print()
    print(f'{datetime.now()} Set to in progress')
    _jc_pipeline_set_in_progress()
    print()
    print(f'{datetime.now()} Set to completed')
    _jc_pipeline_set_completed()
    print()

    print(f'{datetime.now()} Add/update analyses, cram -> gvcf')
    _jc_pipeline_submit_analyses()
    print()
    print(f'{datetime.now()} Set to in progress')
    _jc_pipeline_set_in_progress()
    print()
    print(f'{datetime.now()} Set to completed')
    _jc_pipeline_set_completed()
    print()

    print('Add/update analyses, gvcf -> joint-calling')
    _jc_pipeline_submit_analyses()
    print()
    print(f'{datetime.now()} Set to in progress')
    _jc_pipeline_set_in_progress()
    print()
    print(f'{datetime.now()} Set to completed')
    _jc_pipeline_set_completed()
    print()

    # Checking that after all calls, a 'completed' 'joint-calling' analysis must exist
    # for the initally added samples
    aapi = AnalysisApi()
    analyses = aapi.get_latest_complete_analyses(project=PROJ)
    print(f'Final analyses: {analyses}')
    assert any(
        a['type'] == 'joint-calling'
        and set(sample_ids) & set(sample_id_format(a['sample_ids'])) == set(sample_ids)
        for a in analyses
    ), [
        (a['type'], set(sample_id_format(a['sample_ids'])), set(sample_ids))
        for a in analyses
    ]


def _jc_pipeline_add_samples(test_run_id: str):
    """
    Add 3 samples: one with fastq input, one with CRAM input, one with GVCF input.
    :param test_run_id: to suffix sample names for uniqueness
    """
    sapi = SampleApi()
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
    s2 = NewSample(
        external_id=f'NA12878-from-cram-{test_run_id}',
        type='blood',
        meta={'reads': 'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878.cram'},
    )
    s3 = NewSample(
        external_id=f'NA12878-from-gvcf-{test_run_id}',
        type='blood',
        meta={'reads': 'gs://cpg-seqr-test/batches/NA12878-trio/NA12878.g.vcf.gz'},
    )
    sample_ids = []
    for s in (s1, s2, s3):
        print(f'{datetime.now()}: sapi.create_new_sample() for {s.external_id}')
        sample_ids.append(sapi.create_new_sample(PROJ, s))
    print(f'{datetime.now()}: added samples {", ".join(sample_ids)}')
    return sample_ids


def _jc_pipeline_submit_analyses():
    """
    Add or update analyses. Iterate over completed analyses,
    and submit next-step analyses
    """

    # If there are incomplete analyses, throw an error
    aapi = AnalysisApi()
    print(f'{datetime.now()}: starting _jc_pipeline_submit_analyses()')

    # analyses = aapi.get_incomplete_analyses(project=PROJ)
    # if analyses:
    #     print(f'ERROR: found incomplete or queued analysis: {analyses}')
    #     sys.exit()

    # Get the list of latest complete analyses
    print(f'{datetime.now()}: call aapi.get_latest_complete_analyses()')
    latest_complete_analyses = aapi.get_latest_complete_analyses(project=PROJ)
    print(f'{datetime.now()}: Latest complete analyses: {latest_complete_analyses}')
    latest_by_type_and_sids = defaultdict(list)
    for a in latest_complete_analyses:
        a_s_ids = sample_id_format(a['sample_ids'])
        latest_by_type_and_sids[(a['type'], tuple(set(a_s_ids)))].append(a)

    # Iterate over samples, check latest complete analyses, and add next-step analyses
    print(f'{datetime.now()}: sapi = SampleApi()')
    sapi = SampleApi()
    print(f'{datetime.now()}: sapi.get_all_samples(project=PROJ)')
    samples = sapi.get_all_samples(project=PROJ)

    if latest_by_type_and_sids.get(
        ('joint-calling', tuple(set(s.id for s in samples)))
    ):
        print(
            f'{datetime.now()}: All samples went through joint-calling, nothing to submit'
        )
        return

    print(f'{datetime.now()}: call aapi.get_latest_complete_analyses_by_type')
    latest_complete_gvcf_analyses = aapi.get_latest_complete_analyses_by_type(
        project=PROJ, analysis_type='gvcf'
    )
    print(f'{datetime.now()}: done aapi.get_latest_complete_analyses_by_type')
    sids_with_gvcf = set(
        sample_id_format(a['sample_ids'])[0] for a in latest_complete_gvcf_analyses
    )
    new_sids_with_gvcf = set(s.id for s in samples) - sids_with_gvcf
    if not new_sids_with_gvcf:
        print(
            f'{datetime.now()}: All samples went through variant calling, so can submit joint-calling'
        )
        analysis = AnalysisModel(
            sample_ids=[s.id for s in samples],
            type='joint-calling',
            output='gs://my-bucket/joint-calling/joint-called.g.vcf.gz',
            status='queued',
        )
        print(
            f'{datetime.now()}: Queueing {analysis.type} (call aapi.create_new_analysis)'
        )
        aapi.create_new_analysis(project=PROJ, analysis_model=analysis)
        print(f'{datetime.now()}: done aapi.create_new_analysis')
        return

    for s in [s for s in samples if s.id in new_sids_with_gvcf]:
        print(f'{datetime.now()}: Sample {s.id}')

        if latest_by_type_and_sids.get(('gvcf', (s.id,))):
            print('{datetime.now()}:   Sample has a complete gvcf analysis')

        elif latest_by_type_and_sids.get(('cram', (s.id,))):
            print(
                f'{datetime.now()}:   Sample has a complete CRAM analysis, queueing variant calling'
            )
            analysis = AnalysisModel(
                sample_ids=[s.id],
                type='gvcf',
                output=f'gs://my-bucket/variant-calling/{s.id}.g.vcf.gz',
                status='queued',
            )
            print(f'{datetime.now()}:   call aapi.create_new_analysis')
            aapi.create_new_analysis(project=PROJ, analysis_model=analysis)
            print(f'{datetime.now()}:   done aapi.create_new_analysis')

        else:
            print(
                f'{datetime.now()}:     Sample doesn not have any analysis yet, trying to get "reads" '
                'metadata to submit alignment'
            )
            reads_data = s.meta.get('reads')
            if not reads_data:
                print(f'{datetime.now()}:     ERROR: no "reads" data')
            elif isinstance(reads_data, str):
                if reads_data.endswith('.g.vcf.gz'):
                    analysis = AnalysisModel(
                        sample_ids=[s.id],
                        type='gvcf',
                        output=reads_data,
                        status='completed',
                    )
                    print(f'{datetime.now()}:   call aapi.create_new_analysis')
                    aapi.create_new_analysis(project=PROJ, analysis_model=analysis)
                    print(f'{datetime.now()}:   done aapi.create_new_analysis')
                elif reads_data.endswith('.cram') or reads_data.endswith('.bam'):
                    print(f'{datetime.now()}:     Queueing cram re-alignment analysis')
                    analysis = AnalysisModel(
                        sample_ids=[s.id],
                        type='cram',
                        output=f'gs://my-bucket/realignment/{s.id}.cram',
                        status='queued',
                    )
                    print(f'{datetime.now()}:   call aapi.create_new_analysis')
                    aapi.create_new_analysis(project=PROJ, analysis_model=analysis)
                    print(f'{datetime.now()}:   done aapi.create_new_analysis')
                else:
                    print(
                        f'{datetime.now()}:     ERROR: unrecognised "reads" meta data: {reads_data}'
                    )
            elif isinstance(reads_data, list) and len(reads_data) == 2:
                print(f'{datetime.now()}:     Queueing cram alignment analyses')
                analysis = AnalysisModel(
                    sample_ids=[s.id],
                    type='cram',
                    output=f'gs://my-bucket/alignment/{s.id}.cram',
                    status='queued',
                )
                print(f'{datetime.now()}:   call aapi.create_new_analysis')
                aapi.create_new_analysis(project=PROJ, analysis_model=analysis)
                print(f'{datetime.now()}:   done aapi.create_new_analysis')
            else:
                print(
                    f'{datetime.now()}:     ERROR: can\'t recognise "reads" data: {reads_data}'
                )


def _jc_pipeline_set_in_progress():
    """
    Update existing queued analyses and set their status to in-progress
    """
    aapi = AnalysisApi()
    print(f'{datetime.now()}: call aapi.get_incomplete_analyses')
    analyses = aapi.get_incomplete_analyses(project=PROJ)
    print(f'{datetime.now()}: done aapi.get_incomplete_analyses')
    if analyses:
        for a in analyses:
            print(f'{datetime.now()}: Setting analysis {a} to in-progress')
            aum = AnalysisUpdateModel(status='in-progress')
            print(f'{datetime.now()}: call aapi.update_analysis_status')
            aapi.update_analysis_status(
                analysis_id=a['id'],
                project=PROJ,
                analysis_update_model=aum,
            )
            print(f'{datetime.now()}: done aapi.update_analysis_status')


def _jc_pipeline_set_completed():
    """
    Update existing in-progress analyses and set their status to completed
    """
    aapi = AnalysisApi()
    print(f'{datetime.now()}: call aapi.get_incomplete_analyses')
    analyses = aapi.get_incomplete_analyses(project=PROJ)
    print(f'{datetime.now()}: done aapi.get_incomplete_analyses')
    if analyses:
        for a in analyses:
            print(f'{datetime.now()}: Setting analysis {a} to completed')
            aum = AnalysisUpdateModel(status='completed')
            print(f'{datetime.now()}: call aapi.update_analysis_status')
            aapi.update_analysis_status(
                analysis_id=a['id'],
                project=PROJ,
                analysis_update_model=aum,
            )
            print(f'{datetime.now()}: done aapi.update_analysis_status')


def sample_id_format(sample_id: Union[int, List[int]]):
    """
    Transform raw (int) sample identifier to format (CPGXXXH) where:
        - CPG is the prefix
        - H is the Luhn checksum
        - XXX is the original identifier

    >>> sample_id_format(10)
    'CPG109'

    >>> sample_id_format(12345)
    'CPG123455'
    """

    if isinstance(sample_id, list):
        return [sample_id_format(s) for s in sample_id]

    if isinstance(sample_id, str) and not sample_id.isdigit():
        if sample_id.startswith('CPG'):
            return sample_id
        raise ValueError(f'Unexpected format for sample identifier "{sample_id}"')
    sample_id = int(sample_id)

    return f'CPG{sample_id}{luhn_compute(sample_id)}'


def luhn_compute(n):
    """
    Compute Luhn check digit of number given as string

    >>> luhn_compute(453201511283036)
    6

    >>> luhn_compute(601151443354620)
    1

    >>> luhn_compute(677154949558680)
    2
    """
    m = [int(d) for d in reversed(str(n))]
    result = sum(m) + sum(d + (d >= 5) for d in m[::2])
    return -result % 10


if __name__ == '__main__':
    test_simulate_joint_calling_pipeline()
