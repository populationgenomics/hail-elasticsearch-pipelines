#!/usr/bin/env python3

"""
Simulate communication with the SM server during the workflow
"""

import os
from datetime import datetime
from os.path import basename

from sample_metadata.api import SampleApi
from sample_metadata.models.new_sample import NewSample


PROJ = os.environ.get('SM_DEV_DB_PROJECT', 'seqr-test')

sapi = SampleApi()


def _add_fq_sample(external_id, fqs1):
    s = NewSample(
        external_id=external_id,
        type='blood',
        meta={
            'reads_type': 'fastq',
            'reads': [
                [
                    {
                        'basename': basename(fq1),
                        'class': 'File',
                        'location': fq1,
                    },
                    {
                        'basename': basename(fq1.replace('_R1.', '_R2.')),
                        'class': 'File',
                        'location': fq1.replace('_R1.', '_R2.'),
                    },
                ]
                for fq1 in fqs1
            ],
        },
    )
    return s


def add_for_biobambam2():
    """
    Add 4 full WGS samples from the acute-care project that failed biobambam2 previously
    """
    samples = []
    for sid, fqs1 in [
        (
            'CPG11981',
            [
                'gs://cpg-acute-care-test-upload/cpg_acute_20210728_142236/201201_A00692_0166_ML208851_20W002140-FAM000943_MAN-20201201_NEXTERAFLEXWGS_L001_R1.fastq.gz'
            ],
        ),
        (
            'CPG11817',
            [
                'gs://cpg-acute-care-test-upload/cpg_acute_20210728_142236/210124_A01221_0019_ML210702_21W000122-FAM001099_MAN-20210123_NEXTERAFLEXWGS_L001_R1.fastq.gz',
                'gs://cpg-acute-care-test-upload/cpg_acute_20210728_142236/210124_A01221_0019_ML210702_21W000122-FAM001099_MAN-20210123_NEXTERAFLEXWGS_L002_R1.fastq.gz',
            ],
        ),
        (
            'CPG12302',
            [
                'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211641_21W000384-FAM001172_MAN-20210303_NEXTERAFLEXWGS_L001_R1.fastq.gz',
                'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211641_21W000384-FAM001172_MAN-20210303_NEXTERAFLEXWGS_L002_R1.fastq.gz',
            ],
        ),
        (
            'CPG12229',
            [
                'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210428_A01221_0037_ML212562_21W000847-FAM001327_MAN-20210428_NEXTERAFLEXWGS_L001_R1.fastq.gz',
                'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210428_A01221_0037_ML212562_21W000847-FAM001327_MAN-20210428_NEXTERAFLEXWGS_L002_R1.fastq.gz',
            ],
        ),
    ]:
        samples.append(_add_fq_sample(sid, fqs1))

    sample_ids = []
    for s in samples:
        print(f'{datetime.now()}: sapi.create_new_sample() for {s.external_id}')
        sample_ids.append(sapi.create_new_sample(PROJ, s))
    print(f'{datetime.now()}: added samples {", ".join(sample_ids)}')


def add_large():
    """
    Add 2 full WGS samples (from the acute-care project), one with a cram input,
    another with fastq
    """
    s1 = _add_fq_sample(
        '210303_A00692_0190_ML211637-1',
        [
            'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L001_R1.fastq.gz',
            'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L002_R1.fastq.gz',
        ],
    )

    s2 = NewSample(
        external_id='20W000094-FAM000347-1',
        type='blood',
        meta={
            'reads_type': 'bam',
            'reads': [
                {
                    'basename': '20W000094-FAM000347.cram',
                    'class': 'File',
                    'location': 'gs://cpg-acute-care-test-upload/20W000094-FAM000347.cram',
                }
            ],
        },
    )
    sample_ids = []
    for s in (s1, s2):
        print(f'{datetime.now()}: sapi.create_new_sample() for {s.external_id}')
        sample_ids.append(sapi.create_new_sample(PROJ, s))
    print(f'{datetime.now()}: added samples {", ".join(sample_ids)}')


def add_tiny():
    """
    Add 2 tiny samples, subset to chromosome 21. One with a cram input,
    another with fastq
    """
    s1 = _add_fq_sample(
        'NA12878-fastq',
        [
            'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878_L001_R1.fq',
            'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878_L002_R1.fq',
        ],
    )

    s2 = NewSample(
        external_id='NA12878-cram',
        type='blood',
        meta={
            'reads_type': 'bam',
            'reads': [
                {
                    'basename': 'NA12878.cram',
                    'class': 'File',
                    'location': 'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878.cram',
                }
            ],
        },
    )
    sample_ids = []
    for s in (s1, s2):
        print(f'{datetime.now()}: sapi.create_new_sample() for {s.external_id}')
        sample_ids.append(sapi.create_new_sample(PROJ, s))
    print(f'{datetime.now()}: added samples {", ".join(sample_ids)}')


if __name__ == '__main__':
    add_for_biobambam2()
