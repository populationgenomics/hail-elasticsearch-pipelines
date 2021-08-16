#!/usr/bin/env python3

"""
Simulate communication with the SM server during the workflow
"""

import os
from datetime import datetime
from sample_metadata.api import SampleApi
from sample_metadata.models.new_sample import NewSample


PROJ = os.environ.get('SM_DEV_DB_PROJECT', 'seqr-test')

sapi = SampleApi()


def add_large():
    """
    Add 2 full WGS samples (from the acute-care project), one with a cram input,
    another with fastq
    """
    s1 = NewSample(
        external_id='210303_A00692_0190_ML211637-1',
        type='blood',
        meta={
            'reads_type': 'fastq',
            'reads': [
                [
                    {
                        'basename': '210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L001_R1.fastq.gz',
                        'class': 'File',
                        'location': 'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L001_R1.fastq.gz',
                    },
                    {
                        'basename': '210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L001_R2.fastq.gz',
                        'class': 'File',
                        'location': 'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L001_R2.fastq.gz',
                    },
                ],
                [
                    {
                        'basename': '210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L002_R1.fastq.gz',
                        'class': 'File',
                        'location': 'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L002_R1.fastq.gz',
                    },
                    {
                        'basename': '210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L002_R2.fastq.gzz',
                        'class': 'File',
                        'location': 'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_L002_R2.fastq.gz',
                    },
                ],
            ],
        },
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
    s1 = NewSample(
        external_id='NA12878-fastq',
        type='blood',
        meta={
            'reads_type': 'fastq',
            'reads': [
                [
                    {
                        'basename': 'NA12878_L001_R1.fq',
                        'class': 'File',
                        'location': 'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878_L001_R1.fq',
                    },
                    {
                        'basename': 'NA12878_L001_R2.fq',
                        'class': 'File',
                        'location': 'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878_L001_R2.fq',
                    },
                ],
                [
                    {
                        'basename': 'NA12878_L002_R1.fq',
                        'class': 'File',
                        'location': 'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878_L002_R1.fq',
                    },
                    {
                        'basename': 'NA12878_L002_R2.fq',
                        'class': 'File',
                        'location': 'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878_L002_R2.fq',
                    },
                ],
            ],
        },
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
    add_tiny()
