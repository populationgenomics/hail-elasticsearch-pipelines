#!/usr/bin/env python3
# pylint: skip-file

import hailtop.batch as hb
import os


SEQTK_IMAGE = 'staphb/seqtk'
SAMTOOLS_IMAGE = 'biocontainers/samtools'


def _subset_fq(b: hb.Batch):
    j = b.new_job('Subset fastq')
    j.image(SEQTK_IMAGE)
    j.storage('200G')
    fq_src_base = 'gs://cpg-acute-care-test-upload/cpg_acute_20210729_084835/210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_'
    fq_trg_base = 'gs://cpg-acute-care-test-upload/test-tiny/210303_A00692_0190_ML211637_21W000380-FAM001171_MAN-20210303_NEXTERAFLEXWGS_'
    fq_l1_r1 = b.read_input(f'{fq_src_base}L001_R1.fastq.gz')
    fq_l1_r2 = b.read_input(f'{fq_src_base}L001_R2.fastq.gz')
    fq_l2_r1 = b.read_input(f'{fq_src_base}L002_R1.fastq.gz')
    fq_l2_r2 = b.read_input(f'{fq_src_base}L002_R2.fastq.gz')
    j.command(
        f"""
set -o pipefail
set -ex

seqtk sample -s100 {fq_l1_r1} 1000000 > {j.l1_r1}
seqtk sample -s100 {fq_l1_r2} 1000000 > {j.l1_r2}
seqtk sample -s100 {fq_l2_r1} 1000000 > {j.l2_r1}
seqtk sample -s100 {fq_l2_r2} 1000000 > {j.l2_r2}
    """
    )
    b.write_output(j.l1_r1, f'{fq_trg_base}L001_R1.fastq.gz')
    b.write_output(j.l1_r2, f'{fq_trg_base}L001_R2.fastq.gz')
    b.write_output(j.l2_r1, f'{fq_trg_base}L002_R1.fastq.gz')
    b.write_output(j.l2_r2, f'{fq_trg_base}L002_R2.fastq.gz')
    return j


def _subset_cram(b: hb.Batch):
    j = b.new_job('Subset CRAM')
    j.image(SAMTOOLS_IMAGE)
    j.storage('50G')
    cram_input = b.read_input(
        'gs://cpg-acute-care-test-upload/20W000094-FAM000347.cram'
    )
    j.command(
        f"""
set -o pipefail
set -ex

samtools view -s 0.1 sample -s100 {cram_input} 1000000 > {j.cram}
    """
    )
    b.write_output(
        j.cram, 'gs://cpg-acute-care-test-upload/test-tiny/20W000094-FAM000347.cram'
    )
    return j


billing_project = os.getenv('HAIL_BILLING_PROJECT') or 'seqr'
hail_bucket = os.environ.get('HAIL_BUCKET', 'cpg-seqr-test-tmp')
print(
    f'Starting hail Batch with the project {billing_project}, ' f'bucket {hail_bucket}'
)
backend = hb.ServiceBackend(
    billing_project=billing_project,
    bucket=hail_bucket.replace('gs://', ''),
)
b = hb.Batch(backend=backend, name='test')
_subset_fq(b)
_subset_cram(b)
b.run(open=True)
