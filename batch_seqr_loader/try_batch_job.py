#!/usr/bin/env python3

import hailtop.batch as hb
from os.path import join
import os

# import hailtop.batch as hb
# backend = hb.ServiceBackend()
# batch = hb.Batch(backend=backend, name='test-batch')

# j = batch.new_job(name='test-job')
# j.image('ubuntu:18.04')
# j.command(f'''
# echo "test" > {j.ofile}
# ''')
# batch.write_output(j.ofile, 'gs://playground-us/batch-test/test.done')
# batch.run(open=True, wait=False)
# backend.close()

BAZAM_CONTAINER = f'australia-southeast1-docker.pkg.dev/cpg-common/images/bazam:v2'
REF_BUCKET = 'gs://cpg-reference/hg38/v1/ref-tiny'
REF_FASTA = join(REF_BUCKET, 'Homo_sapiens_assembly38-tiny.fasta')


def _index_bwa_job(
    b: hb.Batch,
    reference: hb.ResourceGroup
):
    exts = ['sa', 'amb', 'bwt', 'ann', 'pac']
    j = b.new_job('Index BWA')
    j.image(BAZAM_CONTAINER)
    total_cpu = 16
    j.cpu(total_cpu)
    j.storage('40G')
    j.declare_resource_group(bwa_index={e: '{root}.' + e for e in exts})
    j.command(
        f"""
set -o pipefail
set -ex

bwa index {reference.base} -p {j.bwa_index}

df -h; pwd; ls | grep -v proc | xargs du -sh
    """
    )
    b.write_output(j.bwa_index, REF_FASTA)
    return j


def _make_realign_jobs(
    b: hb.Batch,
    reference: hb.ResourceGroup,
):
    fq1 = b.read_input('gs://cpg-seqr-test/batches/test/tmp_fq')
    j = b.new_job('Test BWA')
    j.image(BAZAM_CONTAINER)
    total_cpu = 16
    bazam_cpu = 3
    bwa_cpu = 10
    bamsormadup_cpu = 3
    bwa_cpu = 1
    j.memory('highmem')
    j.cpu(total_cpu)
    j.storage('300G')
    j.declare_resource_group(
        output_cram={
            'base': '{root}.cram',
            'crai': '{root}.crai',
        }
    )
    sn = 'TEST'
    rg_line = f'@RG\\tID:{sn}\\tSM:~{sn}'
    # bwa index {reference.base}
    use_bazam = True

    j.command(
        f"""
set -o pipefail
set -ex

(while true; do df -h; pwd; ls | grep -v proc | xargs du -sh; free -m; sleep 300; done) &

bwa mem -K 100000000 {'-p' if use_bazam else ''} -v3 -t{bwa_cpu} -Y \\
-R '{rg_line}' {reference.base} \\
{fq1} - > {j.aligned_sam}

bamsormadup inputformat=sam threads={bamsormadup_cpu} SO=coordinate \\
M={j.duplicate_metrics} outputformat=sam < {j.aligned_sam} | \\
samtools view -T {reference.base} -O cram -o {j.output_cram.base}

samtools index -@{total_cpu} {j.output_cram.base} {j.output_cram.crai}

df -h; pwd; ls | grep -v proc | xargs du -sh
    """
    )

billing_project = os.getenv('HAIL_BILLING_PROJECT') or 'seqr'
hail_bucket = os.environ.get('HAIL_BUCKET')
print(
    f'Starting hail Batch with the project {billing_project}, '
    f'bucket {hail_bucket}'
)
backend = hb.ServiceBackend(
    billing_project=billing_project,
    bucket=hail_bucket.replace('gs://', ''),
)
b = hb.Batch(backend=backend, name='test')
reference = b.read_input_group(
    base=REF_FASTA,
    fai=REF_FASTA + '.fai',
    dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
    + '.dict',
)
_index_bwa_job(b, reference)
# _make_realign_jobs(b, reference)
b.run(open=True)
