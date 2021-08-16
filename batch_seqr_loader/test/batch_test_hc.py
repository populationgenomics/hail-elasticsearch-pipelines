#!/usr/bin/env python3
# pylint: skip-file
import re

import hailtop.batch as hb
from os.path import join, basename, splitext
import os
from typing import Optional, List, Tuple, Set, Dict
from hailtop.batch.job import Job


AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
GATK_VERSION = '4.2.0.0'
GATK_IMAGE = f'{AR_REPO}/gatk:{GATK_VERSION}'
PICARD_IMAGE = f'{AR_REPO}/picard-cloud:2.23.8'

REF_BUCKET = 'gs://cpg-reference/hg38/v1'
REF_FASTA = join(REF_BUCKET, 'Homo_sapiens_assembly38.fasta')
UNPADDED_INTERVALS = join(REF_BUCKET, 'hg38.even.handcurated.20k.intervals')

CRAM_PATH = 'gs://cpg-seqr-test/datasets/acute-care/v1-14/CPG13185.cram'
TMP_BUCKET = 'gs://cpg-seqr-test-tmp/seqr-loader/test-hc'


def _add_split_intervals_job(
    b: hb.Batch,
    interval_list: str,
    scatter_count: int,
    ref_fasta: str,
    method: str = 'INTERVAL_SUBDIVISION',
) -> Job:
    """
    Split genome into intervals to parallelise GnarlyGenotyper.

    Returns: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job(f'Make {scatter_count} intervals, method: {method}')
    j.image(GATK_IMAGE)
    java_mem = 3
    j.memory('standard')  # ~ 4G/core ~ 4G
    j.storage('16G')
    j.declare_resource_group(
        intervals={
            f'interval_{idx}': f'{{root}}/{str(idx).zfill(4)}-scattered.interval_list'
            for idx in range(scatter_count)
        }
    )

    j.command(
        f"""set -e

    # Modes other than INTERVAL_SUBDIVISION will produce an unpredicted number 
    # of intervals. But we have to expect exactly the {scatter_count} number of 
    # output files because our workflow is not dynamic.
    gatk --java-options -Xms{java_mem}g SplitIntervals \\
      -L {interval_list} \\
      -O {j.intervals} \\
      -scatter {scatter_count} \\
      -R {ref_fasta} \\
      -mode {method}
      """
    )
    return j


def _add_haplotype_caller_job(
    b: hb.Batch,
    cram: hb.ResourceGroup,
    interval: hb.ResourceFile,
    reference: hb.ResourceGroup,
    sample_name: str,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    depends_on: Optional[List[Job]] = None,
    output_gvcf_path: Optional[str] = None,
) -> Job:
    """
    Run HaplotypeCaller on an input BAM or CRAM, and output GVCF
    """
    job_name = 'HaplotypeCaller'
    if interval_idx is not None:
        job_name += f', {sample_name} {interval_idx}/{number_of_intervals}'

    j = b.new_job(job_name)
    j.image(GATK_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage('60G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )
    if depends_on:
        j.depends_on(*depends_on)

    j.command(
        f"""set -e
    (while true; do df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m; sleep 300; done) &

    gatk --java-options "-Xms{java_mem}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
      HaplotypeCaller \\
      -R {reference.base} \\
      -I {cram['cram']} \\
      -L {interval} \\
      -O {j.output_gvcf['g.vcf.gz']} \\
      -G AS_StandardAnnotation \\
      -GQB 20 \
      -ERC GVCF \\

    df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m
    """
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, j.output_gvcf)
    return j


def _add_hc_jobs(b: hb.Batch, reference: hb.ResourceGroup) -> List[Job]:
    sample_name = splitext(basename(CRAM_PATH))[0]

    merge_jobs = []

    for method in [
        'INTERVAL_SUBDIVISION',
        'BALANCING_WITHOUT_INTERVAL_SUBDIVISION',
        'BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW',
    ]:
        for n_intervals in [1, 2, 10, 50, 100]:
            gvcf_name = f'{sample_name}_{n_intervals}_{method}'

            intervals_j = _add_split_intervals_job(
                b=b,
                interval_list=UNPADDED_INTERVALS,
                scatter_count=n_intervals,
                ref_fasta=REF_FASTA,
                method=method,
            )
            haplotype_caller_jobs = []
            for idx in range(n_intervals):
                haplotype_caller_jobs.append(
                    _add_haplotype_caller_job(
                        b,
                        cram=b.read_input_group(
                            **{
                                'cram': CRAM_PATH,
                                'crai': re.sub('.cram$', '.crai', CRAM_PATH),
                            }
                        ),
                        interval=intervals_j.intervals[f'interval_{idx}'],
                        reference=reference,
                        sample_name=gvcf_name,
                        interval_idx=idx,
                        number_of_intervals=n_intervals,
                        depends_on=intervals_j,
                    )
                )
            hc_gvcf_path = join(TMP_BUCKET, f'{gvcf_name}.g.vcf.gz')
            merge_j = _add_merge_gvcfs_job(
                b,
                gvcfs=[j.output_gvcf for j in haplotype_caller_jobs],
                output_gvcf_path=hc_gvcf_path,
                sample_name=sample_name,
            )
            merge_jobs.append(merge_j)
    return merge_jobs


def _add_merge_gvcfs_job(
    b: hb.Batch,
    gvcfs: List[hb.ResourceGroup],
    output_gvcf_path: Optional[str],
    sample_name: str,
) -> Job:
    """
    Combine by-interval GVCFs into a single sample GVCF file
    """

    job_name = f'Merge {len(gvcfs)} GVCFs, {sample_name}'
    j = b.new_job(job_name)
    j.image(PICARD_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage(f'{len(gvcfs) * 1.5 + 2}G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )

    input_cmd = ' '.join(f'INPUT={g["g.vcf.gz"]}' for g in gvcfs)

    j.command(
        f"""set -e

    (while true; do df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m; sleep 300; done) &

    java -Xms{java_mem}g -jar /usr/picard/picard.jar \
      MergeVcfs {input_cmd} OUTPUT={j.output_gvcf['g.vcf.gz']}

    df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m
      """
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, output_gvcf_path.replace('.g.vcf.gz', ''))
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
reference = b.read_input_group(
    base=REF_FASTA,
    fai=REF_FASTA + '.fai',
    dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
    + '.dict',
)
_add_hc_jobs(b, reference)
b.run(open=True)
