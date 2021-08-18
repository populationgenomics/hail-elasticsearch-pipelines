#!/usr/bin/env python3
# pylint: skip-file

import hailtop.batch as hb
from os.path import join, splitext, dirname
import os
from sample_metadata import SampleApi
from google.cloud import storage
from typing import Optional, List, Dict, Tuple
from dataclasses import dataclass
import logging

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


BWA_CONTAINER_1 = (
    f'australia-southeast1-docker.pkg.dev/cpg-common/images/biobambam2:debug-asan'
)
BWA_CONTAINER_2 = (
    f'australia-southeast1-docker.pkg.dev/cpg-common/images/biobambam2:debug-tsan'
)
REF_BUCKET = 'gs://cpg-reference/hg38/v1'
REF_FASTA = join(REF_BUCKET, 'Homo_sapiens_assembly38.fasta')


@dataclass
class AlignmentInput:
    """
    Sort of a union type for possible alignment inputs
    """

    bam_or_cram_path: Optional[str] = None
    index_path: Optional[str] = None
    fqs1: Optional[List[str]] = None
    fqs2: Optional[List[str]] = None


def _test_bwa(
    b: hb.Batch,
    reference: hb.ResourceGroup,
    container: str,
    sample_name: str,
    alignment_input: AlignmentInput,
):
    j = b.new_job(f'Test BWA with {container}')
    j.image(container)
    total_cpu = 32
    output_cram_path = f'gs://cpg-seqr-test-tmp/test-biobambam2/{sample_name}-{container.split(":")[1]}.cram'

    if alignment_input.bam_or_cram_path:
        use_bazam = True
        bazam_cpu = 10
        bwa_cpu = 32
        bamsormadup_cpu = 10

        assert alignment_input.index_path
        assert not alignment_input.fqs1 and not alignment_input.fqs2
        cram = b.read_input_group(
            base=alignment_input.bam_or_cram_path, index=alignment_input.index_path
        )
        r1_param = (
            f'<(bazam -Xmx16g -Dsamjdk.reference_fasta={reference.base}'
            f' -n{bazam_cpu} -bam {cram.base})'
        )
        r2_param = '-'
    else:
        assert alignment_input.fqs1 and alignment_input.fqs2
        use_bazam = False
        bwa_cpu = 32
        bamsormadup_cpu = 10
        files1 = [b.read_input(f1) for f1 in alignment_input.fqs1]
        files2 = [b.read_input(f1) for f1 in alignment_input.fqs2]
        r1_param = f'<(cat {" ".join(files1)})'
        r2_param = f'<(cat {" ".join(files2)})'
        logger.info(f'r1_param: {r1_param}')
        logger.info(f'r2_param: {r2_param}')

    j.cpu(total_cpu)
    j.memory('standard')
    j.storage('300G')
    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'crai': '{root}.crai',
        }
    )

    rg_line = f'@RG\\tID:{sample_name}\\tSM:{sample_name}'
    # BWA command options:
    # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
    # -p     smart pairing (ignoring in2.fq)
    # -t16   threads
    # -Y     use soft clipping for supplementary alignments
    # -R     read group header line such as '@RG\tID:foo\tSM:bar'
    command = f"""
set -o pipefail
set -ex

(while true; do df -h; pwd; du -sh $(dirname {j.output_cram.cram}); sleep 600; done) &

bwa mem -K 100000000 {'-p' if use_bazam else ''} -t{bwa_cpu} -Y \\
    -R '{rg_line}' {reference.base} {r1_param} {r2_param} | \\
bamsormadup inputformat=sam threads={bamsormadup_cpu} SO=coordinate \\
    M={j.duplicate_metrics} outputformat=sam \\
    tmpfile=$(dirname {j.output_cram.cram})/bamsormadup-tmp | \\
samtools view -T {reference.base} -O cram -o {j.output_cram.cram}

samtools index -@{total_cpu} {j.output_cram.cram} {j.output_cram.crai}

df -h; pwd; du -sh $(dirname {j.output_cram.cram})
    """
    j.command(command)
    b.write_output(j.output_cram, splitext(output_cram_path)[0])
    return j, output_cram_path


def file_exists(path: str) -> bool:
    """
    Check if the object exists, where the object can be:
        * local file
        * local directory
        * Google Storage object
        * Google Storage URL representing a *.mt or *.ht Hail data,
          in which case it will check for the existence of a
          *.mt/_SUCCESS or *.ht/_SUCCESS file.
    :param path: path to the file/directory/object/mt/ht
    :return: True if the object exists
    """
    if path.startswith('gs://'):
        bucket = path.replace('gs://', '').split('/')[0]
        path = path.replace('gs://', '').split('/', maxsplit=1)[1]
        path = path.rstrip('/')  # '.mt/' -> '.mt'
        if any(path.endswith(f'.{suf}') for suf in ['mt', 'ht']):
            path = os.path.join(path, '_SUCCESS')
        gs = storage.Client()
        return gs.get_bucket(bucket).get_blob(path)
    return os.path.exists(path)


def sm_verify_reads_data(  # pylint: disable=too-many-return-statements
    reads_data: Optional[List],
    reads_type: Optional[str],
) -> Optional[AlignmentInput]:
    """
    Verify the meta.reads object in a sample db entry
    """
    if not reads_data:
        logger.error(f'ERROR: no "meta/reads" field')
        return None
    if not reads_type:
        logger.error(f'ERROR: no "meta/reads_type" field')
        return None
    supported_types = ('fastq', 'bam', 'cram')
    if reads_type not in supported_types:
        logger.error(f'ERROR: "reads_type" is expected to be one of {supported_types}')
        return None

    if reads_type in ('bam', 'cram'):
        if len(reads_data) > 1:
            logger.error('Supporting only single bam/cram input')
            return None

        bam_path = reads_data[0]['location']
        if not (bam_path.endswith('.cram') or bam_path.endswith('.bam')):
            logger.error(
                f'ERROR: expected the file to have an extention .cram or .bam,'
                f'got: {bam_path}'
            )
            return None
        if not file_exists(bam_path):
            logger.error(f'ERROR: index file doesn\'t exist: {bam_path}')
            return None

        # Index:
        if not reads_data[0].get('secondaryFiles'):
            logger.error(
                f'ERROR: bam/cram input is expected to have '
                f'a non-empty list field "secondaryFile" section with indices'
            )
            return None
        index_path = reads_data[0]['secondaryFiles'][0]['location']
        if (
            bam_path.endswith('.cram')
            and not index_path.endswith('.crai')
            or bam_path.endswith('.bai')
            and not index_path.endswith('.bai')
        ):
            logger.error(
                f'ERROR: expected the index file to have an extention '
                f'.crai or .bai, got: {index_path}'
            )
        if not file_exists(index_path):
            logger.error(f'ERROR: index file doesn\'t exist: {index_path}')
            return None

        return AlignmentInput(bam_or_cram_path=bam_path, index_path=index_path)

    else:
        fqs1 = []
        fqs2 = []
        for lane_data in reads_data:
            assert len(lane_data) == 2
            if not file_exists(lane_data[0]['location']):
                logger.error(
                    f'ERROR: read 1 file doesn\'t exist: {lane_data[0]["location"]}'
                )
                return None
            if not file_exists(lane_data[1]['location']):
                logger.error(
                    f'ERROR: read 2 file doesn\'t exist: {lane_data[1]["location"]}'
                )
                return None

            fqs1.append(lane_data[0]['location'])
            fqs2.append(lane_data[1]['location'])
        return AlignmentInput(fqs1=fqs1, fqs2=fqs2)


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

bwa_reference = b.read_input_group(
    base=REF_FASTA,
    fai=REF_FASTA + '.fai',
    dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
    + '.dict',
    sa=REF_FASTA + '.sa',
    amb=REF_FASTA + '.amb',
    bwt=REF_FASTA + '.bwt',
    ann=REF_FASTA + '.ann',
    pac=REF_FASTA + '.pac',
)

sapi = SampleApi()
samples = []
for s in sapi.get_samples(
    body_get_samples_by_criteria_api_v1_sample_post={
        'project_ids': ['seqr-test'],
        'active': True,
    }
):
    if s['external_id'] in ['CPG12229', 'CPG12302', 'CPG11981', 'CPG11817']:
        print(
            f"Processing sample {s['id']}/{s['external_id']}, with metadata {s['meta']}"
        )
        samples.append(s)

for s in samples:
    for cont in [BWA_CONTAINER_1, BWA_CONTAINER_2]:
        alignment_input = sm_verify_reads_data(
            s['meta'].get('reads'), s['meta'].get('reads_type')
        )
        logger.info(f'Submitting {alignment_input}')
        _test_bwa(b, bwa_reference, cont, s['external_id'], alignment_input)
b.run(open=True)
