"""
Utility functions and constants for the seqr loader pipeline
"""

import hashlib
import os
import subprocess
from os.path import join
from typing import Iterable, Tuple, Optional, Dict
import logging
from google.cloud import storage
import hailtop.batch as hb
from hailtop.batch import Batch
from hailtop.batch.job import Job

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
GATK_VERSION = '4.2.2.0-cpgfix00'  # Our fork with a couple of fixes:
# https://github.com/populationgenomics/production-pipelines/tree/initial/dockers/gatk
GATK_IMAGE = f'{AR_REPO}/gatk:{GATK_VERSION}'
PICARD_IMAGE = f'{AR_REPO}/picard-cloud:2.23.8'
ALIGNMENT_IMAGE = f'{AR_REPO}/alignment:v4'
SOMALIER_IMAGE = f'{AR_REPO}/somalier:latest'
PEDDY_IMAGE = f'{AR_REPO}/peddy:0.4.8--pyh5e36f6f_0'
GNARLY_IMAGE = f'{AR_REPO}/gnarly_genotyper:hail_ukbb_300K'
BCFTOOLS_IMAGE = f'{AR_REPO}/bcftools:1.10.2--h4f4756c_2'
SM_IMAGE = f'{AR_REPO}/sm-api:2.0.5'

NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS = 50
NUMBER_OF_GENOMICS_DB_INTERVALS = 50
NUMBER_OF_DATAPROC_WORKERS = 50

REF_BUCKET = 'gs://cpg-reference/hg38'
NOALT_REGIONS = join(REF_BUCKET, 'noalt.bed')
SOMALIER_SITES = join(REF_BUCKET, 'somalier/v0/sites.hg38.vcf.gz')

GATK_REF_BUCKET = f'{REF_BUCKET}/v1'
REF_FASTA = join(GATK_REF_BUCKET, 'Homo_sapiens_assembly38.fasta')
DBSNP_VCF = join(GATK_REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf')
UNPADDED_INTERVALS = join(GATK_REF_BUCKET, 'hg38.even.handcurated.20k.intervals')

DATAPROC_PACKAGES = [
    'seqr-loader',
    'click',
    'google',
    'slackclient',
    'fsspec',
    'sklearn',
    'gcloud',
]


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


def can_reuse(fpath: str, overwrite: bool) -> bool:
    """
    Checks if the file `fpath` exists and we are not overwriting
    """
    if not fpath:
        return False
    if not file_exists(fpath):
        return False
    elif overwrite:
        logger.info(f'File {fpath} exists and will be overwritten')
        return False
    else:
        logger.info(f'Reusing existing {fpath}. Use --overwrite to overwrite')
        return True


def hash_sample_ids(sample_names: Iterable[str]) -> str:
    """
    Return a unique hash string from a set of strings
    :param sample_names: set of strings
    :return: a string hash
    """
    for sn in sample_names:
        assert ' ' not in sn, sn
    return hashlib.sha256(' '.join(sorted(sample_names)).encode()).hexdigest()[:32]


def get_refs(b: hb.Batch) -> Tuple:
    """
    Register reference files
    :param b: batch object
    :return: a tuple of reference objects
    """
    reference = b.read_input_group(
        base=REF_FASTA,
        fai=REF_FASTA + '.fai',
        dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
        + '.dict',
    )
    bwa_reference = b.read_input_group(
        base=REF_FASTA,
        fai=REF_FASTA + '.fai',
        dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
        + '.dict',
        amb=REF_FASTA + '.amb',
        ann=REF_FASTA + '.ann',
        pac=REF_FASTA + '.pac',
        o123=REF_FASTA + '.0123',
        bwa2bit64=REF_FASTA + '.bwt.2bit.64',
    )
    noalt_regions = b.read_input(NOALT_REGIONS)
    return reference, bwa_reference, noalt_regions


def make_sm_in_progress_job(*args, **kwargs) -> Job:
    """
    Creates a job that updates the sample metadata server entry analysis status
    to "in-progress"
    """
    kwargs['status'] = 'in-progress'
    return make_sm_update_status_job(*args, **kwargs)


def make_sm_completed_job(*args, **kwargs) -> Job:
    """
    Creates a job that updates the sample metadata server entry analysis status
    to "completed"
    """
    kwargs['status'] = 'completed'
    return make_sm_update_status_job(*args, **kwargs)


def make_sm_update_status_job(
    b: Batch,
    project: str,
    analysis_id: str,
    analysis_type: str,
    status: str,
    sample_name: Optional[str] = None,
    project_name: Optional[str] = None,
) -> Job:
    """
    Creates a job that updates the sample metadata server entry analysis status.
    """
    assert status in ['in-progress', 'failed', 'completed', 'queued']
    job_name = ''
    if project_name and sample_name:
        job_name += f'{project_name}/{sample_name}: '
    job_name += f'Update SM: {analysis_type} to {status}'
    j = b.new_job(job_name)
    j.image(SM_IMAGE)
    j.command(
        f"""
set -o pipefail
set -ex

export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS
export SM_DEV_DB_PROJECT={project}
export SM_ENVIRONMENT=PRODUCTION

cat <<EOT >> update.py
from sample_metadata.api import AnalysisApi
from sample_metadata import AnalysisUpdateModel
aapi = AnalysisApi()
aapi.update_analysis_status(
    analysis_id='{analysis_id}',
    analysis_update_model=AnalysisUpdateModel(status='{status}'),
)
EOT
python update.py
    """
    )
    return j


def replace_paths_to_test(s: Dict) -> Optional[Dict]:
    """
    Replace paths of all files in -main namespace to -test namespsace,
    and return None if files in -test are not found.
    :param s:
    :return:
    """

    def fix(fpath):
        fpath = fpath.replace('-main-upload/', '-test-upload/')
        if not file_exists(fpath):
            return None
        return fpath

    try:
        reads_type = s['meta']['reads_type']
        if reads_type in ('bam', 'cram'):
            fpath = s['meta']['reads'][0]['location']
            fpath = fix(fpath)
            if not fpath:
                return None
            s['meta']['reads'][0]['location'] = fpath

            fpath = s['meta']['reads'][0]['secondaryFiles'][0]['location']
            fpath = fix(fpath)
            if not fpath:
                return None
            s['meta']['reads'][0]['secondaryFiles'][0]['location'] = fpath

        elif reads_type == 'fastq':
            for li in range(len(s['meta']['reads'])):
                for rj in range(len(s['meta']['reads'][li])):
                    fpath = s['meta']['reads'][li][rj]['location']
                    fpath = fix(fpath)
                    if not fpath:
                        return None
                    s['meta']['reads'][li][rj]['location'] = fpath

        logger.info(f'Found test sample {s["id"]}')
        return s
    except Exception:  # pylint: disable=broad-except
        return None


def gsutil_cp(
    src_path: str,
    dst_path: str,
    disable_check_hashes: bool = False,
    recursive: bool = False,
    quiet: bool = False,
):
    """
    Wrapper around `gsutil cp`

    :param src_path: path to a file to copy from
    :param dst_path: path to copy to
    :param disable_check_hashes:
        Uses the gsutil option `-o GSUtil:check_hashes=never` which is required to
        get around the gsutil integrity checking error, as conda gsutil doesn't use
        CRC32c:
        > Downloading this composite object requires integrity checking with CRC32c,
          but your crcmod installation isn't using the module's C extension, so the
          hash computation will likely throttle download performance.

          To download regardless of crcmod performance or to skip slow integrity
          checks, see the "check_hashes" option in your boto config file.
    :param recursive: to copy a directory
    :param quiet: disable logging of commands and copied files
    """
    cmd = (
        'gsutil '
        + ('-q ' if quiet else '')
        + ('-o GSUtil:check_hashes=never ' if disable_check_hashes else '')
        + 'cp '
        + ('-r ' if recursive else '')
        + f'{src_path} {dst_path}'
    )
    if not quiet:
        logger.info(cmd)
    subprocess.run(cmd, check=False, shell=True)
