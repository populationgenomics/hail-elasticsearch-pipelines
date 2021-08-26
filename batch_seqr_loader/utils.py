"""
Utility functions and constants for the seqr loader pipeline
"""

import hashlib
import os
from os.path import join
from typing import Iterable, Tuple, Optional
import logging
from google.cloud import storage
import hailtop.batch as hb
from hailtop.batch import Batch
from hailtop.batch.job import Job

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
GATK_VERSION = '4.2.1.0'
GATK_IMAGE = f'{AR_REPO}/gatk:{GATK_VERSION}'
PICARD_IMAGE = f'{AR_REPO}/picard-cloud:2.23.8'
BAZAM_IMAGE = f'{AR_REPO}/bazam:v2'
SOMALIER_IMAGE = f'{AR_REPO}/somalier:latest'
PEDDY_IMAGE = f'{AR_REPO}/peddy:0.4.8--pyh5e36f6f_0'
GNARLY_IMAGE = f'{AR_REPO}/gnarly_genotyper:hail_ukbb_300K'
BCFTOOLS_IMAGE = f'{AR_REPO}/bcftools:1.10.2--h4f4756c_2'
SM_IMAGE = f'{AR_REPO}/sm-api:2.0.5'

NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS = 50
NUMBER_OF_GENOMICS_DB_INTERVALS = 50
NUMBER_OF_DATAPROC_WORKERS = 50

REF_BUCKET = 'gs://cpg-reference/hg38/v1'
REF_FASTA = join(REF_BUCKET, 'Homo_sapiens_assembly38.fasta')
DBSNP_VCF = join(REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf')
UNPADDED_INTERVALS = join(REF_BUCKET, 'hg38.even.handcurated.20k.intervals')
SOMALIER_SITES = join(REF_BUCKET, 'sites.hg38.vcf.gz')

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
        sa=REF_FASTA + '.sa',
        amb=REF_FASTA + '.amb',
        bwt=REF_FASTA + '.bwt',
        ann=REF_FASTA + '.ann',
        pac=REF_FASTA + '.pac',
    )
    noalt_regions = b.read_input(join(REF_BUCKET, 'noalt.bed'))
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
