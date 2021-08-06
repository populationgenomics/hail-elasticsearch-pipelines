"""
Utility functions for the seqr loader
"""
import hashlib
import os
from os.path import join
from typing import Iterable, Tuple
import logging
from google.cloud import storage
import hailtop.batch as hb


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


GATK_VERSION = '4.2.0.0'
GATK_IMAGE = (
    f'australia-southeast1-docker.pkg.dev/cpg-common/images/gatk:{GATK_VERSION}'
)
PICARD_IMAGE = f'us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8'
BAZAM_IMAGE = f'australia-southeast1-docker.pkg.dev/cpg-common/images/bazam:v2'
SOMALIER_IMAGE = 'brentp/somalier:latest'
PEDDY_IMAGE = 'quay.io/biocontainers/peddy:0.4.8--pyh5e36f6f_0'
GNARLY_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/gnarly_genotyper:hail_ukbb_300K'
BCFTOOLS_IMAGE = (
    'australia-southeast1-docker.pkg.dev/cpg-common/images/bcftools:1.10.2--h4f4756c_2'
)

NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS = 10
NUMBER_OF_DATAPROC_WORKERS = 50
NUMBER_OF_GENOMICS_DB_INTERVALS = 10

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
    if not file_exists(fpath):
        return False
    elif overwrite:
        logger.info(f'File {fpath} exists and will be overwritten')
        return False
    else:
        logger.info(f'Reusing existing {fpath}. Use --overwrite to overwrite')
        return True


def hash_sample_names(sample_names: Iterable[str]) -> str:
    """
    Return a unique hash string from a from a set of strings
    :param sample_names: set of strings
    :return: a string hash
    """
    for sn in sample_names:
        assert ' ' not in sn, sn
    return hashlib.sha224(' '.join(sorted(sample_names)).encode()).hexdigest()


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
