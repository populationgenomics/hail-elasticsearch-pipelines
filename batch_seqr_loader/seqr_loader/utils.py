"""
Utility functions and constants for the seqr loader pipeline
"""

import hashlib
import os
import subprocess
from dataclasses import dataclass
from os.path import join
from typing import Iterable, Tuple, Optional, List
import logging
from google.cloud import storage

import hailtop.batch as hb

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
GATK_VERSION = '4.2.2.0-cpgfix00'  # Our fork with a couple of fixes:
# https://github.com/populationgenomics/production-pipelines/tree/initial/dockers/gatk
GATK_IMAGE = f'{AR_REPO}/gatk:{GATK_VERSION}'
PICARD_IMAGE = f'{AR_REPO}/picard_samtools:v0'
ALIGNMENT_IMAGE = f'{AR_REPO}/alignment:v3'
SOMALIER_IMAGE = f'{AR_REPO}/somalier:latest'
PEDDY_IMAGE = f'{AR_REPO}/peddy:0.4.8--pyh5e36f6f_0'
GNARLY_IMAGE = f'{AR_REPO}/gnarly_genotyper:hail_ukbb_300K'
BCFTOOLS_IMAGE = f'{AR_REPO}/bcftools:1.10.2--h4f4756c_2'
SM_IMAGE = f'{AR_REPO}/sm-api:3.3.0'

NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS = 50
NUMBER_OF_GENOMICS_DB_INTERVALS = 50
NUMBER_OF_DATAPROC_WORKERS = 50

REF_BUCKET = 'gs://cpg-reference'
NOALT_REGIONS = join(REF_BUCKET, 'noalt.bed')
SOMALIER_SITES = join(REF_BUCKET, 'somalier/v0/sites.hg38.vcf.gz')

BROAD_REF_BUCKET = f'{REF_BUCKET}/hg38/v1'
REF_FASTA = join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.fasta')
DBSNP_VCF = join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf')
UNPADDED_INTERVALS = join(BROAD_REF_BUCKET, 'hg38.even.handcurated.20k.intervals')

DATAPROC_PACKAGES = [
    'seqr-loader',
    'elasticsearch',
    'cpg-gnomad',
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
        ann=REF_FASTA + '.ann',
        pac=REF_FASTA + '.pac',
        bwt=REF_FASTA + '.bwt',
    )
    noalt_regions = b.read_input(NOALT_REGIONS)
    return reference, bwa_reference, noalt_regions


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


@dataclass
class AlignmentInput:
    """
    Sort of a union type for possible alignment inputs
    """

    bam_or_cram_path: Optional[str] = None
    index_path: Optional[str] = None
    fqs1: Optional[List[str]] = None
    fqs2: Optional[List[str]] = None
