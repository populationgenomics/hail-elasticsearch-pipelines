"""
Utility functions for the seqr loader
"""

import os
from os.path import join

from google.cloud import storage


GATK_VERSION = '4.2.0.0'
GATK_CONTAINER = (
    f'australia-southeast1-docker.pkg.dev/cpg-common/images/gatk:{GATK_VERSION}'
)
PICARD_CONTAINER = f'us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8'
BAZAM_CONTAINER = f'australia-southeast1-docker.pkg.dev/cpg-common/images/bazam:v2'
SOMALIER_CONTAINER = 'brentp/somalier:latest'
PEDDY_CONTAINER = 'quay.io/biocontainers/peddy:0.4.8--pyh5e36f6f_0'
GNARLY_DOCKER = 'australia-southeast1-docker.pkg.dev/cpg-common/images/gnarly_genotyper:hail_ukbb_300K'
BCFTOOLS_DOCKER = (
    'australia-southeast1-docker.pkg.dev/cpg-common/images/bcftools:1.10.2--h4f4756c_2'
)

NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS = 50
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
