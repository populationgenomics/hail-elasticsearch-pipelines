#!/usr/bin/env python3

"""
Hail script to submit on a dataproc cluster. 

Annotates the input matrix table with SeqrGenotypesSchema, and loads into ES.
"""

import logging
import math
import subprocess
from os.path import join

import click
import hail as hl
from lib.model.seqr_mt_schema import SeqrGenotypesSchema
from hail_scripts.v02.utils.elasticsearch_client import ElasticsearchClient

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.option(
    '--mt-path',
    'mt_path',
    required=True,
)
@click.option(
    '--bucket',
    'work_bucket',
    required=True,
)
@click.option(
    '--es-host',
    'es_host',
    help='Elasticsearch host',
)
@click.option('--es-port', 'es_port', type=click.STRING, help='Elasticsearch port')
@click.option(
    '--es-username', 'es_username', type=click.STRING, help='Elasticsearch username'
)
@click.option(
    '--es-password', 'es_password', type=click.STRING, help='Elasticsearch password'
)
@click.option(
    '--es-index',
    'es_index',
    type=click.STRING,
    help='Elasticsearch index. Usually the dataset name. Will be lowercased',
    required=True,
)
@click.option(
    '--es-index-min-num-shards',
    'es_index_min_num_shards',
    default=1,
    help='Number of shards for the index will be the greater of this value '
    'and a calculated value based on the matrix.',
)
@click.option(
    '--sample-list',
    'sample_list_fpath',
    help='File with a comma-separated list of samples to subset the input matrix table to',
)
@click.option(
    '--prod', is_flag=True, help='Run under the production ES credentials instead'
)
@click.option('--genome-version', 'genome_version', default='GRCh38')
@click.option(
    '--make-checkpoints',
    'make_checkpoints',
    is_flag=True,
    help='Create checkpoints for intermediate matrix tables',
)
def main(
    mt_path: str,
    work_bucket: str,
    es_host: str,
    es_port: str,
    es_username: str,
    es_password: str,
    es_index: str,
    es_index_min_num_shards: int,
    sample_list_fpath: str,
    prod: bool,  # pylint: disable=unused-argument
    genome_version: str,
    make_checkpoints: bool = False,
):  # pylint: disable=missing-function-docstring
    hl.init(default_reference=genome_version)

    if not all([es_host, es_port, es_username, es_password]):
        if any([es_host, es_port, es_username, es_password]):
            raise click.BadParameter(
                f'Either none, or all ES configuration parameters '
                f'must be specified: --es-host, --es-port, --es-username, --es-password. '
                f'If none are specified, defaults for the CPG are used'
            )
        es_host = 'elasticsearch.es.australia-southeast1.gcp.elastic-cloud.com'
        es_port = '9243'
        es_username = 'seqr'
        es_password = _read_es_password()

    es = ElasticsearchClient(
        host=es_host,
        port=str(es_port),
        es_username=es_username,
        es_password=es_password,
        es_use_ssl=(es_host != 'localhost'),
    )

    mt = hl.read_matrix_table(mt_path)

    # Subsetting to requested samples
    with hl.hadoop_open(sample_list_fpath, 'r') as f:
        sample_list = f.read().strip().split(',')
    mt = mt.filter_cols(hl.set(sample_list).contains(mt.s))

    logger.info('Annotating genotypes')
    mt = _compute_genotype_annotated_vcf(mt)
    if make_checkpoints:
        mt = mt.checkpoint(join(work_bucket, 'annotated_genotype.mt'), overwrite=True)

    row_table = SeqrGenotypesESSchema.elasticsearch_row(mt)

    es_shards = _mt_num_shards(mt, es_index_min_num_shards)
    logger.info('Exporting to ES')
    es.export_table_to_elasticsearch(
        row_table,
        index_name=es_index.lower(),
        num_shards=es_shards,
        write_null_values=True,
    )
    _cleanup(es, es_index, es_shards)


def _read_es_password(
    project_id='seqr-308602',
    secret_id='seqr-es-password',
    version_id='latest',
) -> str:
    """
    Read a GCP secret storing the ES password
    """
    cmd = f'gcloud secrets versions access {version_id} --secret {secret_id} --project {project_id}'
    logger.info(cmd)
    return subprocess.check_output(cmd, shell=True).decode()


def _mt_num_shards(mt, es_index_min_num_shards):
    """
    The greater of the user specified min shards and calculated based on the variants
    and samples
    """
    denominator = 1.4 * 10 ** 9
    calculated_num_shards = math.ceil((mt.count_rows() * mt.count_cols()) / denominator)
    return max(es_index_min_num_shards, calculated_num_shards)


def _cleanup(es, es_index, es_shards):
    es.route_index_off_temp_es_cluster(es_index)
    # Current disk configuration requires the previous index to be deleted prior to large indices, ~1TB, transferring off loading nodes
    if es_shards < 25:
        es.wait_for_shard_transfer(es_index)


class SeqrGenotypesESSchema(SeqrGenotypesSchema):
    """
    Modified version of SeqrVariantsAndGenotypesSchema to just base
    on SeqrGenotypesSchema, as we apply SeqrVariantSchema separately
    on the cohort level
    """

    @staticmethod
    def elasticsearch_row(ds):
        """
        Prepares the mt to export using ElasticsearchClient V02.
        - Flattens nested structs
        - drops locus and alleles key
        """
        # Converts a mt to the row equivalent.
        if isinstance(ds, hl.MatrixTable):
            ds = ds.rows()
        # Converts nested structs into one field, e.g. {a: {b: 1}} => a.b: 1
        table = ds.drop('vep').flatten()
        # When flattening, the table is unkeyed, which causes problems because our locus and alleles should not
        # be normal fields. We can also re-key, but I believe this is computational?
        table = table.drop(table.locus, table.alleles)
        return table


def _compute_genotype_annotated_vcf(
    mt, schema_cls=SeqrGenotypesESSchema
) -> hl.MatrixTable:
    r"""
                   BaseMTSchema
                     /       \
            SeqrSchema        |
                    |         |
      SeqrVariantSchema     SeqrGenotypesSchema
                    |         |
      SeqrVariantASSchema   SeqrGenotypesESSchema
                    \        /
            SeqrVariantsAndGenotypesSchema
            
    SeqrVariantASSchema is applied on the cohort level separately.
    """
    annotation_schema = schema_cls(mt)
    mt = annotation_schema.annotate_all(overwrite=True).mt
    return mt


if __name__ == '__main__':
    main()  # pylint: disable=E1120
