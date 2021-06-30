#!/usr/bin/env python3

"""
Hail script to submit on a dataproc cluster. Converts input multi-sample VCFs
into a matrix table, which annotates and prepares to load into ES.
"""

import logging
import math
import click
import hail as hl

from lib.model.seqr_mt_schema import SeqrVariantsAndGenotypesSchema
from hail_scripts.v02.utils.elasticsearch_client import ElasticsearchClient

logger = logging.getLogger('load_to_es')
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.option(
    '--mt-path',
    'mt_path',
    required=True,
)
@click.option(
    '--es-host',
    'es_host',
    default='elasticsearch.es.australia-southeast1.gcp.elastic-cloud.com',
    help='ElasticSearch host',
)
@click.option('--es-port', 'es_port', default=9243, help='ElasticSearch port')
@click.option(
    '--es-use-ssl', 'es_use_ssl', is_flag=True, help='Use SSL for ElasticSearch'
)
@click.option(
    '--es-username', 'es_username', help='ElasticSearch username', default='seqr'
)
@click.option(
    '--es-index', 'es_index', help='ElasticSearch index. Usually the dataset name'
)
@click.option(
    '--es-index-min-num-shards',
    'es_index_min_num_shards',
    default=1,
    help='Number of shards for the index will be the greater of this value '
    'and a calculated value based on the matrix.',
)
@click.option('--genome-version', 'genome_version', default='GRCh38')
def main(
    mt_path: str,
    es_host: str,
    es_port: int,
    es_use_ssl: bool,
    es_username: str,
    es_index: str,
    es_index_min_num_shards: int,
    genome_version: str,
):  # pylint: disable=missing-function-docstring
    hl.init(default_reference=genome_version)

    logger.info('Starting the seqr_load pipeline')
    mt = hl.read_matrix_table(mt_path)
    row_table = SeqrVariantsAndGenotypesSchema.elasticsearch_row(mt)
    es_password = '<TODO: load from a secret?>'
    es = ElasticsearchClient(
        host=es_host,
        port=str(es_port),
        es_username=es_username,
        es_password=es_password,
        es_use_ssl=es_use_ssl,
    )
    es.export_table_to_elasticsearch(
        row_table,
        index_name=es_index,
        num_shards=_mt_num_shards(mt, es_index_min_num_shards),
        write_null_values=True,
    )


def _mt_num_shards(mt, es_index_min_num_shards):
    """
    The greater of the user specified min shards and calculated based on the variants
    and samples
    """
    denominator = 1.4 * 10 ** 9
    calculated_num_shards = math.ceil((mt.count_rows() * mt.count_cols()) / denominator)
    return max(es_index_min_num_shards, calculated_num_shards)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
