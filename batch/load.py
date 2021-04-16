import logging
from collections import defaultdict
from typing import Optional, List, Union

import hail as hl
import hailtop.batch as hb

from hail_scripts.v02.utils.elasticsearch_client import ElasticsearchClient
from luigi_pipeline.lib.model.seqr_mt_schema import SeqrVariantsAndGenotypesSchema

logger = logging.getLogger()


class ElasticSearchCredentials:
    def __init__(self, host="localhost", port=9200, index="data", username="pipeline", password=None, index_min_num_shards=6, use_ssl=None):
        self.host = host
        self.port = port
        self.index = index
        self.username = username
        self.password = password
        self.index_min_num_shards = index_min_num_shards
        self.use_ssl = use_ssl
        if self.use_ssl is None:
            self.use_ssl = host != "localhost"


def parse_and_annotate_vcfs(
    source_paths: Union[str, List[str]],
    output_path,
    genome_version,
    sample_type,
    validate,
    vep_config_json_path,
    reference_ht_path,
    clinvar_ht_path,
    hgmd_ht_path,
):
    mt = import_vcf(source_paths, genome_version)
    mt = annotate_old_and_split_multi_hts(mt)
    if validate:
        pass

    mt = run_vep(
        mt,
        genome_version=genome_version,
        vep_config_json_path=vep_config_json_path
    )

    ref_data = hl.read_table(reference_ht_path)
    clinvar = hl.read_table(clinvar_ht_path)
    hgmd = hl.read_table(hgmd_ht_path) if hgmd_ht_path else None

    mt = compute_annotated_vcf(mt, ref_data=ref_data, clinvar=clinvar, hgmd=hgmd)

    mt = mt.annotate_globals(
        sourceFilePath=",".join(source_paths),
        genomeVersion=genome_version,
        sampleType=sample_type,
        hail_version=hl.version(),
    )

    mt.describe()
    mt.write(output_path, stage_locally=True, overwrite=True)


def dump_to_estask(mt, es_credentials: ElasticSearchCredentials):
    row_table = elasticsearch_row(mt)
    es = ElasticsearchClient(
        host=es_credentials.host, port=str(es_credentials.port),
        es_username=es_credentials.username, es_password=es_credentials.password,
        es_use_ssl=es_credentials.use_ssl
    )
    es.export_table_to_elasticsearch(row_table)

    return True


def import_vcf(source_paths, genome_version):
    # https://github.com/populationgenomics/hail-elasticsearch-pipelines/blob/e41582d4842bc0d2e06b1d1e348eb071e00e01b3/luigi_pipeline/lib/hail_tasks.py#L77-L89
    # Import the VCFs from inputs. Set min partitions so that local pipeline execution takes advantage of all CPUs.
    recode = {}
    if genome_version == "38":
        recode = {f"{i}": f"chr{i}" for i in (list(range(1, 23)) + ["X", "Y"])}
    elif genome_version == "37":
        recode = {f"chr{i}": f"{i}" for i in (list(range(1, 23)) + ["X", "Y"])}

    return hl.import_vcf(
        [vcf_file for vcf_file in source_paths],
        reference_genome="GRCh" + genome_version,
        skip_invalid_loci=True,
        contig_recoding=recode,
        force_bgz=True,
        min_partitions=500,
    )


def annotate_old_and_split_multi_hts(mt):
    """
    https://github.com/populationgenomics/hail-elasticsearch-pipelines/blob/e41582d4842bc0d2e06b1d1e348eb071e00e01b3/luigi_pipeline/seqr_loading.py#L89-L96

    Saves the old allele and locus because while split_multi does this, split_multi_hts drops this. Will see if
    we can add this to split_multi_hts and then this will be deprecated.
    :return: mt that has pre-annotations
    """
    # Named `locus_old` instead of `old_locus` because split_multi_hts drops `old_locus`.
    return hl.split_multi_hts(
        mt.annotate_rows(locus_old=mt.locus, alleles_old=mt.alleles)
    )


def run_vep(
    mt: hl.MatrixTable,
    genome_version: str,
    name: str = "vep",
    block_size: int = 1000,
    vep_config_json_path=None,
) -> hl.MatrixTable:
    """Runs VEP.

    https://github.com/populationgenomics/hail-elasticsearch-pipelines/blob/e41582d4842bc0d2e06b1d1e348eb071e00e01b3/hail_scripts/v02/utils/hail_utils.py#L103-L129

    :param MatrixTable mt: MT to annotate with VEP
    :param str genome_version: "37" or "38"
    :param str name: Name for resulting row field
    :param int block_size: Number of rows to process per VEP invocation.
    :param str vep_config_json_path: JSON config path for VEP
    :return: annotated MT
    :rtype: MatrixTable
    """
    if vep_config_json_path is not None:
        config = vep_config_json_path
        mt = mt.annotate_globals(gencodeVersion="unknown")
    else:
        if genome_version not in ["37", "38"]:
            raise ValueError(f"Invalid genome version: {genome_version}")

        # I think this is a default config when you start the VEP cluster on spark
        config = "file:///vep_data/vep-gcloud.json"

    mt = hl.vep(mt, config=config, name=name, block_size=block_size)

    logger.info("==> Done with VEP")
    return mt


def compute_annotated_vcf(
    mt, ref_data, clinvar, hgmd, schema_cls=SeqrVariantsAndGenotypesSchema
):
    """
    Returns a matrix table with an annotated rows where each row annotation is a previously called
    annotation (e.g. with the corresponding method or all in case of `annotate_all`).
    :return: a matrix table
    """
    # Annotations are declared as methods on the schema_cls.
    # There's a strong coupling between the @row_annotation decorator
    # and the BaseMTSchema that's impractical to refactor, so we'll just leave it.
    #
    #   class SomeClass(BaseMTSchema):
    #       @row_annotation()
    #       def a(self):
    #           return 'a_val'
    #
    # This loops through all the @row_annotation decorated methods
    # on `schema_cls` and applies them all.
    #
    # See https://user-images.githubusercontent.com/22381693/113672350-f9b04000-96fa-11eb-91fe-e45d34c165c0.png
    # for a rough overview of the structure and methods declared on:
    #
    #               BaseMTSchema
    #                 /       \
    #        SeqrSchema     SeqrGenotypesSchema
    #                |         |
    #  SeqrVariantSchema       |
    #                \        /
    #        SeqrVariantsAndGenotypesSchema
    #
    # we can call the annotation on this class in two steps:
    annotation_schema = schema_cls(
        mt, ref_data=ref_data, clinvar_data=clinvar, hgmd_data=hgmd
    )

    mt = annotation_schema.annotate_all(overwrite=True).select_annotated_mt()

    return mt


def elasticsearch_row(ds):
    """
    Copied from: https://github.com/populationgenomics/hail-elasticsearch-pipelines/blob/e41582d4842bc0d2e06b1d1e348eb071e00e01b3/luigi_pipeline/lib/model/seqr_mt_schema.py#L269-L290


    Prepares the mt to export using ElasticsearchClient V02.
    - Flattens nested structs
    - drops locus and alleles key

    TODO:
    - Call validate
    - when all annotations are here, whitelist fields to send instead of blacklisting.
    :return:
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