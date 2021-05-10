#!/usr/bin/env python3

import logging
from typing import Optional, List

import click
import hail as hl

from hail_scripts.v02.utils.elasticsearch_client import ElasticsearchClient
from model.seqr_mt_schema import SeqrVariantsAndGenotypesSchema

logger = logging.getLogger()


@click.command()
@click.option(
    '--source-path', 'source_paths', multiple=True, required=True,
)
@click.option(
    '--dest-mt-path', 'dest_path', required=True,
)
@click.option(
    '--bucket', 'bucket', required=True,
)
@click.option('--genome-version', 'genome_version', default='GRCh38')
@click.option(
    '--reference-ht', 'reference_path', required=True,
    help='Path to the Hail table storing the reference variants.',
    default='gs://seqr-reference-data/GRCh38/all_reference_data/v2/combined_reference_data_grch38-2.0.2.ht',
)
@click.option(
    '--clinvar-ht', 'clinvar_path', required=True,
    help='Path to the Hail table storing the clinvar variants.',
    default='gs://seqr-reference-data/GRCh38/clinvar/clinvar.GRCh38.ht',
)
@click.option(
    '--hgmd-ht', 'hgmd_path',
    help='Path to the Hail table storing the hgmd variants.',
    # default='gs://seqr-reference-data-private/GRCh38/HGMD/hgmd_pro_2018.4_hg38_without_db_field.vds'
)
@click.option(
    '--disable-validation', 'disable_validation', is_flag=True
)
@click.option(
    '--sample-type', 'sample_type', type=click.Choice(['WGS', 'WES']), default='WGS'
)
@click.option(
    '--dataset-type', 'dataset_type', type=click.Choice(['VARIANTS', 'SV']), 
    default='VARIANTS'
)
@click.option(
    '--remap-tsv', 'remap_path',
    help='Path to a TSV file with two columns: s and seqr_id.'
)
@click.option(
    '--subset-tsv', 'subset_path',
    help='Path to a TSV file with one column of sample IDs: s.'
)
@click.option(
    '--vep-config', 'vep_config_json_path',
    help='Path of hail vep config .json file',
    #default='gs://cpg-reference/hg38/v0/vep-config.json',
)
@click.option(
    '--vep-block-size', 'vep_block_size'
)
def main(
    source_paths: List[str],
    dest_path: str,
    bucket: str,
    genome_version: str,
    reference_path: str,
    clinvar_path: str,
    hgmd_path: Optional[str],
    disable_validation: bool,
    sample_type: str,
    dataset_type: str,  # pylint: disable=unused-argument
    remap_path: str = None,  # pylint: disable=unused-argument
    subset_path: str = None,  # pylint: disable=unused-argument
    vep_config_json_path: Optional[str] = None,
    vep_block_size: Optional[int] = None,
):
    print('Running')
    genome_version = genome_version.replace('GRCh', '')
    mt = import_vcf(source_paths, genome_version)
    mt = annotate_old_and_split_multi_hts(mt)
    mt.write(f'{bucket}/annotate_old_and_split_multi_hts.mt', overwrite=True)
    if not disable_validation:
        validate_mt(mt, sample_type)
    if remap_path:
        mt = remap_sample_ids(mt, remap_path)
    if subset_path:
        mt = subset_samples_and_variants(mt, subset_path)
    if genome_version == '38':
        mt = add_37_coordinates(mt)
        mt.write(f'{bucket}/add_37_coordinates.mt', overwrite=True)
    mt = run_vep(
        mt,
        genome_version=genome_version,
        vep_config_json_path=vep_config_json_path,
        block_size=vep_block_size,
    )
    mt.write(f'{bucket}/run_vep.mt', overwrite=True)

    ref_data = hl.read_table(reference_path)
    clinvar = hl.read_table(clinvar_path)
    hgmd = hl.read_table(hgmd_path) if hgmd_path else None

    mt = compute_annotated_vcf(mt, ref_data=ref_data, clinvar=clinvar, hgmd=hgmd)
    mt.write(f'{bucket}/compute_annotated_vcf.mt', overwrite=True)

    mt = mt.annotate_globals(
        sourceFilePath=",".join(source_paths),
        genomeVersion=genome_version,
        sampleType=sample_type,
        hail_version=hl.version(),
    )

    mt.describe()
    mt.write(dest_path, stage_locally=True, overwrite=True)


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


def run_vep(
    mt: hl.MatrixTable,
    genome_version: str,
    name: str = "vep",
    vep_config_json_path: str = None,
    block_size: Optional[int] = None,
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

        # Default config when you start the VEP cluster on spark
        config = "file:///vep_data/vep-gcloud.json"

    mt = hl.vep(mt, config=config, name=name, block_size=block_size or 1000)

    logger.info("==> Done with VEP")
    return mt


def get_sample_type_stats(mt, threshold=0.3):
    """
    Calculate stats for sample type by checking against a list of common coding and non-coding variants.
    If the match for each respective type is over the threshold, we return a match.

    :param mt: Matrix Table to check
    :param threshold: if the matched percentage is over this threshold, we classify as match
    :return: a dict of coding/non-coding to dict with 'matched_count', 'total_count' and 'match' boolean.
    """
    stats = {}
    types_to_ht_path = {
        'noncoding': 'gs://seqr-reference-data/GRCh38/validate_ht/common_noncoding_variants.grch38.ht',
        'coding': 'gs://seqr-reference-data/GRCh38/validate_ht/common_coding_variants.grch38.ht'
    }
    for sample_type, ht_path in types_to_ht_path.items():
        ht = hl.read_table(ht_path)
        stats[sample_type] = ht_stats = {
            'matched_count': mt.semi_join_rows(ht).count_rows(),
            'total_count': ht.count(),

        }
        ht_stats['match'] = (ht_stats['matched_count']/ht_stats['total_count']) >= threshold
    return stats


class SeqrValidationError(Exception):
    pass


def validate_mt(mt, sample_type):
    """
    Validate the mt by checking against a list of common coding and non-coding variants given its
    genome version. This validates genome_version, variants, and the reported sample type.

    :param mt: mt to validate
    :param genome_version: reference genome version
    :param sample_type: WGS or WES
    :return: True or Exception
    """
    sample_type_stats = get_sample_type_stats(mt)

    for name, stat in sample_type_stats.items():
        logger.info('Table contains %i out of %i common %s variants.' %
                    (stat['matched_count'], stat['total_count'], name))

    has_coding = sample_type_stats['coding']['match']
    has_noncoding = sample_type_stats['noncoding']['match']

    if not has_coding and not has_noncoding:
        # No common variants detected.
        raise SeqrValidationError(
            'Genome version validation error: dataset specified but doesn\'t contain the expected number of common variants')
    elif has_noncoding and not has_coding:
        # Non coding only.
        raise SeqrValidationError(
            'Sample type validation error: Dataset contains noncoding variants but is missing common coding '
            'variants. Please verify that the dataset contains coding variants.'
        )
    elif has_coding and not has_noncoding:
        # Only coding should be WES.
        if sample_type != 'WES':
            raise SeqrValidationError(
                f'Sample type validation error: dataset sample-type is specified as {sample_type} but appears to be '
                'WGS because it contains many common coding variants'
            )
    elif has_noncoding and has_coding:
        # Both should be WGS.
        if sample_type != 'WGS':
            raise SeqrValidationError(
                f'Sample type validation error: dataset sample-type is specified as {sample_type} but appears to be '
                'WES because it contains many common non-coding variants'
            )
    return True


def remap_sample_ids(mt, remap_path):
    """
    Remap the MatrixTable's sample ID, 's', field to the sample ID used within seqr, 'seqr_id'
    If the sample 's' does not have a 'seqr_id' in the remap file, 's' becomes 'seqr_id'
    :param mt: MatrixTable from VCF
    :param remap_path: Path to a file with two columns 's' and 'seqr_id'
    :return: MatrixTable remapped and keyed to use seqr_id
    """
    remap_ht = hl.import_table(remap_path, key='s')
    missing_samples = remap_ht.anti_join(mt.cols()).collect()
    remap_count = remap_ht.count()

    if len(missing_samples) != 0:
        raise Exception(
            f'Only {remap_ht.semi_join(mt.cols()).count()} out of {remap_count} '
            'remap IDs matched IDs in the variant callset.\n'
            f'IDs that aren\'t in the callset: {missing_samples}\n'
            f'All callset sample IDs:{mt.s.collect()}', missing_samples
        )

    mt = mt.annotate_cols(**remap_ht[mt.s])
    remap_expr = hl.cond(hl.is_missing(mt.seqr_id), mt.s, mt.seqr_id)
    mt = mt.annotate_cols(seqr_id=remap_expr, vcf_id=mt.s)
    mt = mt.key_cols_by(s=mt.seqr_id)
    logger.info(f'Remapped {remap_count} sample ids...')
    return mt


def subset_samples_and_variants(mt, subset_path):
    """
    Subset the MatrixTable to the provided list of samples and to variants present in those samples
    :param mt: MatrixTable from VCF
    :param subset_path: Path to a file with a single column 's'
    :return: MatrixTable subsetted to list of samples
    """
    subset_ht = hl.import_table(subset_path, key='s')
    subset_count = subset_ht.count()
    anti_join_ht = subset_ht.anti_join(mt.cols())
    anti_join_ht_count = anti_join_ht.count()

    if anti_join_ht_count != 0:
        missing_samples = anti_join_ht.s.collect()
        raise Exception(
            f'Only {subset_count-anti_join_ht_count} out of {subset_count} '
            'subsetting-table IDs matched IDs in the variant callset.\n'
            f'IDs that aren\'t in the callset: {missing_samples}\n'
            f'All callset sample IDs:{mt.s.collect()}', missing_samples
        )

    mt = mt.semi_join_cols(subset_ht)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    logger.info(f'Finished subsetting samples. Kept {subset_count} '
                f'out of {mt.count()} samples in vds')
    return mt


def add_37_coordinates(mt):
    """Annotates the GRCh38 MT with 37 coordinates using hail's built-in liftover
    :param mt: MatrixTable from VCF
    :return: MatrixTable annotated with GRCh37 coordinates
    """
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg38.add_liftover('gs://hail-common/references/grch38_to_grch37.over.chain.gz', rg37)
    mt = mt.annotate_rows(rg37_locus=hl.liftover(mt.locus, 'GRCh37'))
    return mt


if __name__ == '__main__':
    main()  # pylint: disable=E1120
