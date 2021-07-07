import click
import hail as hl

@click.command()
@click.option('--input_path', help='GCS path to the annotated Hail MT to read', required=True)
@click.option('--output_path', help='GCS path to the Parquet file to write', required=True)
def main(input_path, output_path):
    hl.init(default_reference='GRCh38')
    mt = hl.read_matrix_table(input_path)
    # See https://github.com/populationgenomics/hail-elasticsearch-pipelines/blob/91f6adf07d4708391340dd2f568fca7d166c5735/luigi_pipeline/lib/model/seqr_mt_schema.py#L304.
    ht = mt.rows()
    ht = ht.drop(ht.vep)
    ht = ht.flatten()
    ht = ht.drop(ht.locus, ht.alleles)
    df = ht.to_spark()
    df.write.parquet(output_path)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
