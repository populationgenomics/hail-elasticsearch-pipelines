"""Converts from Hail MT to Parquet."""

import click
import hail as hl


@click.command()
@click.option(
    '--input', help='Input path for annotated Hail MatrixTable', required=True
)
@click.option('--output', help='Output path for Parquet partitions', required=True)
def mt_to_parquet(input, output):
    """Script entry point."""

    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(input)

    # Convert MT to HT, dropping VEP. See `elasticsearch_row` in the ES-based pipeline.
    ht = mt.rows()
    ht = ht.drop(ht.vep)
    ht = ht.flatten()
    ht = ht.drop(ht.locus, ht.alleles)

    # Convert HT to Parquet.
    df = ht.to_spark()
    df.write.parquet(output)


if __name__ == '__main__':
    mt_to_parquet()  # pylint: disable=no-value-for-parameter
