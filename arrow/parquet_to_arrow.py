"""Converts from Parquet to Arrow IPC."""

import click
import math
import google.cloud.storage as gcs
import pyarrow as pa
import pyarrow.parquet as pq


@click.command()
@click.option('--input', help='Input path for Parquet partitions', required=True)
@click.option(
    '--shard_index', help='Shard index over partitions', type=int, required=True
)
@click.option(
    '--shard_count', help='Shard count over partitions', type=int, required=True
)
def parquet_to_arrow(input, shard_index, shard_count):
    gcs_client = gcs.Client()

    input_parts = input.split('/')
    input_bucket = input_parts[2]
    input_blob_prefix = '/'.join(input_parts[3:])

    all_files = list(gcs_client.list_blobs(input_bucket, prefix=input_blob_prefix))
    num_files = len(all_files)
    per_shard = math.ceil(num_files / shard_count)
    shard_files = all_files[shard_index * per_shard : (shard_index + 1) * per_shard]
    for file in shard_files:
        print(file)


if __name__ == '__main__':
    parquet_to_arrow()  # pylint: disable=no-value-for-parameter
