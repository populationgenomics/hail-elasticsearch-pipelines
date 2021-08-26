"""Converts from Parquet to Arrow IPC."""

import click
import math
import google.cloud.storage as gcs
import pyarrow as pa
import pyarrow.parquet as pq


@click.command()
@click.option('--input', help='Input path for Parquet files', required=True)
@click.option('--output', help='Output path for Arrow files', required=True)
@click.option(
    '--shard_index', help='Shard index for input files', type=int, required=True
)
@click.option(
    '--shard_count', help='Shard count for input files', type=int, required=True
)
def parquet_to_arrow(input, output, shard_index, shard_count):
    gcs_client = gcs.Client()

    def bucket_and_name(gcs_path):
        parts = gcs_path.split('/')
        return parts[2], '/'.join(parts[3:])

    input_bucket, input_dir = bucket_and_name(input)
    output_bucket, output_dir = bucket_and_name(output)

    all_blobs = list(gcs_client.list_blobs(input_bucket, prefix=input_dir))
    parquet_blobs = [blob for blob in all_blobs if blob.name.endswith('.parquet')]
    num_blobs = len(parquet_blobs)
    per_shard = math.ceil(num_blobs / shard_count)
    shard_blobs = parquet_blobs[shard_index * per_shard : (shard_index + 1) * per_shard]
    for input_blob in shard_blobs:
        print(f'Reading {input_blob.name}...')
        bytes = input_blob.download_as_bytes()
        buffer_reader = pa.BufferReader(bytes)
        pq_file = pq.ParquetFile(buffer_reader)
        table = pq_file.read()

        print('Converting to Arrow format...')
        output_buffer_stream = pa.BufferOutputStream()
        ipc_options = pa.ipc.IpcWriteOptions(compression='zstd')
        with pa.ipc.RecordBatchFileWriter(
            output_buffer_stream, table.schema, options=ipc_options
        ) as ipc_writer:
            ipc_writer.write_table(table)
        buffer = output_buffer_stream.getvalue()

        output_name = input_blob.name.split('/')[-1].replace('.parquet', '.arrow')
        output_blob = gcs.Blob(f'{output_dir}/{output_name}', output_bucket)
        print(f'Writing {output_blob.name}...')
        output_blob.upload_from_string(buffer)


if __name__ == '__main__':
    parquet_to_arrow()  # pylint: disable=no-value-for-parameter
