import hail as hl
from analysis_runner import bucket_path, output_path

hl.init(default_reference='GRCh38')

mt = hl.read_matrix_table(bucket_path('seqr-loader-tmp/zornitza-stark-kccg-gvcf/v1-0/jointly-called/680e476abb7a1e4f9ba35a97d38d273cb4e68ba1390a51bf3ce1b893/annotated.mt', 'tmp'))

# See https://github.com/populationgenomics/hail-elasticsearch-pipelines/blob/91f6adf07d4708391340dd2f568fca7d166c5735/luigi_pipeline/lib/model/seqr_mt_schema.py#L304.
ht = mt.rows()
ht = ht.drop(ht.vep)
ht = ht.flatten()
ht = ht.drop(ht.locus, ht.alleles)

df = ht.to_spark()
df.write.parquet(output_path('table.parquet', 'tmp'))
