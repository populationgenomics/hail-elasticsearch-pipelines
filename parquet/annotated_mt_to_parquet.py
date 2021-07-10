import hail as hl
from analysis_runner import bucket_path, output_path

hl.init(default_reference='GRCh38')

mt = hl.read_matrix_table(
    bucket_path(
        'seqr-loader-tmp/acute-care/v1-0/jointly-called/680e476abb7a1e4f9ba35a97d38d273cb4e68ba1390a51bf3ce1b893/annotated.mt',
        'tmp',
    )
)


def maybe_call(f, filename):
    """Calls f(filename) iff filename doesn't exist yet."""
    if not hl.hadoop_exists(filename):
        f(filename)


# Discard VEP, which was only used to derive other annotations.
mt = mt.drop(mt.vep)
maybe_call(mt.write, output_path('no_vep.mt', 'tmp'))

# See https://github.com/populationgenomics/hail-elasticsearch-pipelines/blob/91f6adf07d4708391340dd2f568fca7d166c5735/luigi_pipeline/lib/model/seqr_mt_schema.py#L304.
ht = mt.rows().flatten()
ht = ht.drop(ht.locus, ht.alleles)
maybe_call(ht.write, output_path('no_vep.ht', 'tmp'))

df = ht.to_spark()
maybe_call(
    df.write.option('compression', 'none').parquet,
    output_path('original_uncompressed.parquet', 'tmp'),
)
maybe_call(
    df.write.option('compression', 'zstd').parquet,
    output_path('original.parquet', 'tmp'),
)

# Drop all sample-specific fields.
mt = mt.drop(
    mt.genotypes, mt.samples_no_call, mt.samples_num_alt, mt.samples_gq, mt.samples_ab
)
ht = mt.rows().flatten()
ht = ht.drop(ht.locus, ht.alleles)
df = ht.to_spark()
maybe_call(
    df.write.option('compression', 'zstd').parquet,
    output_path('sites_only.parquet', 'tmp'),
)

# Annotate entries with genotype information instead.
mt = mt.select_entries(
    num_alt=mt.GT.n_alt_alleles(),
    gq=mt.GQ,
    ab=hl.bind(
        lambda total: hl.or_missing(
            (total > 0) & (hl.len(mt.AD) > 1), mt.AD[1] / total
        ),
        hl.sum(mt.AD),
    ),
    dp=mt.DP,
)
ht = mt.make_table()
ht = ht.key_by()
ht = ht.drop(ht.locus, ht.alleles)
df = ht.to_spark()
maybe_call(
    df.write.option('compression', 'zstd').parquet,
    output_path('genotype_columns.parquet', 'tmp'),
)
