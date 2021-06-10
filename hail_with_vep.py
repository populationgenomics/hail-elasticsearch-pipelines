import hail as hl

hl.init(default_reference='GRCh38')

import subprocess

cmd = 'ls /'
print(f'Running {cmd}')
subprocess.run(cmd, shell=True, check=False)
print()

cmd = 'ls /vep_data'
print(f'Running {cmd}')
subprocess.run(cmd, shell=True, check=False)
print()

vcf_path = 'gs://playground-au/seqr/vcf/seqr.vcf.gz'
mt = hl.import_vcf(
    [vcf_path],
    reference_genome='GRCh38',
    skip_invalid_loci=True,
    force_bgz=True,
    min_partitions=100,
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


mt = annotate_old_and_split_multi_hts(mt)

mt = hl.vep(mt, block_size=1000, config='/vep_data/vep-gcloud.json')
mt.entries().show()
