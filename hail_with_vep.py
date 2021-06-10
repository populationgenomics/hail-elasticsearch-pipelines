import hail as hl
hl.init(default_reference='GRCh38')

import subprocess
res1 = subprocess.check_output('ls /', shell=True)
res2 = subprocess.check_output('ls /vep_data', shell=True)
print(res1.decode().split('\n'), res2.decode().split('\n'))

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

mt = hl.vep(mt, block_size=1000)
mt.entries().show()