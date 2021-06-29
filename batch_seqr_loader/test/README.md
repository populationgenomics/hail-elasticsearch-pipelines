# Generate test data

* Clone the `fewgenomes` repository, and generate the test dataset for 2 families:

```bash
git clone https://github.com/populationgenomics/fewgenomes
pushd fewgenomes
snakemake -s prep_warp_inputs.smk -j2 -p --config families=2 input_type=gvcf dataset_name=seqr-gvcf
popd
```

* Copy the data from `fewgenomes` buckets to the `seqr` buckets

```bash
PED=fewgenomes/datasets/seqr-gvcf/samples.ped
for fam in $(cut -f1 $PED | tail -n+2 | sort -u); do
    echo $fam
    # Split the 2-families PED file into 2 separate PED files to test extending genomicsdbs
    head -n1 $PED > $fam.ped
    grep $fam $PED >> $fam.ped
    gsutil cp $fam.ped gs://cpg-seqr-test/datasets/$fam/
    # Copy GVCFs
    for sam in $(grep $fam $PED | cut -f2); do 
        echo $sam
        gsutil cp gs://cpg-fewgenomes-main/gvcf/batch0/$sam.g.vcf.gz gs://cpg-seqr-test/datasets/$fam/
        gsutil cp gs://cpg-fewgenomes-main/gvcf/batch0/$sam.g.vcf.gz.tbi gs://cpg-seqr-test/datasets/$fam/
    done
done
```

* Create a PED file version with mismatched sex and relationships:

```bash
# Make a mismatched dataset to test pedigree checks:
PED=$(ls *.ped | head -n1)
FAM=${PED//.ped/}
PED_MISMATCHED=${PED//.ped/-mismatched.ped}
cat $PED | head -n1 > $PED_MISMATCHED
cat $PED | tail -n+2 | grep father | sed $'s/\t1\t/\t2\t/g' | sed 's/father/grandmother/' >> $PED_MISMATCHED
cat $PED | tail -n+2 | grep father | tail -n+2 >> $PED_MISMATCHED
cat $PED | tail -n+2 | grep -v father >> $PED_MISMATCHED
gsutil cp $PED_MISMATCHED gs://cpg-seqr-test/datasets/$FAM/
```
