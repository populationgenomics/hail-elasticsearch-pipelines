# Batch loading pipeline

A loading pipeline, but formatted for Hail Batch, ultimately for automated ingestion of data into Seqr.

The input is a set of GVCFs. It iteratively adds them into a [GenomicsDB](https://github.com/Intel-HLS/GenomicsDB/wiki) located on a bucket specified with `--genomicsdb-bucket`, jointly re-genotypes, converts into a Hail MatrixTable and annotates.

To run, use the analysis runner:

```sh
VERSION=v1.01
analysis-runner \
    --dataset seqr \
    --access-level test \
    --output-dir "gs://cpg-seqr-test-tmp/seqr_${VERSION}/hail" \
    --description "test seqr loader - batch small1" \
batch_seqr_loader/batch_workflow.py \
    --gvcf-bucket gs://cpg-seqr-test/gvcf/small1 \
    --ped-file gs://cpg-seqr-test/gvcf/small1/samples.ped \
    --dataset seqr \
    --work-bucket "gs://cpg-seqr-test/seqr_${VERSION}/work" \
    --dest-mt-path "gs://cpg-seqr-test/seqr_${VERSION}/output/annotated.mt" \
    --genomicsdb-bucket gs://cpg-seqr-test/seqr_${VERSION}/genomicsdb
```

This script will automatically start a Dataproc cluster where relevant.

## Customisation

This pipeline should be fairly configurable in a few ways:

- Upload to ElasticSearch if elastic search credentials have been provided.

Later additions:

- Load pedigree data if provided ([using this process](https://centrepopgen.slack.com/archives/C01R7CKJGHM/p1618551394039300))
