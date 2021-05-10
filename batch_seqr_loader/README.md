# Batch loading pipeline

Loading pipeline, but formatted for Hail Batch, ultimately for automated ingestion of data into SEQR.

For the purposes of this, we'll say you can only run it through the analysis-runner.

Example:

```
analysis-runner \
    --dataset seqr \
    --output-dir "gs://cpg-seqr-temporary/loader-test" \
    --description "test seqr loader" \
    --access-level test batch_seqr_loader/batch_workflow.py \
batch_workflow.py \
    --gvcf-bucket gs://playground-au/seqr/gvcf/ \
    --ped-file gs://playground-au/seqr/samples.ped \
    --dest-mt-path gs://cpg-seqr-temporary/loader-test/NA12878_trio.mt \
    --work-bucket gs://cpg-seqr-temporary/loader-test/work \
    --keep-scratch
```

This script will automatically start a Dataproc cluster where relevant.

## Customisation

This pipeline should be fairly configurable in a few ways:

- Start from the sparse matrix instead of a well-formed VCF
- Upload to ElasticSearch if elastic search credentials have been provided.

Later additions:

- Load pedigree data if provided ([using this process](https://centrepopgen.slack.com/archives/C01R7CKJGHM/p1618551394039300))
