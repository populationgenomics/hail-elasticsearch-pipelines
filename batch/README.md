# Batch loading pipeline

Loading pipeline, but formatted for Hail Batch, ultimately for automated ingestion of data into SEQR.

For the purposes of this, we'll say you can only run it through the analysis-runner.

Example:

```
analysis-runner \
    --dataset <dataset> \
    batch/load.py --
```

This script will automatically start a Dataproc cluster where relevant.

## Customisation

This pipeline should be fairly configurable in a few ways:

- Start from the sparse matrix instead of a well-formed VCF
- Upload to ElasticSearch if elastic search credentials have been provided.

Later additions:

- Load pedigree data if provided ([using this process](https://centrepopgen.slack.com/archives/C01R7CKJGHM/p1618551394039300))