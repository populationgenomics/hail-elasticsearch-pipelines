# Batch loading pipeline

A loading pipeline, but formatted for Hail Batch, ultimately for automated ingestion of data into Seqr.

The inputs are buckets with GVCFs, CRAMs or BAMs, as well as an optional PED file. It calls the variants if needed, and iteratively adds them into a [GenomicsDB](https://github.com/Intel-HLS/GenomicsDB/wiki), located on a bucket `gs://cpg-seqr-main/datasets`. Then it jointly re-genotypes all variants in a DB together, converts them into a Hail MatrixTable, and annotates.

Analysis runner must be used to run the workflow. For examples of commands, see the [Makefile](../Makefile) at the root of the repository. 


## Further work

* Upload to Elasticsearch if corresponding credentials have been provided.
* Load pedigree data if provided ([using this process](https://centrepopgen.slack.com/archives/C01R7CKJGHM/p1618551394039300)).
* Support re-aligning input CRAMs/BAMs/FASTQs if needed.
