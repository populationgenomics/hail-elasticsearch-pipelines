VERSION := v0

.PHONY: commit
commit:
	git all
	git commit -m "$(msg)" --no-verify
	git push

.PHONY: run_test
run_test:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir "gs://cpg-seqr-temporary/batch-loader-output/hail" \
	--description "test seqr loader" \
	batch_seqr_loader/batch_workflow.py \
	--gvcf-bucket gs://cpg-seqr-temporary/gvcf \
	--ped-file gs://cpg-seqr-temporary/gvcf/pedigree.ped \
	--dataset seqr \
	--work-bucket "gs://cpg-seqr-temporary/batch-loader-output/$(VERSION)/work" \
	--dest-mt-path "gs://cpg-seqr-temporary/batch-loader-output/$(VERSION)/annotated.mt" \
        --keep-scratch
