VERSION := v0

.PHONY: package
package:
	rm -rf dist/*
	python setup.py sdist bdist_wheel
	twine upload dist/*

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
	--output-dir "gs://cpg-seqr-temporary/seqr/hail" \
	--description "test seqr loader" \
	batch_seqr_loader/batch_workflow.py \
	--gvcf-bucket gs://cpg-seqr-temporary/gvcf \
	--ped-file gs://cpg-seqr-temporary/gvcf/pedigree.ped \
	--dataset seqr \
	--work-bucket "gs://cpg-seqr-temporary/seqr2/work" \
	--dest-mt-path "gs://cpg-seqr-temporary/seqr2/output/annotated.mt" \
	--genomicsdb-bucket gs://cpg-seqr-temporary/seqr2/genomicsdb \
	--keep-scratch \
	--reuse-scratch-run-id a628bc
