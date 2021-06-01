VERSION := v0

.PHONY: package
package:
	rm -rf dist/*
	python setup.py sdist bdist_wheel
	twine upload dist/*

.PHONY: patch
patch:
	bump2version patch

.PHONY: minor
minor:
	bump2version minor

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
	--output-dir "gs://cpg-seqr-test/seqr_$(VERSION)/hail" \
	--description "test seqr loader" \
	batch_seqr_loader/batch_workflow.py \
	--gvcf-bucket gs://cpg-seqr-temporary/gvcf \
	--ped-file gs://cpg-seqr-test/gvcf/samples.ped \
	--dataset seqr \
	--work-bucket "gs://cpg-seqr-test/seqr_$(VERSION)/work" \
	--dest-mt-path "gs://cpg-seqr-test/seqr_$(VERSION)/output/annotated.mt" \
	--genomicsdb-bucket gs://cpg-seqr-test/seqr_$(VERSION)/genomicsdb \
	--keep-scratch
