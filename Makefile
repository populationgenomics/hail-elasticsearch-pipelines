VERSION := v1_2

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
	--output-dir "gs://cpg-seqr-test-tmp/seqr_$(VERSION)/hail" \
	--description "test seqr loader - batch small1" \
	batch_seqr_loader/batch_workflow.py \
	--gvcf-bucket "gs://cpg-seqr-test/gvcf/small1" \
	--cram-bucket "gs://cpg-seqr-test/cram/small1" \
	--ped-file    "gs://cpg-seqr-test/gvcf/small1/samples.ped" \
	--data-bucket "gs://cpg-seqr-test/data/test-$(VERSION)" \
	--work-bucket "gs://cpg-seqr-test-tmp/work/loader-$(VERSION)/work-small1-cram" \
	--keep-scratch \
	--reuse

.PHONY: run_test_extend
run_test_extend:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir "gs://cpg-seqr-test-tmp/seqr_$(VERSION)/hail" \
	--description "test seqr loader - extend with batch small2" \
	batch_seqr_loader/batch_workflow.py \
	--gvcf-bucket  "gs://cpg-seqr-test/gvcf/small2" \
	--ped-file     "gs://cpg-seqr-test/gvcf/small2/samples.ped" \
	--data-bucket  "gs://cpg-seqr-test/data/test-$(VERSION)" \
	--work-bucket  "gs://cpg-seqr-test-tmp/work/loader-$(VERSION)/work-small2" \
	--keep-scratch \
	--reuse

.PHONY: run_full_family
run_test:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir "gs://cpg-seqr-test-tmp/seqr_$(VERSION)/hail" \
	--description "test seqr loader - families" \
	batch_seqr_loader/batch_workflow.py \
	--gvcf-bucket "gs://cpg-seqr-test/gvcf/families" \
	--ped-file    "gs://cpg-seqr-test/gvcf/families/samples.ped" \
	--data-bucket "gs://cpg-seqr-test/data/families-$(VERSION)" \
	--work-bucket "gs://cpg-seqr-test-tmp/work/loader-$(VERSION)/work-families" \
	--keep-scratch \
	--reuse
