VERSION := v1-0

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
	--output-dir  "gs://cpg-seqr-test/data/test-$(VERSION)" \
	--description "test seqr loader - test" \
	batch_seqr_loader/batch_workflow.py \
	--namespace test \
	--gvcf-bucket "gs://cpg-seqr-test/datasets/BB01" \
	--ped-file    "gs://cpg-seqr-test/datasets/BB01/BB01.ped" \
	--data-bucket "gs://cpg-seqr-test/data/test-$(VERSION)" \
	--work-bucket "gs://cpg-seqr-test-tmp/work/seqr-loader-$(VERSION)" \
	--keep-scratch \
	--reuse

.PHONY: run_test_extend
run_test_extend:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir  "gs://cpg-seqr-test/data/test-$(VERSION)" \
	--description "test seqr loader - extend" \
	batch_seqr_loader/batch_workflow.py \
	--namespace test \
	--gvcf-bucket "gs://cpg-seqr-test/datasets/BB02" \
	--ped-file    "gs://cpg-seqr-test/datasets/BB02/BB02.ped" \
	--data-bucket "gs://cpg-seqr-test/data/test-$(VERSION)" \
	--work-bucket "gs://cpg-seqr-test-tmp/work/seqr-loader-$(VERSION)-extend" \
	--keep-scratch \
	--reuse

.PHONY: run_test_mismatched
run_test_mismatched:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir  "gs://cpg-seqr-test/data/test-$(VERSION)" \
	--description "test seqr loader - mismatched" \
	batch_seqr_loader/batch_workflow.py \
	--namespace test \
	--gvcf-bucket "gs://cpg-seqr-test/datasets/BB01" \
	--ped-file    "gs://cpg-seqr-test/datasets/BB01/BB01-mismatched.ped" \
	--data-bucket "gs://cpg-seqr-test/data/test-$(VERSION)" \
	--work-bucket "gs://cpg-seqr-test-tmp/work/seqr-loader-$(VERSION)-mismatched" \
	--keep-scratch \
	--reuse

.PHONY: run_zornitza-stark
run_zornitza-stark:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir  "gs://cpg-seqr-test/data/zornitza-stark-v0" \
	--description "test seqr loader - zornitza-stark" \
	batch_seqr_loader/batch_workflow.py \
	--namespace test \
	--bam-bucket  "gs://cpg-seqr-upload-zornitza-stark" \
	--ped-file    "gs://cpg-seqr-upload-zornitza-stark/cpg_acute.ped" \
	--data-bucket "gs://cpg-seqr-test/data/zornitza-stark-v0" \
	--work-bucket "gs://cpg-seqr-test-tmp/work/seqr-loader-$(VERSION)/work-zornitza-stark" \
	--reuse
