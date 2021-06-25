TEST_VERSION := v1-0

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
	--output-dir  "gs://cpg-seqr-test-tmp/hail" \
	--description "test seqr loader - test" \
	batch_seqr_loader/batch_workflow.py \
	--namespace   test \
	--version     $(TEST_VERSION) \
	--dataset     "BB01-BB02" \
	--gvcf-bucket "gs://cpg-seqr-test/batches/BB01" \
	--ped-file    "gs://cpg-seqr-test/batches/BB01/BB01.ped" \
	--keep-scratch \
	--reuse

.PHONY: run_test_extend
run_test_extend:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir  "gs://cpg-seqr-test-tmp/hail" \
	--description "test seqr loader - extend" \
	batch_seqr_loader/batch_workflow.py \
	--namespace   test \
	--version     $(TEST_VERSION) \
	--dataset     "BB01-BB02" \
	--gvcf-bucket "gs://cpg-seqr-test/batches/BB02" \
	--ped-file    "gs://cpg-seqr-test/batches/BB02/BB02.ped" \
	--keep-scratch \
	--reuse

.PHONY: run_test_mismatched
run_test_mismatched:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir  "gs://cpg-seqr-test-tmp/hail" \
	--description "test seqr loader - mismatched" \
	batch_seqr_loader/batch_workflow.py \
	--namespace   test \
	--version     $(TEST_VERSION) \
	--dataset     "BB01-BB02-mismatched" \
	--gvcf-bucket "gs://cpg-seqr-test/batches/BB01" \
	--ped-file    "gs://cpg-seqr-test/batches/BB01/BB01-mismatched.ped" \
	--keep-scratch \
	--reuse

.PHONY: run_test_cram
run_test_cram:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir  "gs://cpg-seqr-test-tmp/hail" \
	--description "test seqr loader - cram" \
	batch_seqr_loader/batch_workflow.py \
	--namespace   test \
	--version     $(TEST_VERSION) \
	--dataset     "NA12878-cram" \
	--cram-bucket "gs://cpg-seqr-test/batches/NA12878-cram" \
	--keep-scratch \
	--reuse

.PHONY: run_zornitza-stark
run_zornitza-stark:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir  "gs://cpg-seqr-test-tmp/hail" \
	--description "test seqr loader - zornitza-stark" \
	batch_seqr_loader/batch_workflow.py \
	--namespace   test \
	--version     "v1-0" \
	--dataset     "zornitza-stark" \
	--bam-bucket  "gs://cpg-seqr-upload-zornitza-stark" \
	--ped-file    "gs://cpg-seqr-upload-zornitza-stark/cpg_acute-fixed.ped" \
	--reuse

.PHONY: run_zornitza-stark-kccg-gvcf
run_zornitza-stark-kccg-gvcf:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir  "gs://cpg-seqr-test-tmp/hail" \
	--description "seqr loader - zornitza-stark KCCG GVCFs" \
	batch_seqr_loader/batch_workflow.py \
	--namespace   test \
	--version     "v1-0" \
	--dataset     "zornitza-stark-kccg-gvcf" \
	--gvcf-bucket "gs://cpg-seqr-upload-zornitza-stark" \
	--ped-file    "gs://cpg-seqr-upload-zornitza-stark/cpg_acute-fixed.ped" \
	--reuse
