TEST_VERSION := v1-3
TEST_DATASET := NA12878-trio

ZORNITA_STARK_DATASET := zornitza-stark
ZORNITA_STARK_VERSION := v1-0

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
	--output-dir   "datasets/$(TEST_DATASET)/$(TEST_VERSION)" \
	--description  "seqr loader - $(TEST_DATASET) test - from GVCFs" \
	batch_seqr_loader/batch_workflow.py \
	--namespace    test \
	--version      $(TEST_VERSION) \
	--seqr-dataset $(TEST_DATASET) \
	--gvcf         'gs://cpg-seqr-test/batches/NA12878-trio/*.g.vcf.gz' \
	--ped-file     "gs://cpg-seqr-test/batches/NA12878-trio/NA12878-trio.ped" \
	--keep-scratch \
	--reuse

.PHONY: run_test_mismatched
run_test_mismatched:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "datasets/$(TEST_DATASET)/$(TEST_VERSION)" \
	--description  "seqr loader - $(TEST_DATASET) test - mismatched PED" \
	batch_seqr_loader/batch_workflow.py \
	--namespace    test \
	--version      $(TEST_VERSION) \
	--seqr-dataset $(TEST_DATASET)-mismatched \
	--gvcf         'gs://cpg-seqr-test/batches/NA12878-trio/*.g.vcf.gz' \
	--cram         "gs://cpg-seqr-test/batches/NA12878-trio/SS6004470.cram" \
	--ped-file     "gs://cpg-seqr-test/batches/NA12878-trio/NA12878-trio-mismatched.ped" \
	--keep-scratch \
	--reuse

.PHONY: run_test_extend
run_test_extend:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "datasets/$(TEST_DATASET)/$(TEST_VERSION)" \
	--description  "seqr loader - $(TEST_DATASET) test - extend" \
	batch_seqr_loader/batch_workflow.py \
	--namespace    test \
	--version      $(TEST_VERSION) \
	--seqr-dataset $(TEST_DATASET)-extend \
	--cram         "gs://cpg-seqr-test/batches/NA12878-trio/SS6004470.cram" \
	--ped-file     "gs://cpg-seqr-test/batches/NA12878-trio/NA12878-trio.ped" \
	--keep-scratch \
	--reuse

.PHONY: run_test_from_cram
run_test_from_cram:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "datasets/$(TEST_DATASET)/$(TEST_VERSION)" \
	--description  "seqr loader - $(TEST_DATASET) test - from CRAM" \
	batch_seqr_loader/batch_workflow.py \
	--namespace    test \
	--version      $(TEST_VERSION) \
	--seqr-dataset $(TEST_DATASET)-from-cram \
	--data-to-realign 'gs://cpg-seqr-test/batches/NA12878-trio/NA12878.cram' \
	--cram         'gs://cpg-seqr-test/batches/NA12878-trio/NA1289*.cram' \
	--cram         'gs://cpg-seqr-test/batches/NA12878-trio/SS6004470.cram' \
	--ped-file     "gs://cpg-seqr-test/batches/NA12878-trio/NA12878-trio.ped" \
	--keep-scratch \
	--reuse

.PHONY: run_zornitza-stark
run_zornitza-stark:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "datasets/$(ZORNITA_STARK_DATASET)/$(ZORNITA_STARK_VERSION)" \
	--description  "seqr loader - $(ZORNITA_STARK_DATASET)" \
	batch_seqr_loader/batch_workflow.py \
	--namespace    test \
	--version      $(ZORNITA_STARK_VERSION) \
	--seqr-dataset $(ZORNITA_STARK_DATASET) \
	--cram         'gs://cpg-seqr-upload-zornitza-stark/*.bam' \
	--ped-file     "gs://cpg-seqr-upload-zornitza-stark/cpg_acute-fixed.ped" \
	--reuse

.PHONY: run_zornitza-stark-kccg-gvcf
run_zornitza-stark-kccg-gvcf:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "datasets/$(ZORNITA_STARK_DATASET)/$(ZORNITA_STARK_VERSION)" \
	--description  "seqr loader - $(ZORNITA_STARK_DATASET) from KCCG GVCFs" \
	batch_seqr_loader/batch_workflow.py \
	--namespace    test \
	--version      $(ZORNITA_STARK_VERSION) \
	--seqr-dataset $(ZORNITA_STARK_DATASET)-kccg-gvcf \
	--gvcf         'gs://cpg-seqr-upload-zornitza-stark/*.g.vcf.gz' \
	--ped-file     "gs://cpg-seqr-upload-zornitza-stark/cpg_acute-fixed.ped" \
	--reuse
