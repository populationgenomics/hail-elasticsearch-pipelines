TEST_VERSION := v1-7
TEST_DATASET := na12878-trio

ACUTE_CARE_DATASET := acute-care
ACUTE_CARE_VERSION := v1-0

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
	--cram         'gs://cpg-seqr-test/batches/NA12878-trio/NA1289*.cram' \
	--cram         'gs://cpg-seqr-test/batches/NA12878-trio/SS6004470.cram' \
	--ped-file     "gs://cpg-seqr-test/batches/NA12878-trio/NA12878-trio.ped" \
	--keep-scratch \
	--reuse

.PHONY: run_test_from_cram_with_realignment
run_test_from_cram_with_realignment:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "datasets/$(TEST_DATASET)/$(TEST_VERSION)" \
	--description  "seqr loader - $(TEST_DATASET) test - from CRAM with realignment" \
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

.PHONY: run_test_from_cram_tiny
run_test_from_cram_tiny:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "datasets/$(TEST_DATASET)/$(TEST_VERSION)" \
	--description  "seqr loader - $(TEST_DATASET) test - from CRAM - tiny" \
	batch_seqr_loader/batch_workflow.py \
	--namespace    test \
	--version      $(TEST_VERSION) \
	--seqr-dataset $(TEST_DATASET)-from-cram-tiny \
	--data-to-realign gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878.cram \
	--ped-file     "gs://cpg-seqr-test/batches/NA12878-trio/NA12878-trio.ped" \
	--keep-scratch \
	--reuse

.PHONY: run_test_from_fq_tiny
run_test_from_fq_tiny:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "datasets/$(TEST_DATASET)/$(TEST_VERSION)" \
	--description  "seqr loader - $(TEST_DATASET) test - from FQ - tiny" \
	batch_seqr_loader/batch_workflow.py \
	--namespace    test \
	--version      $(TEST_VERSION) \
	--seqr-dataset $(TEST_DATASET)-from-fq-tiny \
	--data-to-realign 'gs://cpg-seqr-test/batches/NA12878-trio-tiny/NA12878*.fq' \
	--ped-file     "gs://cpg-seqr-test/batches/NA12878-trio/NA12878-trio.ped" \
	--keep-scratch \
	--reuse

.PHONY: run_acute_care
run_acute_care:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "datasets/$(ACUTE_CARE_DATASET)/$(ACUTE_CARE_VERSION)" \
	--description  "seqr loader - $(ACUTE_CARE_DATASET) from KCCG GVCFs and Simons fastqs (1 sample)" \
	batch_seqr_loader/batch_workflow.py \
	--namespace    test \
	--version      $(ACUTE_CARE_VERSION) \
	--seqr-dataset $(ACUTE_CARE_DATASET) \
	--gvcf         'gs://cpg-seqr-upload-zornitza-stark/*.g.vcf.gz' \
	--data-to-realign 'gs://cpg-seqr-upload-zornitza-stark/cpg_acute_20210727_185421/200721_A00692_0122_ML206418_20W001106-FAM000553_MAN-20200721_NEXTERAFLEXWGS_*.fastq.gz' \
	--ped-file     "gs://cpg-seqr-upload-zornitza-stark/cpg_acute-fixed.ped" \
	--reuse

.PHONY: run_test_sm_workflow
run_test_sm_workflow:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "datasets/run_test_sm_workflow" \
	--description  "test SM workflow" \
	batch_seqr_loader/test/test_communicate_with_sm.py
