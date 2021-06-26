TEST_VERSION := v1-2

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
	--description "seqr loader - NA12878 trio test - from GVCFs" \
	batch_seqr_loader/batch_workflow.py \
	--namespace   test \
	--version     $(TEST_VERSION) \
	--dataset     "NA12878-trio" \
	--gvcf        'gs://cpg-seqr-test/batches/NA12878-trio/*.g.vcf.gz' \
	--ped-file    "gs://cpg-seqr-test/batches/NA12878-trio/NA12878-trio.ped" \
	--keep-scratch \
	--reuse

.PHONY: run_test_extend_with_cram
run_test_extend_with_cram:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir  "gs://cpg-seqr-test-tmp/hail" \
	--description "seqr loader - NA12878 trio test - extend with CRAM" \
	batch_seqr_loader/batch_workflow.py \
	--namespace   test \
	--version     $(TEST_VERSION) \
	--dataset     "NA12878-trio" \
	--cram        "gs://cpg-seqr-test/batches/NA12878-trio/SS6004470.cram" \
	--ped-file    "gs://cpg-seqr-test/batches/NA12878-trio/NA12878-trio.ped" \
	--keep-scratch \
	--reuse

.PHONY: run_test_mismatched
run_test_mismatched:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir  "gs://cpg-seqr-test-tmp/hail" \
	--description "seqr loader - NA12878 trio test - mismatched PED" \
	batch_seqr_loader/batch_workflow.py \
	--namespace   test \
	--version     $(TEST_VERSION) \
	--dataset     "NA12878-trio-mismatched" \
	--gvcf        'gs://cpg-seqr-test/batches/NA12878-trio/*.g.vcf.gz' \
	--cram        "gs://cpg-seqr-test/batches/NA12878-trio/SS6004470.cram" \
	--ped-file    "gs://cpg-seqr-test/batches/NA12878-trio/NA12878-trio-mismatched.ped" \
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
	--cram        'gs://cpg-seqr-upload-zornitza-stark/*.bam' \
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
	--gvcf        'gs://cpg-seqr-upload-zornitza-stark/*.g.vcf.gz' \
	--ped-file    "gs://cpg-seqr-upload-zornitza-stark/cpg_acute-fixed.ped" \
	--reuse
