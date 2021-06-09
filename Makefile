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
	--output-dir "gs://cpg-seqr-test-tmp/seqr_$(VERSION)/hail" \
	--description "test seqr loader - batch small1" \
	batch_seqr_loader/batch_workflow.py \
	--gvcf-bucket gs://cpg-seqr-test/gvcf/small1 \
	--ped-file gs://cpg-seqr-test/gvcf/small1/samples.ped \
	--dataset seqr \
	--work-bucket "gs://cpg-seqr-test/seqr_$(VERSION)/work" \
	--dest-mt-path "gs://cpg-seqr-test/seqr_$(VERSION)/output/annotated.mt" \
	--genomicsdb-bucket gs://cpg-seqr-test/seqr_$(VERSION)/genomicsdb \
	--keep-scratch \
	--disable-validation --skip-vep

.PHONY: run_test_extend
run_test_extend:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir "gs://cpg-seqr-test/seqr_$(VERSION)/hail" \
	--description "test seqr loader - extend with batch small2" \
	batch_seqr_loader/batch_workflow.py \
	--gvcf-bucket gs://cpg-seqr-test/gvcf/small2 \
	--ped-file gs://cpg-seqr-test/gvcf/small2/samples.ped \
	--dataset seqr \
	--work-bucket "gs://cpg-seqr-test/seqr_$(VERSION)/work-withsmall2" \
	--dest-mt-path "gs://cpg-seqr-test/seqr_$(VERSION)/output/annotated-withsmall2.mt" \
	--genomicsdb-bucket gs://cpg-seqr-test/seqr_$(VERSION)/genomicsdb \
	--keep-scratch \
	--disable-validation --skip-vep
