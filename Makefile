TEST_VERSION := v2-3
PROD_VERSION := v3-0

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

.PHONY: run_test_sm_workflow
run_test_sm_workflow:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "datasets/run_test_sm_workflow" \
	--description  "test SM workflow" \
	batch_seqr_loader/test/batch_test_simulate_sm_workfow.py

.PHONY: run_seqr_loader_test
run_seqr_loader_test:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "seqr_loader_test" \
	--description  "seqr loader - test $(TEST_VERSION)" \
	batch_seqr_loader/batch_workflow.py \
	--namespace test \
	--analysis-project seqr-test \
	--input-project acute-care \
	--input-project perth-neuro \
	--output-version $(TEST_VERSION) \
	--reuse \
	--keep-scratch

.PHONY: run_seqr_loader_prod
run_seqr_loader_prod:
	analysis-runner \
	--dataset seqr \
	--access-level standard \
	--output-dir   "seqr_loader" \
	--description  "seqr loader $(PROD_VERSION)" \
	batch_seqr_loader/batch_workflow.py \
	--namespace main \
	--analysis-project seqr \
	--input-project acute-care \
	--input-project perth-neuro \
	--output-version $(PROD_VERSION) \
	--reuse \
	--keep-scratch \
	--skip-check-inputs-existence \
	--hc-shards-num 10

.PHONY: run_seqr_loader_prod_test
run_seqr_loader_prod_test:
	analysis-runner \
	--dataset seqr \
	--access-level test \
	--output-dir   "seqr_loader" \
	--description  "seqr loader $(PROD_VERSION)" \
	batch_seqr_loader/batch_workflow.py \
	--namespace test \
	--analysis-project seqr \
	--input-project acute-care \
	--output-version $(PROD_VERSION) \
	-S CPG11783 \
    --start-from-stage joint_calling \
	--reuse \
	--keep-scratch
