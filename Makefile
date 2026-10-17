PYTHON ?= python
SNAKEMAKE ?= snakemake
CORES ?= 4
USE_CONDA ?= --use-conda

.PHONY: check-env smoke-test mock-test test dry-run run

check-env:
	VACDES_USE_CONDA=$(if $(strip $(USE_CONDA)),1,0) $(PYTHON) scripts/check_environment.py

# The smoke test builds an isolated temporary config and checks the complete DAG.
# No workflow jobs or external predictors are executed.
smoke-test:
	$(PYTHON) scripts/smoke_test.py

mock-test:
	VACDES_USE_CONDA=0 VACDES_CHECK_SCOPE=mock-test $(PYTHON) scripts/check_environment.py
	$(PYTHON) scripts/mock_predictor_test.py

test:
	$(PYTHON) -m pytest -q

dry-run: check-env
	$(SNAKEMAKE) $(USE_CONDA) --cores $(CORES) --dry-run

run: check-env
	$(SNAKEMAKE) $(USE_CONDA) --cores $(CORES)
