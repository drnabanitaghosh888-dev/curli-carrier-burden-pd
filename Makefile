PYTHON := python

.PHONY: env-check phase1-fetch-dry-run phase1-fetch phase1-curate phase1-align phase1-hmms phase1-validate phase1-report phase1-all phase1-clean-temp

env-check:
	$(PYTHON) scripts/phase1_00_check_environment.py

phase1-fetch-dry-run:
	$(PYTHON) scripts/phase1_01_fetch_curli_sequences.py --dry-run

phase1-fetch:
	$(PYTHON) scripts/phase1_01_fetch_curli_sequences.py

phase1-curate:
	$(PYTHON) scripts/phase1_02_filter_and_curate.py

phase1-align:
	$(PYTHON) scripts/phase1_03_build_alignments.py

phase1-hmms:
	$(PYTHON) scripts/phase1_04_build_hmms.py

phase1-validate:
	$(PYTHON) scripts/phase1_05_validate_hmms.py

phase1-report:
	$(PYTHON) scripts/phase1_06_make_report.py

phase1-all: env-check phase1-fetch phase1-curate phase1-align phase1-hmms phase1-validate phase1-report
	@echo "Phase 1 complete."

phase1-clean-temp:
	find data/phase1 -type f -name "*.tmp*" -delete
