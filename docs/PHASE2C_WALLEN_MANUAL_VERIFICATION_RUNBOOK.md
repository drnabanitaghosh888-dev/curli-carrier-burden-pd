# PHASE 2C-W1 Wallen Manual Verification Runbook

## Scope

This runbook is for manual verification of processed-table and metadata provenance
for the Wallen discovery anchor cohort (`PRJNA834801`) only.

## Required Manual Sources to Inspect

1. Wallen Nature Communications article
2. Zenodo study record (if available for this cohort)
3. NCBI BioProject `PRJNA834801`
4. Supplementary tables and/or repository files linked from the publication

## Core Safety Rules

- Processed tables may be downloaded manually **only if small and clearly licensed**.
- FASTQ/SRA files must **not** be downloaded in this phase.
- No raw-read processing should be performed in this phase.

## Verification Workflow

1. Open `metadata/phase2c_wallen_manual_source_verification.tsv`.
2. For each item, manually record:
   - `verified_source_url`
   - `verified_file_name_at_source`
   - update `verification_status`
3. Update `metadata/phase2c_wallen_file_renaming_plan.tsv` with source file names.
4. Place manually approved processed tables into:
   - `data/phase2b/processed_tables/wallen_prjna834801/`
5. Run local existence checker:
   - `python scripts/phase2c_01_check_wallen_local_files.py`

## Readiness Rules

Wallen cannot be marked ready for **Phase 2D** until:
- `species_abundance.tsv` is verified and present
- `sample_metadata.tsv` is verified and present
- PD/control labels are verified in metadata

Wallen cannot be marked ready for **Phase 3** until:
- `run_accession_mapping.tsv` is verified and present

## Notes

This is provenance and table-availability verification only.  
Do not perform metagenomic compute in this step.
