# PHASE 2A Dataset Acquisition Plan

## Scope

Phase 2A defines a metadata-first acquisition plan for PD gut metagenomic datasets.
No raw read download or heavy compute is executed in this phase.

## Why Wallen PRJNA834801 Is First Priority

- Large PD/control stool shotgun metagenomics cohort.
- Strong candidate for both processed-table analysis and selected raw-read follow-up.
- Practical first anchor dataset for curli-oriented downstream screening design.
- Wallen et al. enrolled 490 PD and 234 controls (724 total), but exact downloadable
  run/sample counts must be verified from accession metadata before any raw-read step.

## Why Processed Tables First

- Processed tables provide rapid hypothesis triage at low compute cost.
- They enable early filtering of candidate taxa/signatures before expensive raw-read workflows.
- This reduces unnecessary download/storage burden and limits premature heavy compute.

## Why Only a Selected Subset Should Be Downloaded Later

- Raw shotgun files are large and expensive to process.
- Phase 2B should preselect biologically informative samples/cohorts.
- Targeted subset selection keeps compute, storage, and runtime within MacBook-safe limits.

## How Phase 2B Will Select Samples For Targeted csg Screening

- Prioritize shotgun PD/control cohorts with clear accession provenance.
- Use clinical metadata availability and phenotype relevance when present.
- Prefer cohorts with compatible metadata fields for body-first/vagal analyses.
- Defer uncertain cohorts as `candidate_pending_manual_verification`.

## Compute-Safety Notes (MacBook)

- Keep Phase 2A metadata-only and validation-only.
- Do not run SRA retrieval, metagenome profiling, or HMM scanning in Phase 2A.
- Keep thread counts conservative (default 4) when heavy steps start in later phases.
- Start later raw-read processing with a small pilot subset before scaling.
