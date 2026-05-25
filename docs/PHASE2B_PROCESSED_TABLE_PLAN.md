# PHASE 2B Processed-Table Plan

## Scope

Phase 2B prepares processed-table acquisition and validation scaffolding only.
No FASTQ/SRA download and no raw metagenomic compute are executed in this phase.

## Why Processed Tables Before FASTQ

- Lower compute and storage burden.
- Faster multi-cohort harmonization and hypothesis triage.
- Better prioritization for later targeted raw-read csg follow-up.

## Why Wallen Is Discovery Anchor

- Largest immediately actionable anchor cohort in this project context.
- Supports discovery-stage processed-table analyses while accessions and metadata are audited.

## Why Replication Cohorts Are Required

- Final study cannot rely on Wallen-only findings.
- Replication cohorts provide geographic and cohort-design robustness.

## Candidate Cohort Verification Status

Replication cohorts remain `candidate_pending_manual_verification` until:
- accession/source is confirmed,
- processed tables are confirmed available,
- metadata suitability is confirmed.

## 16S-only Cohorts

16S-only cohorts are excluded from raw csg detection.
They may be used as broad ecological/taxonomic context only.

Scientific rule: Only shotgun metagenomic datasets are eligible for the core
curli/csg analysis. 16S datasets lack gene-level resolution and are excluded
from curli gene detection, csg operon analysis, raw csg screening, strain-level
reconstruction, and primary association testing. They may be retained only as
optional background ecological context and should not be used in the main
analytical workflow.

## Compute-Safety Rules

- No heavy tools in Phase 2B (HMMER, Kraken, HUMAnN, MetaPhlAn, Bowtie2, etc.).
- No SRA retrieval in Phase 2B.
- Keep this phase metadata/schema/validation focused.
