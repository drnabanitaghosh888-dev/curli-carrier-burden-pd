# PHASE 2H Nishiwaki Source Verification Runbook

## Scope

- Target cohort: `NISHIWAKI_MULTICOUNTRY_PD`
- Phase: `PHASE 2H-A`
- Objective: source verification and processed-table readiness planning only.

## Hard Constraints

- No downloads.
- No FASTQ/SRA/ENA/Qiita or raw-read processing.
- No Curli Carrier Burden analysis or replication statistics in this phase.
- No Phase 3 operations.

## Required Verification Items

1. Publication/source identity verification
2. Data repository provenance verification
3. Processed species abundance table verification
4. Sample metadata and PD/control label verification
5. Data dictionary / README verification
6. License / access terms verification
7. Optional run accession mapping check (future raw-read dependency)

## Required Processed Inputs for Phase 2H progression

- `sample_metadata.tsv`
- `species_abundance.tsv`
- `README_or_data_dictionary.txt`
- `license_or_access_terms.txt`

## Classification Logic

- `READY_FOR_MANUAL_DOWNLOAD`: source identity, license, species abundance, sample metadata, and PD/control labels verified.
- `PARTIAL_EVIDENCE_NEEDS_MANUAL_REVIEW`: some evidence present but critical items unresolved.
- `SOURCE_IDENTITY_UNCLEAR`: source provenance or dataset identity remains uncertain.
- `NOT_READY_SOURCE_NOT_VERIFIED`: core requirements (species metadata/license/labels) not yet verified.

## Current Starting Point (from existing Phase 2E metadata)

- Manual check status for Nishiwaki items is pending.
- Source blocking issue is `source_not_checked`.
- Working Phase 2H-A classification is therefore conservative: `NOT_READY_SOURCE_NOT_VERIFIED`.
