# PHASE2C Manual Processed-Table Verification

## Scope

Phase 2C is metadata/provenance verification only.  
No FASTQ/SRA download and no raw-read processing are performed.

## Why Verification Before Analysis

- Prevents invalid comparisons due to missing/ambiguous metadata.
- Ensures PD/control labels can be linked to sample IDs.
- Confirms table provenance and access terms before downstream analysis.

## Cohort Positioning

- Wallen (`WALLEN_PRJNA834801`) remains discovery anchor.
- Replication cohorts remain candidates until source tables and metadata are verified.

## 16S Exclusion Rule

16S-only datasets are excluded from active Phase 2C verification and core curli/csg workflow
because they lack gene-level resolution for curli gene detection and operon-level analysis.

## Evidence Required Before Phase 2D

- Species-abundance table available
- Sample metadata available
- PD/control labels available
- Acceptable access/license terms
- Cohort-level provenance documentation
