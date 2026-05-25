# Phase 2D Wallen Discovery Analysis Plan

## Scope
Phase 2D is a processed-table-only discovery framework for the Wallen discovery anchor cohort.

## Inputs
- `sample_metadata.tsv`
- `clinical_metadata.tsv`
- `species_abundance.tsv`
- `genus_abundance.tsv`
- `gene_family_abundance.tsv`
- `pathway_abundance.tsv`
- `metadata/curli_candidate_taxa_from_phase1.tsv`

## Outputs
- Input validation reports
- Clean analysis metadata table
- Curli candidate species-matching table
- Sample-level Curli Carrier Burden table
- Discovery statistics summary table
- Final Phase 2D summary markdown and manifest

## Key Rules
- Processed-table analysis only.
- No FASTQ/SRA or raw-read operations in Phase 2D.
- No csg-gene presence claims from Phase 2D.
- No operon-integrity claims from Phase 2D.
- No Phase 3 claims while run accession mapping is unavailable.

## Wallen Role
Wallen is the discovery anchor for processed-table discovery. Replication cohorts remain separate follow-up work.
