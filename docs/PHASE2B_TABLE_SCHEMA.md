# PHASE 2B Table Schema Guidance

## Required Table Formats

- UTF-8 TSV preferred (CSV accepted if documented).
- One header row, no merged cells, explicit sample IDs.

## Sample ID Harmonization Rules

- `sample_id` must be stable across abundance and metadata tables.
- If alternate IDs exist, provide mapping table fields in sample metadata.
- Preserve cohort-specific IDs but add canonical harmonized ID when possible.

## Abundance Table Expectations

- Species/genus abundance tables: `sample_id`, `taxon`, `abundance`.
- Gene family/pathway tables: `sample_id`, `gene_family/pathway`, `abundance`.
- Include normalization metadata when available.

## Metadata Expectations

- Minimum: PD/control status, sequencing type, sample type, cohort identifier.
- Preferred: age, sex, geography, batch, medication, disease duration.

## Clinical Phenotype Expectations

Preferred vagal/body-first proxy fields:
- constipation
- REM_sleep_behavior_disorder
- hyposmia
- autonomic_symptoms
- disease_duration
- UPDRS
- medication
- age
- sex
- geography
- batch

## Handling Missing Vagal/Body-First Metadata

- Do not exclude the cohort automatically.
- Mark unavailable fields explicitly.
- Use available fields for partial analyses and sensitivity reporting.

## Core Eligibility Rule

Only shotgun metagenomic datasets are eligible for the core curli/csg analysis.
16S datasets lack gene-level resolution and are excluded from curli gene
detection, csg operon analysis, raw csg screening, strain-level reconstruction,
and primary association testing. They may be retained only as optional
background ecological context and should not be used in the main analytical
workflow.
