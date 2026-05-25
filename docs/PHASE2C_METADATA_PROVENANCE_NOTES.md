# PHASE2C Metadata Provenance Notes

## Sample ID Harmonization

- Define one canonical `sample_id` per cohort.
- Preserve original IDs in mapping columns.
- Require explicit mapping where run accessions differ from sample IDs.

## Accession Mapping Rules

- `run_accession_mapping` must resolve sample-to-run relationships unambiguously.
- Ambiguous or missing mappings block Phase 3 readiness.

## Provenance Requirements

- Keep source repository/publication reference for each processed table.
- Keep access/license terms with timestamped verification notes.

## Clinical Variable Uncertainty

- Missing clinical fields should be explicitly marked as unavailable.
- Do not infer values from secondary narrative text.

## Vagal/Body-First Phenotype Uncertainty

- Treat constipation/RBD/hyposmia/autonomic fields as optional unless verified.
- Mark uncertain fields as pending and exclude from confirmatory claims.

## No-Overclaiming Rule

No cohort should be considered analysis-ready without explicit table+metadata evidence.
No missing provenance should be silently assumed.
