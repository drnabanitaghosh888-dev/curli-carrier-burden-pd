# PHASE2C Cohort Decision Rules

## Ready for Phase 2D Processed-Table Analysis

A cohort is ready only if all are true:
- `sequencing_type = shotgun_metagenomics`
- `processed_species_table_status = available`
- `sample_metadata_status = available`
- `pd_control_label_status = available`
- `access_restriction` is acceptable for analysis

## Ready for Phase 3 Raw-Read Subset Selection

A cohort is ready only if all are true:
- `run_accession_mapping_status = available`
- `sample_metadata_status = available`
- `pd_control_label_status = available`
- access/reuse terms permit download/use
- Phase 2D processed-table analysis suggests biologically informative strata

## Not Ready Conditions

A cohort is not ready if any applies:
- processed tables are missing
- metadata cannot be linked to samples
- access is restricted and not resolvable
- sequencing type is 16S-only
