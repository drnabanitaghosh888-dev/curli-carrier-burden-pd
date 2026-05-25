# PHASE2E Manual Source Verification Runbook

## Purpose
After Wallen discovery, replication is required to evaluate consistency of processed-table proxy signals in independent cohorts.

## Scope
- Manual source verification and triage only.
- Manual download planning only.
- No FASTQ/SRA download.
- No raw-read processing.

## Required Evidence Per Cohort
- Species abundance table availability
- Sample metadata availability
- PD/control label availability
- Data dictionary/README availability
- License/access terms clarity
- Run-accession mapping availability (for future Phase 3 only)

## Manual Source Inspection Steps
1. Review publication and supplementary materials.
2. Check repository records (Zenodo/Figshare/GitHub/OSF/other cited source).
3. Verify dataset identifiers and metadata provenance.
4. Record license/access terms explicitly.
5. Update:
   - `metadata/phase2e_replication_source_triage.tsv`
   - `metadata/phase2e_replication_manual_download_plan.tsv`

## Download Policy
- Only small processed tables may be downloaded.
- Download allowed only after source and license verification.
- Do not download FASTQ/SRA in Phase 2E-B.

## Readiness Decision
A cohort can be marked ready for processed-table replication only if:
- shotgun metagenomics
- species table available
- sample metadata + PD/control labels available
- access/reuse terms acceptable

Raw-read readiness remains pending until run-accession mapping is verified.
