# PHASE2E Cohort-by-Cohort Source Verification Runbook

## Objective
Verify replication cohort source evidence one cohort at a time before any manual processed-table download.

## Safety Rules
- Do not download FASTQ/SRA.
- Do not run raw-read processing.
- Do not run metagenomic profilers or sequence search tools.
- Only small processed tables are eligible for manual download after source and license verification.

## Cohort-First Workflow
For each cohort:
1. Inspect publication page.
2. Inspect supplementary tables/files.
3. Inspect cited data repository pages.
4. Record exact URL and file/resource name.
5. Verify license/access/reuse terms.
6. Check whether sample IDs can link abundance tables to metadata.
7. Update verification form rows.

## Required Ready Criteria for Processed-Table Replication
A cohort can be marked ready only if all are verified:
- species abundance table
- sample metadata
- PD/control labels
- acceptable access/reuse terms

## If Only Raw SRA Exists
- Mark not ready for processed-table replication.
- Do not download raw files in this phase.
- Keep ready_for_raw_read_validation as pending until mapping and access are verified.
