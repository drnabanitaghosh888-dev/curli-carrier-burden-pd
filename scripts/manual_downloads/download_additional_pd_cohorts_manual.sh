#!/usr/bin/env bash
set -euo pipefail

if [[ "${RUN_MANUAL_DOWNLOAD:-}" != "YES" ]]; then
  echo "[WARNING] RUN_MANUAL_DOWNLOAD is not YES. Exiting without download."
  exit 0
fi

echo "[INFO] Manual template only for additional independent PD cohorts."
echo "[INFO] Candidate cohorts require manual accession/source verification first."

# ------------------------------------------------------------------
# TEMPLATE (commented intentionally)
# 1) Confirm accession/source from publications or repositories.
# 2) Confirm sequencing type is shotgun metagenomics for raw csg screening.
# 3) Exclude 16S-only cohorts from raw csg detection workflows.
# 4) Download only selected sample subsets for compute-safe targeted analysis.
#
# Example placeholders (DO NOT EXECUTE AS-IS):
# # prefetch <CANDIDATE_SRR_ID>
# # fasterq-dump <CANDIDATE_SRR_ID> --split-files --threads 4
# # gzip <CANDIDATE_SRR_ID>_1.fastq <CANDIDATE_SRR_ID>_2.fastq
# ------------------------------------------------------------------

echo "[DONE] Manual candidate plan acknowledged. No download performed."
