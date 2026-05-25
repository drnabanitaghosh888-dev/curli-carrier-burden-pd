#!/usr/bin/env bash
set -euo pipefail

if [[ "${RUN_MANUAL_DOWNLOAD:-}" != "YES" ]]; then
  echo "[WARNING] RUN_MANUAL_DOWNLOAD is not YES. Exiting without download."
  exit 0
fi

echo "[INFO] Manual template only. No download commands are executed by default."
echo "[INFO] Dataset: Wallen et al. PD shotgun metagenomics cohort (PRJNA834801)"

# ------------------------------------------------------------------
# TEMPLATE (commented intentionally)
# 1) Verify access/consent constraints for PRJNA834801.
# 2) Use local institutional workflow for controlled or open retrieval.
# 3) Record retrieved run accessions into metadata/provenance logs.
# 4) Download only selected subset approved for Phase 2B targeted screening.
#
# Example placeholders (DO NOT EXECUTE AS-IS):
# # prefetch <SRR_ID>
# # fasterq-dump <SRR_ID> --split-files --threads 4
# # gzip <SRR_ID>_1.fastq <SRR_ID>_2.fastq
# ------------------------------------------------------------------

echo "[DONE] Manual plan acknowledged. No download performed."
