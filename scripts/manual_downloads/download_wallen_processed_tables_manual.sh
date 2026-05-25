#!/usr/bin/env bash
set -euo pipefail

if [ "${RUN_MANUAL_DOWNLOAD:-NO}" != "YES" ]; then
  echo "Manual download guard active. Set RUN_MANUAL_DOWNLOAD=YES to run."
  exit 1
fi

echo "[INFO] Processed-table manual download template only."
echo "[INFO] Dataset: WALLEN_PRJNA834801"

# Placeholder commands only (commented intentionally):
# # Verify accession metadata and processed table availability.
# # Create local directory:
# # mkdir -p data/phase2b/processed_tables/wallen_prjna834801
# # Manually place:
# #   sample_metadata.tsv
# #   species_abundance.tsv
# #   gene_family_abundance.tsv
# #   pathway_abundance.tsv
# #   run_accession_mapping.tsv

echo "[DONE] No active download executed."
