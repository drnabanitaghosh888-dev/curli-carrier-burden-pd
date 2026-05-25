#!/usr/bin/env bash
set -euo pipefail

if [ "${RUN_MANUAL_DOWNLOAD:-NO}" != "YES" ]; then
  echo "Manual download guard active. Set RUN_MANUAL_DOWNLOAD=YES to run."
  exit 1
fi

echo "[INFO] Replication processed-table manual template only."

# Placeholder commands only (commented intentionally):
# # Verify exact accession/source and license before any transfer.
# # Cohorts:
# #   palacios_prodromal_pd
# #   nishiwaki_multicountry_pd
# #   mao_central_china_pd
# #   integrated_us_multicohort_pd
# # Create directories under data/phase2b/processed_tables/<cohort>/
# # Manually place processed tables and metadata once verified.

echo "[DONE] No active download executed."
