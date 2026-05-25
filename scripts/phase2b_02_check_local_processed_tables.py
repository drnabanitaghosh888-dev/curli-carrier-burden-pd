#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


COHORTS = [
    "wallen_prjna834801",
    "palacios_prodromal_pd",
    "nishiwaki_multicountry_pd",
    "mao_central_china_pd",
    "integrated_us_multicohort_pd",
]
FILES = [
    "sample_metadata.tsv",
    "species_abundance.tsv",
    "gene_family_abundance.tsv",
    "pathway_abundance.tsv",
    "run_accession_mapping.tsv",
]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    base = root / "data/phase2b/processed_tables"
    out = root / "data/phase2b/reports/local_processed_table_inventory.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for cohort in COHORTS:
        for fname in FILES:
            p = base / cohort / fname
            rows.append(
                {
                    "cohort_id": cohort,
                    "table_file": fname,
                    "expected_path": str(p.relative_to(root)),
                    "exists": "yes" if p.exists() else "no",
                }
            )

    with out.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["cohort_id", "table_file", "expected_path", "exists"], delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"[DONE] wrote: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
