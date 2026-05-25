#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


EXPECTED = [
    "sample_metadata.tsv",
    "species_abundance.tsv",
    "genus_abundance.tsv",
    "gene_family_abundance.tsv",
    "pathway_abundance.tsv",
    "run_accession_mapping.tsv",
]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    base = root / "data/phase2b/processed_tables/wallen_prjna834801"
    out = root / "data/phase2c/reports/phase2c_wallen_local_file_check.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for fname in EXPECTED:
        p = base / fname
        rows.append(
            {
                "expected_file": fname,
                "local_path": str(p.relative_to(root)),
                "exists": "yes" if p.exists() else "no",
            }
        )

    with out.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["expected_file", "local_path", "exists"], delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"[DONE] wrote: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
