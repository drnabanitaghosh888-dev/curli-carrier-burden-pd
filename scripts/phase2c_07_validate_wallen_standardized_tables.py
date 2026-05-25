#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


BASE = "data/phase2b/processed_tables/wallen_prjna834801"


def exists(path: Path) -> bool:
    return path.exists() and path.is_file()


def table_columns(path: Path) -> list[str]:
    if not exists(path):
        return []
    df = pd.read_csv(path, sep="\t", nrows=0)
    return [str(c) for c in df.columns]


def sample_col_count(path: Path) -> int:
    cols = table_columns(path)
    if "clade_name" in cols:
        return max(0, len(cols) - 1)
    if "Gene Family" in cols:
        return max(0, len(cols) - 1)
    if "Pathway" in cols:
        return max(0, len(cols) - 1)
    return len(cols)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    base = root / BASE
    report_dir = root / "data/phase2c/reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = report_dir / "phase2c_wallen_standardized_table_validation.tsv"
    out_status = report_dir / "phase2c_wallen_standardized_table_status.txt"

    sample = base / "sample_metadata.tsv"
    species = base / "species_abundance.tsv"
    genefam = base / "gene_family_abundance.tsv"
    pathway = base / "pathway_abundance.tsv"
    runmap = base / "run_accession_mapping.tsv"

    rows = []
    status = "PASS"

    for label, p in [
        ("sample_metadata.tsv", sample),
        ("species_abundance.tsv", species),
        ("gene_family_abundance.tsv", genefam),
        ("pathway_abundance.tsv", pathway),
    ]:
        ok = exists(p)
        if not ok:
            status = "FAIL"
        rows.append(
            {
                "check": f"exists::{label}",
                "result": "pass" if ok else "fail",
                "details": str(p.relative_to(root)),
            }
        )

    scols = table_columns(sample)
    rows.append(
        {
            "check": "sample_metadata_has_sample_name",
            "result": "pass" if "sample_name" in scols else "fail",
            "details": ",".join(scols[:20]),
        }
    )
    if "sample_name" not in scols:
        status = "FAIL"
    rows.append(
        {
            "check": "sample_metadata_has_Case_status",
            "result": "pass" if "Case_status" in scols else "fail",
            "details": ",".join(scols[:20]),
        }
    )
    if "Case_status" not in scols:
        status = "FAIL"

    spcols = table_columns(species)
    rows.append(
        {
            "check": "species_abundance_has_clade_name",
            "result": "pass" if "clade_name" in spcols else "fail",
            "details": ",".join(spcols[:20]),
        }
    )
    if "clade_name" not in spcols:
        status = "FAIL"

    gfcols = table_columns(genefam)
    rows.append(
        {
            "check": "gene_family_abundance_has_Gene Family",
            "result": "pass" if "Gene Family" in gfcols else "fail",
            "details": ",".join(gfcols[:20]),
        }
    )
    if "Gene Family" not in gfcols:
        status = "FAIL"

    pwcols = table_columns(pathway)
    rows.append(
        {
            "check": "pathway_abundance_has_Pathway",
            "result": "pass" if "Pathway" in pwcols else "fail",
            "details": ",".join(pwcols[:20]),
        }
    )
    if "Pathway" not in pwcols:
        status = "FAIL"

    # broad compatibility: abundance sample columns roughly compatible with metadata rows
    try:
        n_meta = len(pd.read_csv(sample, sep="\t")) if exists(sample) else 0
    except Exception:
        n_meta = 0
    for label, p in [
        ("species_abundance", species),
        ("gene_family_abundance", genefam),
        ("pathway_abundance", pathway),
    ]:
        n_cols = sample_col_count(p)
        compatible = (n_meta == 0) or (n_cols == 0) or (abs(n_cols - n_meta) <= max(5, int(0.1 * max(n_meta, 1))))
        rows.append(
            {
                "check": f"sample_count_compatibility::{label}",
                "result": "pass" if compatible else "warn",
                "details": f"metadata_rows={n_meta}; abundance_sample_columns={n_cols}",
            }
        )

    runmap_exists = exists(runmap)
    rows.append(
        {
            "check": "run_accession_mapping_presence",
            "result": "pass" if runmap_exists else "warn",
            "details": "missing blocks Phase 3 but not Phase 2D",
        }
    )

    with out_tsv.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["check", "result", "details"], delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    if status == "PASS" and not runmap_exists:
        status = "PASS_WITH_WARNINGS"
    out_status.write_text(status + "\n", encoding="utf-8")
    print(f"[DONE] wrote: {out_tsv}")
    print(f"[DONE] wrote: {out_status}")
    print(f"[STATUS] {status}")
    return 0 if status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
