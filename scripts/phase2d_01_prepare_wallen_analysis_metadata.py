#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd
import yaml


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare Wallen analysis metadata.")
    p.add_argument("--execute", action="store_true", help="Write output metadata table.")
    return p.parse_args()


def _load_cfg(root: Path) -> dict:
    with (root / "config/phase2d_config.yml").open("r", encoding="utf-8") as h:
        return yaml.safe_load(h)


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg = _load_cfg(root)
    sm_path = root / cfg["input_files"]["sample_metadata"]
    out_results = root / cfg["output_dirs"]["results"]
    out_reports = root / cfg["output_dirs"]["reports"]
    out_results.mkdir(parents=True, exist_ok=True)
    out_reports.mkdir(parents=True, exist_ok=True)
    out_meta = out_results / "wallen_analysis_metadata.tsv"
    out_summary = out_reports / "phase2d_wallen_metadata_summary.tsv"

    if not sm_path.exists():
        print(f"[ERROR] Missing sample metadata: {sm_path}")
        return 1

    sm = pd.read_csv(sm_path, sep="\t", dtype=str)
    sample_col = cfg["sample_id_column"]
    case_col = cfg["case_control_column"]

    wanted = [sample_col, case_col] + cfg["covariates"]["core"] + cfg["covariates"]["exploratory"]
    keep = [c for c in wanted if c in sm.columns]
    analysis = sm[keep].copy()
    if case_col in analysis.columns:
        analysis["PD_binary"] = analysis[case_col].map({cfg["case_label"]: 1, cfg["control_label"]: 0})

    rows = []
    for c in analysis.columns:
        missing_n = int(analysis[c].isna().sum())
        rows.append({"column": c, "missing_n": str(missing_n), "missing_fraction": f"{missing_n / max(len(analysis), 1):.4f}"})

    with out_summary.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["column", "missing_n", "missing_fraction"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    if not args.execute:
        print("[INFO] Dry-run only. Use --execute to write wallen_analysis_metadata.tsv")
        print(f"[DONE] {out_summary}")
        return 0

    analysis.to_csv(out_meta, sep="\t", index=False)
    print(f"[DONE] {out_meta}")
    print(f"[DONE] {out_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
