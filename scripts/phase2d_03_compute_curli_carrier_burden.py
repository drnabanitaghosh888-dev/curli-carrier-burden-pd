#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd
import yaml


WEIGHTS = {"high": 1.0, "medium": 0.5, "low": 0.25, "unknown": 0.25}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute Curli Carrier Burden from processed species abundance.")
    p.add_argument("--execute", action="store_true", help="Write burden output table.")
    return p.parse_args()


def _load_cfg(root: Path) -> dict:
    with (root / "config/phase2d_config.yml").open("r", encoding="utf-8") as h:
        return yaml.safe_load(h)


def _to_weight(value: str) -> float:
    s = str(value).strip().lower()
    if s in WEIGHTS:
        return WEIGHTS[s]
    try:
        return float(s)
    except ValueError:
        return WEIGHTS["unknown"]


def build_primary_species_set(matches: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    required_cols = {
        "candidate_taxon_name",
        "curli_candidate_confidence",
        "metaphlan_clade_name",
        "metaphlan_species_name",
        "matching_method",
        "include_for_primary_burden",
    }
    missing = [c for c in required_cols if c not in matches.columns]
    if missing:
        raise ValueError(f"matches table missing required columns: {missing}")

    base = matches.copy()
    base["matching_method"] = base["matching_method"].astype(str)
    base["include_for_primary_burden"] = base["include_for_primary_burden"].astype(str).str.lower()
    base["candidate_taxon_name"] = base["candidate_taxon_name"].astype(str)
    base["curli_candidate_confidence"] = base["curli_candidate_confidence"].astype(str).str.lower()
    base["confidence_weight"] = base["curli_candidate_confidence"].map(_to_weight)

    pre_primary = base[base["include_for_primary_burden"] == "yes"].copy()
    pre_primary = pre_primary[pre_primary["matching_method"].isin(["exact_species_match", "binomial_species_match"])].copy()
    pre_primary = pre_primary[pre_primary["metaphlan_clade_name"].astype(str) != ""].copy()

    def _agg(group: pd.DataFrame) -> pd.Series:
        idx = group["confidence_weight"].astype(float).idxmax()
        best_conf = str(group.loc[idx, "curli_candidate_confidence"])
        collapsed = sorted(set(group["candidate_taxon_name"].tolist()))
        return pd.Series(
            {
                "metaphlan_species_name": group["metaphlan_species_name"].iloc[0],
                "confidence_weight": float(group["confidence_weight"].max()),
                "best_curli_candidate_confidence": best_conf,
                "n_candidate_taxa_collapsed": int(group.shape[0]),
                "collapsed_candidate_taxa": ";".join(collapsed),
            }
        )

    species_set = (
        pre_primary.groupby("metaphlan_clade_name", as_index=False)
        .apply(_agg)
        .reset_index(drop=True)
    )
    species_set = species_set[
        [
            "metaphlan_clade_name",
            "metaphlan_species_name",
            "confidence_weight",
            "best_curli_candidate_confidence",
            "n_candidate_taxa_collapsed",
            "collapsed_candidate_taxa",
        ]
    ].copy()

    stats = {
        "primary_candidate_rows_before_dedup": int(pre_primary.shape[0]),
        "unique_burden_species_after_dedup": int(species_set.shape[0]),
        "duplicated_candidate_rows_collapsed": int(pre_primary.shape[0] - species_set.shape[0]),
        "exact_or_binomial_rows_used": int(pre_primary.shape[0]),
        "genus_fallback_rows_excluded": int((base["matching_method"] == "genus_fallback_exploratory").sum()),
        "no_match_rows_excluded": int((base["matching_method"] == "no_match").sum()),
    }
    return species_set, stats


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg = _load_cfg(root)
    species_path = root / cfg["input_files"]["species_abundance"]
    matches_path = root / cfg["output_dirs"]["results"] / "wallen_curli_candidate_species_matches.tsv"
    meta_path = root / cfg["output_dirs"]["results"] / "wallen_analysis_metadata.tsv"

    out_results = root / cfg["output_dirs"]["results"]
    out_reports = root / cfg["output_dirs"]["reports"]
    out_results.mkdir(parents=True, exist_ok=True)
    out_reports.mkdir(parents=True, exist_ok=True)
    out_species_set = out_results / "wallen_curli_burden_species_set.tsv"
    out_burden = out_results / "wallen_curli_carrier_burden.tsv"
    out_summary = out_reports / "phase2d_curli_carrier_burden_summary.tsv"
    out_status = out_reports / "phase2d_curli_carrier_burden_status.txt"

    missing = [str(p) for p in [species_path, matches_path, meta_path] if not p.exists()]
    if missing:
        print("[ERROR] Missing required input(s):")
        for m in missing:
            print(m)
        return 1

    sp = pd.read_csv(species_path, sep="\t")
    mt = pd.read_csv(matches_path, sep="\t", dtype=str).fillna("")
    md = pd.read_csv(meta_path, sep="\t", dtype=str).fillna("")
    sample_col = cfg["sample_id_column"]

    species_set, stats = build_primary_species_set(mt)
    if species_set.empty:
        with out_summary.open("w", encoding="utf-8", newline="") as h:
            w = csv.DictWriter(h, fieldnames=["metric", "value"], delimiter="\t")
            w.writeheader()
            w.writerow({"metric": "status", "value": "FAIL"})
            for k, v in stats.items():
                w.writerow({"metric": k, "value": str(v)})
            w.writerow({"metric": "interpretation_note", "value": "no primary burden-eligible species; cannot compute burden"})
        out_status.write_text("FAIL\n", encoding="utf-8")
        print(f"[DONE] {out_summary}")
        print(f"[DONE] {out_status}")
        return 1

    weight_map = dict(zip(species_set["metaphlan_clade_name"], species_set["confidence_weight"]))

    sp = sp.copy()
    sp = sp[sp["clade_name"].astype(str).isin(weight_map)].copy()
    if sp.empty:
        print("[ERROR] No matched species found for burden computation.")
        return 1

    sample_cols = [c for c in sp.columns if c != "clade_name"]
    for c in sample_cols:
        sp[c] = pd.to_numeric(sp[c], errors="coerce").fillna(0.0)
    sp["weight"] = sp["clade_name"].map(weight_map).astype(float)

    weighted = sp[sample_cols].multiply(sp["weight"], axis=0)
    burden = weighted.sum(axis=0).rename("CurliCarrierBurden").reset_index().rename(columns={"index": sample_col})
    merged = md[[sample_col, cfg["case_control_column"], "PD_binary"]].merge(burden, on=sample_col, how="left")
    merged["CurliCarrierBurden"] = merged["CurliCarrierBurden"].fillna(0.0)

    summary_rows = [
        {"metric": "status", "value": "PASS"},
        {"metric": "primary_candidate_rows_before_dedup", "value": str(stats["primary_candidate_rows_before_dedup"])},
        {"metric": "unique_burden_species_after_dedup", "value": str(stats["unique_burden_species_after_dedup"])},
        {"metric": "duplicated_candidate_rows_collapsed", "value": str(stats["duplicated_candidate_rows_collapsed"])},
        {"metric": "exact_or_binomial_rows_used", "value": str(stats["exact_or_binomial_rows_used"])},
        {"metric": "genus_fallback_rows_excluded", "value": str(stats["genus_fallback_rows_excluded"])},
        {"metric": "no_match_rows_excluded", "value": str(stats["no_match_rows_excluded"])},
        {"metric": "n_samples", "value": str(len(merged))},
        {"metric": "burden_min", "value": f"{merged['CurliCarrierBurden'].min():.6f}"},
        {"metric": "burden_median", "value": f"{merged['CurliCarrierBurden'].median():.6f}"},
        {"metric": "burden_max", "value": f"{merged['CurliCarrierBurden'].max():.6f}"},
        {"metric": "interpretation_note", "value": "processed-table proxy only; not gene/operon proof"},
    ]
    with out_summary.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["metric", "value"], delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)

    if not args.execute:
        status = "PASS" if stats["exact_or_binomial_rows_used"] > 0 else "FAIL"
        out_status.write_text(status + "\n", encoding="utf-8")
        print("[INFO] Dry-run only. Use --execute to write wallen_curli_carrier_burden.tsv")
        print(f"[DONE] {out_summary}")
        print(f"[DONE] {out_status}")
        print(f"[STATUS] {status}")
        return 0

    species_set.to_csv(out_species_set, sep="\t", index=False)
    merged.to_csv(out_burden, sep="\t", index=False)
    out_status.write_text("PASS\n", encoding="utf-8")
    print(f"[DONE] {out_species_set}")
    print(f"[DONE] {out_burden}")
    print(f"[DONE] {out_summary}")
    print(f"[DONE] {out_status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
