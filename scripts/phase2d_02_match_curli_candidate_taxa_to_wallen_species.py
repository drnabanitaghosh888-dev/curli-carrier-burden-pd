#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
import yaml

WEIGHTS = {"high": 1.0, "medium": 0.5, "low": 0.25, "unknown": 0.25}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Match curli candidate taxa to Wallen species rows.")
    p.add_argument("--execute", action="store_true", help="Write species matches table.")
    return p.parse_args()


def _load_cfg(root: Path) -> dict:
    with (root / "config/phase2d_config.yml").open("r", encoding="utf-8") as h:
        return yaml.safe_load(h)


def _normalize_text(x: str) -> str:
    s = str(x).strip().lower().replace("_", " ")
    s = re.sub(r"[^\w\s]", " ", s)
    s = re.sub(r"\s+", " ", s)
    return s


def _normalize_candidate_taxon(taxon_name: str) -> dict[str, str]:
    full = _normalize_text(taxon_name)
    tokens = full.split()
    binomial = " ".join(tokens[:2]) if len(tokens) >= 2 else full
    genus = tokens[0] if tokens else ""
    return {
        "candidate_full_normalized": full,
        "candidate_binomial_normalized": binomial,
        "candidate_genus": genus,
    }


def _parse_metaphlan_clade(clade_name: str) -> dict[str, str]:
    raw = str(clade_name)
    parts = raw.split("|")
    species_raw = ""
    for p in parts:
        if p.startswith("s__"):
            species_raw = p
            break

    species_name = species_raw.replace("s__", "").replace("_", " ").strip()
    species_normalized = _normalize_text(species_name)
    tokens = species_normalized.split()
    binomial = " ".join(tokens[:2]) if len(tokens) >= 2 else species_normalized
    genus = tokens[0] if tokens else ""
    return {
        "metaphlan_clade_name": raw,
        "metaphlan_species_raw": species_raw,
        "metaphlan_species_name": species_name,
        "metaphlan_species_normalized": species_normalized,
        "metaphlan_binomial_normalized": binomial,
        "metaphlan_genus_normalized": genus,
    }


def _confidence_weight(conf: str) -> float:
    c = _normalize_text(conf)
    return WEIGHTS.get(c, WEIGHTS["unknown"])


def _match_candidate_to_species(cand: dict[str, str], species_row: dict[str, str]) -> tuple[str, str]:
    if cand["candidate_full_normalized"] == species_row["metaphlan_species_normalized"]:
        return "exact_species_match", "yes"
    if cand["candidate_binomial_normalized"] and cand["candidate_binomial_normalized"] == species_row["metaphlan_binomial_normalized"]:
        return "binomial_species_match", "yes"
    if cand["candidate_genus"] and cand["candidate_genus"] == species_row["metaphlan_genus_normalized"]:
        return "genus_fallback_exploratory", "no"
    return "no_match", "no"


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg = _load_cfg(root)
    cand_path = root / cfg["input_files"]["curli_candidate_taxa"]
    sp_path = root / cfg["input_files"]["species_abundance"]
    out_results = root / cfg["output_dirs"]["results"]
    out_reports = root / cfg["output_dirs"]["reports"]
    out_results.mkdir(parents=True, exist_ok=True)
    out_reports.mkdir(parents=True, exist_ok=True)
    out_match = out_results / "wallen_curli_candidate_species_matches.tsv"
    out_summary = out_reports / "phase2d_curli_taxon_matching_summary.tsv"
    out_preview = out_reports / "phase2d_curli_taxon_matching_preview.tsv"
    out_status = out_reports / "phase2d_curli_taxon_matching_status.txt"

    if not cand_path.exists() or not sp_path.exists():
        print("[ERROR] Missing required input(s) for taxon matching.")
        return 1

    cand = pd.read_csv(cand_path, sep="\t", dtype=str).fillna("")
    sp = pd.read_csv(sp_path, sep="\t", dtype=str, usecols=["clade_name"]).fillna("")
    if "taxon_name" not in cand.columns:
        print("[ERROR] curli candidate taxa file missing taxon_name column.")
        return 1

    species_rows = [_parse_metaphlan_clade(x) for x in sp["clade_name"].tolist()]

    rows = []
    for _, r in cand.iterrows():
        cand_name = r["taxon_name"]
        cand_norm = _normalize_candidate_taxon(cand_name)
        conf = r.get("curli_candidate_confidence", "unknown")
        weight = _confidence_weight(conf)

        best_method = "no_match"
        best_row = {
            "metaphlan_clade_name": "",
            "metaphlan_species_name": "",
            "metaphlan_species_normalized": "",
            "metaphlan_binomial_normalized": "",
            "metaphlan_genus_normalized": "",
        }
        rank = {
            "exact_species_match": 3,
            "binomial_species_match": 2,
            "genus_fallback_exploratory": 1,
            "no_match": 0,
        }
        for srow in species_rows:
            method, include_primary = _match_candidate_to_species(cand_norm, srow)
            if rank[method] > rank[best_method]:
                best_method = method
                best_row = srow
                best_include = include_primary
                if best_method == "exact_species_match":
                    break

        if best_method == "no_match":
            best_include = "no"

        rows.append(
            {
                "candidate_taxon_name": cand_name,
                "candidate_full_normalized": cand_norm["candidate_full_normalized"],
                "candidate_binomial_normalized": cand_norm["candidate_binomial_normalized"],
                "candidate_genus": cand_norm["candidate_genus"],
                "curli_candidate_confidence": conf,
                "metaphlan_clade_name": best_row["metaphlan_clade_name"],
                "metaphlan_species_name": best_row["metaphlan_species_name"],
                "metaphlan_species_normalized": best_row["metaphlan_species_normalized"],
                "metaphlan_binomial_normalized": best_row["metaphlan_binomial_normalized"],
                "metaphlan_genus_normalized": best_row["metaphlan_genus_normalized"],
                "matching_method": best_method,
                "include_for_primary_burden": best_include,
                "confidence_weight": f"{weight:.2f}",
                "notes": "processed-table proxy only; no strain-level or operon-level inference",
            }
        )

    out_df = pd.DataFrame(
        rows,
        columns=[
            "candidate_taxon_name",
            "candidate_full_normalized",
            "candidate_binomial_normalized",
            "candidate_genus",
            "curli_candidate_confidence",
            "metaphlan_clade_name",
            "metaphlan_species_name",
            "metaphlan_species_normalized",
            "metaphlan_binomial_normalized",
            "metaphlan_genus_normalized",
            "matching_method",
            "include_for_primary_burden",
            "confidence_weight",
            "notes",
        ],
    )

    summary = out_df["matching_method"].value_counts().rename_axis("matching_method").reset_index(name="count")
    summary.to_csv(out_summary, sep="\t", index=False)

    non_no = out_df[out_df["matching_method"] != "no_match"]
    preview = non_no.head(50) if not non_no.empty else out_df[out_df["matching_method"] == "no_match"].head(50)
    preview.to_csv(out_preview, sep="\t", index=False)

    if ((out_df["matching_method"] == "exact_species_match") | (out_df["matching_method"] == "binomial_species_match")).any():
        status = "PASS"
    elif (out_df["matching_method"] == "genus_fallback_exploratory").any():
        status = "PASS_WITH_WARNINGS"
    else:
        status = "FAIL"
    out_status.write_text(status + "\n", encoding="utf-8")

    if not args.execute:
        print("[INFO] Dry-run only. Use --execute to write matches table.")
        print(f"[DONE] {out_summary}")
        print(f"[DONE] {out_preview}")
        print(f"[DONE] {out_status}")
        print(f"[STATUS] {status}")
        return 0

    out_df.to_csv(out_match, sep="\t", index=False)
    print(f"[DONE] {out_match}")
    print(f"[DONE] {out_summary}")
    print(f"[DONE] {out_preview}")
    print(f"[DONE] {out_status}")
    print(f"[STATUS] {status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
