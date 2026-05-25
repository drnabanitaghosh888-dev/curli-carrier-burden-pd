#!/usr/bin/env python3
"""Phase 2U Duru processed-table Curli Carrier Burden analysis.

This script is intended to be run manually. It uses existing local processed
relative-abundance rows and local phenotype metadata prepared during Phase 2T.
It does not retrieve data or perform sequence-level processing.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

META = ROOT / "data/phase2t/manual_sources/duru_asapmac/duruIC_2024_uuid_phenotype.tsv"
ABUNDANCE = (
    ROOT / "data/phase2t/manual_sources/duru_asapmac/"
    / "duruIC_2024_relative_abundance_species_like.tsv"
)
CANDIDATES = ROOT / "metadata/curli_candidate_taxa_from_phase1.tsv"

REPORT_DIR = ROOT / "data/phase2u/reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)

OUT_BURDEN = REPORT_DIR / "phase2u_duru_ccb_sample_burden.tsv"
OUT_MATCHED_TAXA = REPORT_DIR / "phase2u_duru_ccb_matched_taxa.tsv"
OUT_TEST = REPORT_DIR / "phase2u_duru_ccb_test_summary.tsv"
OUT_STATUS = REPORT_DIR / "phase2u_duru_ccb_status.txt"
OUT_DECISION = REPORT_DIR / "phase2u_duru_final_decision_report.txt"

WEIGHTS = {
    "high": 1.0,
    "medium": 0.5,
    "low": 0.25,
}


def normalize_name(value: object) -> str:
    text = "" if pd.isna(value) else str(value)
    text = text.replace("_", " ")
    text = re.sub(r"\[[^\]]+\]", " ", text)
    text = re.sub(r"\b(strain|subsp|subspecies|serovar|serotype|group|complex)\b", " ", text, flags=re.I)
    text = re.sub(r"[^A-Za-z0-9 ]+", " ", text)
    return re.sub(r"\s+", " ", text).strip().lower()


def binomial(value: object) -> str:
    norm = normalize_name(value)
    parts = norm.split()
    if len(parts) >= 2:
        return " ".join(parts[:2])
    return norm


def clean_species_from_clade(row: pd.Series) -> str:
    species = row.get("clade_name_species", "")
    if pd.notna(species) and str(species).strip():
        label = str(species).strip()
    else:
        label = str(row.get("clade_name", "")).split("|")[-1]
    label = re.sub(r"^[a-z]__", "", label)
    return re.sub(r"\s+", " ", label.replace("_", " ")).strip()


def load_candidates() -> pd.DataFrame:
    df = pd.read_csv(CANDIDATES, sep="\t", dtype=str)
    if "taxon_name" not in df.columns:
        raise RuntimeError("Candidate table is missing taxon_name")
    confidence_col = "curli_candidate_confidence"
    if confidence_col not in df.columns:
        df[confidence_col] = "high"
    out = df[["taxon_name", confidence_col]].copy()
    out = out.rename(columns={confidence_col: "confidence"})
    out["confidence"] = out["confidence"].str.lower().map(lambda x: x if x in WEIGHTS else "high")
    out["weight"] = out["confidence"].map(WEIGHTS).astype(float)
    out["candidate_full_norm"] = out["taxon_name"].map(normalize_name)
    out["candidate_binomial_norm"] = out["taxon_name"].map(binomial)
    out = out[out["candidate_binomial_norm"].str.contains(" ", regex=False, na=False)].copy()
    return out.drop_duplicates(["candidate_binomial_norm", "confidence"])


def candidate_lookup(candidates: pd.DataFrame) -> dict[str, dict[str, Any]]:
    lookup: dict[str, dict[str, Any]] = {}
    for _, row in candidates.iterrows():
        key = row["candidate_binomial_norm"]
        record = {
            "candidate_taxon": row["taxon_name"],
            "candidate_binomial_norm": key,
            "confidence": row["confidence"],
            "weight": float(row["weight"]),
        }
        if key not in lookup or record["weight"] > lookup[key]["weight"]:
            lookup[key] = record
    return lookup


def collect_matched_taxa(candidates: pd.DataFrame) -> pd.DataFrame:
    lookup = candidate_lookup(candidates)
    keep_cols = [
        "clade_name",
        "clade_name_species",
        "clade_name_genus",
        "uuid",
        "relative_abundance",
    ]
    matched_records: dict[str, dict[str, Any]] = {}

    for chunk in pd.read_csv(ABUNDANCE, sep="\t", dtype=str, chunksize=100_000, usecols=lambda c: c in keep_cols):
        if "clade_name" not in chunk.columns or "relative_abundance" not in chunk.columns:
            raise RuntimeError("Abundance table is missing clade_name or relative_abundance")
        for _, row in chunk.drop_duplicates("clade_name").iterrows():
            clade = str(row["clade_name"])
            if clade in matched_records:
                continue
            species_name = clean_species_from_clade(row)
            species_norm = normalize_name(species_name)
            binom = binomial(species_name)
            if binom in lookup:
                candidate = lookup[binom]
                matched_records[clade] = {
                    "clade_name": clade,
                    "metaphlan_species_name": species_name,
                    "metaphlan_species_normalized": species_norm,
                    "metaphlan_binomial_normalized": binom,
                    "matched_candidate_taxon": candidate["candidate_taxon"],
                    "candidate_binomial_normalized": candidate["candidate_binomial_norm"],
                    "confidence": candidate["confidence"],
                    "confidence_weight": candidate["weight"],
                    "matching_method": "binomial_species_match",
                    "include_for_primary_burden": "yes",
                }
    return pd.DataFrame(matched_records.values())


def compute_burden(matched: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    weights = matched.set_index("clade_name")["confidence_weight"].astype(float).to_dict()
    wanted = set(weights)
    totals: dict[str, float] = {uuid: 0.0 for uuid in meta["uuid"].astype(str)}

    for chunk in pd.read_csv(
        ABUNDANCE,
        sep="\t",
        dtype=str,
        chunksize=100_000,
        usecols=["clade_name", "uuid", "relative_abundance"],
    ):
        chunk = chunk[chunk["clade_name"].isin(wanted)]
        if chunk.empty:
            continue
        chunk["relative_abundance"] = pd.to_numeric(chunk["relative_abundance"], errors="coerce").fillna(0.0)
        chunk["weighted_abundance"] = chunk["relative_abundance"] * chunk["clade_name"].map(weights).astype(float)
        grouped = chunk.groupby("uuid")["weighted_abundance"].sum()
        for uuid, value in grouped.items():
            uuid = str(uuid)
            if uuid in totals:
                totals[uuid] += float(value)

    burden = meta.copy()
    burden["CurliCarrierBurden"] = burden["uuid"].map(totals).fillna(0.0).astype(float)
    burden["CurliCarrierBurden_log1p"] = burden["CurliCarrierBurden"].map(math.log1p)
    burden["CurliCarrierPresence"] = (burden["CurliCarrierBurden"] > 0).astype(int)
    return burden


def mann_whitney(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    try:
        from scipy.stats import mannwhitneyu  # type: ignore

        result = mannwhitneyu(x, y, alternative="two-sided")
        return float(result.statistic), float(result.pvalue)
    except Exception:
        return float("nan"), float("nan")


def cliffs_delta(x: pd.Series, y: pd.Series) -> float:
    xv = pd.to_numeric(x, errors="coerce").dropna().tolist()
    yv = pd.to_numeric(y, errors="coerce").dropna().tolist()
    if not xv or not yv:
        return float("nan")
    gt = sum(1 for a in xv for b in yv if a > b)
    lt = sum(1 for a in xv for b in yv if a < b)
    return float((gt - lt) / (len(xv) * len(yv)))


def summarize_tests(burden: pd.DataFrame, n_matched: int) -> pd.DataFrame:
    pd_vals = burden.loc[burden["phenotype"] == "PD", "CurliCarrierBurden"]
    hc_vals = burden.loc[burden["phenotype"] == "HC", "CurliCarrierBurden"]
    u_stat, p_value = mann_whitney(pd_vals, hc_vals)
    return pd.DataFrame(
        [
            {
                "analysis_set": "DuruIC_2024_processed_species_like",
                "n_samples": len(burden),
                "n_PD": int((burden["phenotype"] == "PD").sum()),
                "n_HC": int((burden["phenotype"] == "HC").sum()),
                "n_matched_curli_taxa": int(n_matched),
                "PD_mean_CCB": float(pd_vals.mean()),
                "HC_mean_CCB": float(hc_vals.mean()),
                "PD_median_CCB": float(pd_vals.median()),
                "HC_median_CCB": float(hc_vals.median()),
                "effect_direction_mean": "PD>HC" if pd_vals.mean() > hc_vals.mean() else "PD<=HC",
                "effect_direction_median": "PD>HC" if pd_vals.median() > hc_vals.median() else "PD<=HC",
                "mannwhitney_U": u_stat,
                "mannwhitney_p": p_value,
                "cliffs_delta": cliffs_delta(pd_vals, hc_vals),
            }
        ]
    )


def decide(test: pd.Series) -> str:
    mean_support = test["effect_direction_mean"] == "PD>HC"
    median_support = test["effect_direction_median"] == "PD>HC"
    p_value = pd.to_numeric(pd.Series([test["mannwhitney_p"]]), errors="coerce").iloc[0]
    if mean_support and median_support and pd.notna(p_value) and float(p_value) < 0.05:
        return "SUPPORTIVE_REPLICATION"
    if mean_support or median_support:
        return "DIRECTIONALLY_CONSISTENT_NOT_SIGNIFICANT"
    return "NO_SUPPORT_OR_OPPOSITE_DIRECTION"


def write_decision(status: str, decision: str, test: pd.DataFrame | None = None, notes: str = "") -> None:
    lines = [
        "PHASE2U_DURU_CCB_ANALYSIS",
        f"status={status}",
        f"decision={decision}",
        f"metadata_file={META.relative_to(ROOT)}",
        f"abundance_file={ABUNDANCE.relative_to(ROOT)}",
        f"candidate_file={CANDIDATES.relative_to(ROOT)}",
        "input_type=local_processed_relative_abundance",
        "sequence_level_processing_performed=no",
    ]
    if test is not None and not test.empty:
        row = test.iloc[0].to_dict()
        lines.extend(f"{key}={value}" for key, value in row.items())
    if notes:
        lines.append(f"notes={notes}")
    OUT_STATUS.write_text(status + "\n", encoding="utf-8")
    OUT_DECISION.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    missing = [path for path in [META, ABUNDANCE, CANDIDATES] if not path.exists()]
    if missing:
        write_decision(
            "INPUT_MISSING",
            "NOT_ANALYZED",
            notes="missing_files:" + ";".join(str(path.relative_to(ROOT)) for path in missing),
        )
        print("STATUS=INPUT_MISSING")
        return

    meta = pd.read_csv(META, sep="\t", dtype=str)
    if not {"uuid", "phenotype_CCB"}.issubset(meta.columns):
        write_decision("FAIL", "NOT_ANALYZED", notes="metadata missing uuid or phenotype_CCB")
        print("STATUS=FAIL")
        return
    meta = meta.rename(columns={"phenotype_CCB": "phenotype"}).copy()
    meta = meta[meta["phenotype"].isin(["PD", "HC"])].drop_duplicates("uuid")

    candidates = load_candidates()
    matched = collect_matched_taxa(candidates)
    if matched.empty:
        matched.to_csv(OUT_MATCHED_TAXA, sep="\t", index=False)
        write_decision("NO_CURLI_TAXA_MATCHED", "NOT_ANALYZED", notes="no species-like processed taxa matched candidate binomials")
        print("STATUS=NO_CURLI_TAXA_MATCHED")
        return

    burden = compute_burden(matched, meta)
    test = summarize_tests(burden, len(matched))
    decision = decide(test.iloc[0])

    matched.to_csv(OUT_MATCHED_TAXA, sep="\t", index=False)
    burden.to_csv(OUT_BURDEN, sep="\t", index=False)
    test.to_csv(OUT_TEST, sep="\t", index=False)
    write_decision("PASS", decision, test)
    print("STATUS=PASS")
    print(f"decision={decision}")
    print(f"n_samples={len(burden)}")
    print(f"n_matched_curli_taxa={len(matched)}")


if __name__ == "__main__":
    main()
