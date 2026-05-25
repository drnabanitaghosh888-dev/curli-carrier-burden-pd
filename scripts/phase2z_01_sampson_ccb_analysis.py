#!/usr/bin/env python3
from __future__ import annotations

import math
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = Path(".")
META_FILE = ROOT / "data/phase2t/manual_sources/duru_asapmac/asapmac_sampleMetadata_from_source.tsv"
CANDIDATE_FILES = [
    ROOT / "metadata/curli_candidate_taxa_from_phase1.tsv",
    ROOT / "data/000_phase2o_curli_reference/000_curli_candidate_species.tsv",
]
OUTDIR = ROOT / "data/phase2z/reports"
OUTDIR.mkdir(parents=True, exist_ok=True)

REMOTE_RELATIVE_ABUNDANCE = "https://huggingface.co/datasets/waldronlab/metagenomics_mac/resolve/main/relative_abundance_uuid.parquet"
COHORT = "SampsonTR_2025"


def norm_taxon(x: str) -> str:
    x = str(x or "").strip()
    x = re.sub(r"\[[^\]]+\]", "", x)
    x = re.sub(r"\b[skpcofg]__", "", x)
    x = x.replace("|", " ")
    x = x.replace("_", " ")
    x = re.sub(r"\s+", " ", x)
    return x.strip().lower()


def choose_candidate_file() -> Path:
    for p in CANDIDATE_FILES:
        if p.exists():
            return p
    raise FileNotFoundError("No curated curli candidate taxa file found.")


def candidate_taxon_column(df: pd.DataFrame) -> str:
    preferred = [
        "candidate_taxon", "candidate_species", "species", "taxon",
        "matched_candidate_taxon", "curli_candidate_taxon", "scientific_name",
        "organism", "name"
    ]
    for c in preferred:
        if c in df.columns:
            return c
    for c in df.columns:
        cl = c.lower()
        if "species" in cl or "taxon" in cl or "organism" in cl or "name" in cl:
            return c
    raise ValueError(f"Could not identify candidate taxon column in {list(df.columns)}")


def infer_weight(row: pd.Series) -> float:
    for c in row.index:
        if c.lower() == "weight":
            try:
                return float(row[c])
            except Exception:
                pass
    confidence = ""
    for c in row.index:
        if "confidence" in c.lower():
            confidence = str(row[c]).lower()
            break
    if "high" in confidence:
        return 1.0
    if "medium" in confidence:
        return 0.5
    if "low" in confidence:
        return 0.25
    return 0.5


def load_candidates() -> pd.DataFrame:
    p = choose_candidate_file()
    cand = pd.read_csv(p, sep="\t", dtype=str)
    tax_col = candidate_taxon_column(cand)
    cand = cand.copy()
    cand["candidate_taxon_raw"] = cand[tax_col].astype(str)
    cand["candidate_taxon_norm"] = cand["candidate_taxon_raw"].map(norm_taxon)
    cand["candidate_genus_norm"] = cand["candidate_taxon_norm"].str.split().str[0]
    cand["weight"] = cand.apply(infer_weight, axis=1)
    cand["confidence"] = [
        "high" if w >= 1 else "medium" if w >= 0.5 else "low"
        for w in cand["weight"]
    ]
    cand = cand[cand["candidate_taxon_norm"].str.len() > 0].drop_duplicates("candidate_taxon_norm")
    return cand


def load_sampson_metadata() -> pd.DataFrame:
    meta = pd.read_csv(META_FILE, sep="\t", dtype=str, low_memory=False)

    cohort = meta[meta["study_name"].astype(str).str.fullmatch(COHORT, case=False, na=False)].copy()
    if cohort.empty:
        cohort = meta[meta["study_name"].astype(str).str.contains("SampsonTR_2025|SampsonTR|Sampson", case=False, regex=True, na=False)].copy()

    if cohort.empty:
        raise ValueError("No SampsonTR_2025 rows found in ASAP-MAC sampleMetadata.")

    # Phase 2W already classified SampsonTR_2025 as fecal material.
    # In the exported sampleMetadata, body_site is missing for this cohort,
    # so body_site absence must not remove otherwise valid Sampson rows.
    cohort = cohort[cohort["uuid"].notna() & (cohort["uuid"].astype(str).str.strip() != "")].copy()

    ctrl = cohort.get("control", pd.Series(index=cohort.index, dtype=str)).astype(str).str.strip().str.lower()
    targ = cohort.get("target_condition", pd.Series(index=cohort.index, dtype=str)).astype(str).str.strip().str.lower()

    phenotype = np.where(
        targ.str.contains("parkinson", na=False) & ctrl.eq("case"),
        "PD",
        np.where(
            targ.str.contains("parkinson", na=False) & ctrl.str.contains("control", na=False),
            "HC",
            "UNMAPPED",
        ),
    )
    cohort["phenotype_CCB"] = phenotype
    cohort = cohort[cohort["phenotype_CCB"].isin(["PD", "HC"])].copy()

    out_cols = []
    for c in ["uuid", "study_name", "phenotype_CCB", "body_site", "target_condition", "control", "disease", "country", "age", "sex", "sample_id"]:
        if c in cohort.columns:
            out_cols.append(c)

    manifest = cohort[out_cols].drop_duplicates("uuid").copy()
    manifest.to_csv(OUTDIR / "phase2z_sampson_uuid_phenotype.tsv", sep="\t", index=False)
    return manifest


def fetch_relative_abundance(uuids: pd.DataFrame) -> pd.DataFrame:
    uuid_file = OUTDIR / "phase2z_sampson_uuid_phenotype.tsv"

    con = duckdb.connect()
    con.execute("INSTALL httpfs;")
    con.execute("LOAD httpfs;")

    query = f"""
    SELECT ra.*
    FROM read_parquet('{REMOTE_RELATIVE_ABUNDANCE}') AS ra
    INNER JOIN read_csv_auto('{uuid_file}', delim='\\t', header=true) AS m
    ON ra.uuid = m.uuid
    """

    rel = con.execute(query).fetchdf()
    rel.to_csv(OUTDIR / "phase2z_sampson_relative_abundance_long.tsv", sep="\t", index=False)

    if "clade_name_species" in rel.columns:
        sp = rel[rel["clade_name_species"].notna() & (rel["clade_name_species"].astype(str).str.strip() != "")].copy()
    elif "clade_name" in rel.columns:
        sp = rel[rel["clade_name"].astype(str).str.contains("s__", na=False)].copy()
    else:
        sp = rel.copy()

    sp.to_csv(OUTDIR / "phase2z_sampson_relative_abundance_species_like.tsv", sep="\t", index=False)
    return sp


def taxon_label(row: pd.Series) -> str:
    for c in ["clade_name_species", "clade_name_terminal", "clade_name"]:
        if c in row.index and pd.notna(row[c]) and str(row[c]).strip():
            return str(row[c])
    return ""


def match_curli_taxa(sp: pd.DataFrame, cand: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sp = sp.copy()
    sp["sampson_taxon"] = sp.apply(taxon_label, axis=1)
    sp["taxon_norm"] = sp["sampson_taxon"].map(norm_taxon)
    sp["genus_norm"] = sp["taxon_norm"].str.split().str[0]

    exact_map = dict(zip(cand["candidate_taxon_norm"], cand["candidate_taxon_raw"]))
    weight_map = dict(zip(cand["candidate_taxon_norm"], cand["weight"]))
    conf_map = dict(zip(cand["candidate_taxon_norm"], cand["confidence"]))

    genus_best = {}
    for _, r in cand.iterrows():
        g = r["candidate_genus_norm"]
        if g and g not in genus_best:
            genus_best[g] = r

    matched_rows = []
    for idx, r in sp.iterrows():
        tn = r["taxon_norm"]
        g = r["genus_norm"]
        if tn in exact_map:
            matched_rows.append((idx, exact_map[tn], "species_exact", conf_map[tn], weight_map[tn]))
        elif g in genus_best:
            cr = genus_best[g]
            matched_rows.append((idx, cr["candidate_taxon_raw"], "genus_fallback", cr["confidence"], cr["weight"]))

    if not matched_rows:
        return sp.iloc[0:0].copy(), pd.DataFrame()

    midx = [x[0] for x in matched_rows]
    matched = sp.loc[midx].copy()
    matched["matched_candidate_taxon"] = [x[1] for x in matched_rows]
    matched["match_level"] = [x[2] for x in matched_rows]
    matched["confidence"] = [x[3] for x in matched_rows]
    matched["weight"] = [float(x[4]) for x in matched_rows]
    matched["weighted_abundance"] = matched["relative_abundance"].astype(float) * matched["weight"]

    taxa_report = (
        matched.groupby(["sampson_taxon", "matched_candidate_taxon", "match_level", "confidence", "weight"], dropna=False)
        .agg(
            n_samples_present=("uuid", "nunique"),
            mean_relative_abundance=("relative_abundance", "mean"),
            max_relative_abundance=("relative_abundance", "max"),
        )
        .reset_index()
        .sort_values(["n_samples_present", "mean_relative_abundance"], ascending=False)
    )
    taxa_report.to_csv(OUTDIR / "phase2z_sampson_matched_curli_taxa.tsv", sep="\t", index=False)
    return matched, taxa_report


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x)
    y = np.asarray(y)
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    gt = sum((xi > y).sum() for xi in x)
    lt = sum((xi < y).sum() for xi in x)
    return (gt - lt) / (len(x) * len(y))


def analyze(sample_meta: pd.DataFrame, matched: pd.DataFrame, n_taxa: int) -> tuple[str, str, pd.DataFrame]:
    if sample_meta.empty or matched.empty:
        status = "DATA_NOT_READY"
        decision = "DATA_NOT_READY"
        summary = pd.DataFrame([{
            "analysis_set": "SampsonTR_2025_processed_species_like",
            "n_samples": len(sample_meta),
            "n_PD": int((sample_meta["phenotype_CCB"] == "PD").sum()) if "phenotype_CCB" in sample_meta else 0,
            "n_HC": int((sample_meta["phenotype_CCB"] == "HC").sum()) if "phenotype_CCB" in sample_meta else 0,
            "n_matched_curli_taxa": n_taxa,
        }])
        return status, decision, summary

    ccb = matched.groupby("uuid", as_index=False).agg(
        CCB=("weighted_abundance", "sum"),
        n_matched_curli_taxa_present=("matched_candidate_taxon", "nunique"),
    )

    samples = sample_meta[["uuid", "phenotype_CCB"]].drop_duplicates("uuid").merge(ccb, on="uuid", how="left")
    samples["CCB"] = samples["CCB"].fillna(0.0)
    samples["n_matched_curli_taxa_present"] = samples["n_matched_curli_taxa_present"].fillna(0).astype(int)
    samples["log1p_CCB"] = np.log1p(samples["CCB"])
    samples.to_csv(OUTDIR / "phase2z_sampson_sample_ccb.tsv", sep="\t", index=False)

    pd_vals = samples.loc[samples["phenotype_CCB"] == "PD", "CCB"].astype(float).to_numpy()
    hc_vals = samples.loc[samples["phenotype_CCB"] == "HC", "CCB"].astype(float).to_numpy()

    if len(pd_vals) < 2 or len(hc_vals) < 2:
        status = "DATA_NOT_READY"
        decision = "DATA_NOT_READY"
        u_stat = pval = delta = float("nan")
    else:
        u_stat, pval = mannwhitneyu(pd_vals, hc_vals, alternative="two-sided")
        delta = cliffs_delta(pd_vals, hc_vals)
        pd_mean = float(np.mean(pd_vals))
        hc_mean = float(np.mean(hc_vals))
        pd_med = float(np.median(pd_vals))
        hc_med = float(np.median(hc_vals))
        if pd_mean > hc_mean and pd_med > hc_med and pval < 0.05:
            decision = "STRONG_SUPPORT"
        elif pd_mean > hc_mean and pd_med > hc_med:
            decision = "DIRECTIONALLY_CONSISTENT_NOT_SIGNIFICANT"
        else:
            decision = "MIXED_OR_NO_SUPPORT"
        status = "PASS"

    summary = pd.DataFrame([{
        "analysis_set": "SampsonTR_2025_processed_species_like",
        "n_samples": len(samples),
        "n_PD": int((samples["phenotype_CCB"] == "PD").sum()),
        "n_HC": int((samples["phenotype_CCB"] == "HC").sum()),
        "n_matched_curli_taxa": n_taxa,
        "PD_mean_CCB": float(np.mean(pd_vals)) if len(pd_vals) else math.nan,
        "HC_mean_CCB": float(np.mean(hc_vals)) if len(hc_vals) else math.nan,
        "PD_median_CCB": float(np.median(pd_vals)) if len(pd_vals) else math.nan,
        "HC_median_CCB": float(np.median(hc_vals)) if len(hc_vals) else math.nan,
        "effect_direction_mean": "PD>HC" if len(pd_vals) and len(hc_vals) and np.mean(pd_vals) > np.mean(hc_vals) else "HC>=PD",
        "effect_direction_median": "PD>HC" if len(pd_vals) and len(hc_vals) and np.median(pd_vals) > np.median(hc_vals) else "HC>=PD",
        "mannwhitney_U": float(u_stat),
        "mannwhitney_p": float(pval),
        "cliffs_delta": float(delta),
        "decision": decision,
    }])
    summary.to_csv(OUTDIR / "phase2z_sampson_ccb_test_summary.tsv", sep="\t", index=False)
    return status, decision, summary


def write_reports(status: str, decision: str, summary: pd.DataFrame, n_rel_rows: int, n_species_rows: int):
    row = summary.iloc[0].to_dict()
    (OUTDIR / "phase2z_sampson_ccb_status.txt").write_text(status + "\n")

    lines = [
        "PHASE2Z_SAMPSONTR_2025_CCB_ANALYSIS",
        f"status={status}",
        f"decision={decision}",
        "cohort=SampsonTR_2025",
        "selection_rule=largest_not_yet_analyzed_eligible_fecal_PD_control_ASAP_MAC_cohort_from_Phase_2W",
        f"metadata_file={META_FILE}",
        f"processed_abundance_resource={REMOTE_RELATIVE_ABUNDANCE}",
        "analysis_type=processed_relative_abundance_only",
        "sequence_level_processing_performed=no",
        "causality_claimed=no",
        "clinical_marker_claimed=no",
    ]
    for k, v in row.items():
        lines.append(f"{k}={v}")
    lines += [
        f"n_relative_abundance_rows={n_rel_rows}",
        f"n_species_like_relative_abundance_rows={n_species_rows}",
        "interpretation=SampsonTR_2025 was analyzed using the locked CCB workflow without raw-read processing or method changes after cohort selection.",
    ]
    (OUTDIR / "phase2z_sampson_final_decision_report.txt").write_text("\n".join(lines) + "\n")
    print(f"STATUS={status}")
    print(f"decision={decision}")


def main():
    sample_meta = load_sampson_metadata()
    rel_sp = fetch_relative_abundance(sample_meta)
    cand = load_candidates()
    matched, taxa_report = match_curli_taxa(rel_sp, cand)
    status, decision, summary = analyze(sample_meta, matched, len(taxa_report))
    write_reports(status, decision, summary, n_rel_rows=-1, n_species_rows=len(rel_sp))


if __name__ == "__main__":
    main()
