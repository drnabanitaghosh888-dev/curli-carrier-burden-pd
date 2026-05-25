#!/usr/bin/env python3
"""Phase 2X NishiwakiH_2024 processed-abundance CCB analysis.

Run manually only. This script uses exported ASAP-MAC sample metadata and the
processed relative_abundance_uuid.parquet resource. It does not perform
sequence-level processing.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
METADATA = ROOT / "data/phase2t/manual_sources/duru_asapmac/asapmac_sampleMetadata_from_source.tsv"
HF_PARQUET_PATHS = ROOT / "data/phase2t/manual_sources/duru_asapmac/hf_parquet_paths.txt"
CURATED_CANDIDATES = ROOT / "metadata/curli_candidate_taxa_from_phase1.tsv"
FALLBACK_CANDIDATES = ROOT / "data/000_phase2o_curli_reference/000_curli_candidate_species.tsv"

RELATIVE_ABUNDANCE_PARQUET = (
    "https://huggingface.co/datasets/waldronlab/metagenomics_mac/resolve/main/"
    "relative_abundance_uuid.parquet"
)

REPORT_DIR = ROOT / "data/phase2x/reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)

OUT_MANIFEST = REPORT_DIR / "phase2x_nishiwaki_uuid_phenotype.tsv"
OUT_LONG = REPORT_DIR / "phase2x_nishiwaki_relative_abundance_long.tsv"
OUT_SPECIES_LIKE = REPORT_DIR / "phase2x_nishiwaki_relative_abundance_species_like.tsv"
OUT_SAMPLE_CCB = REPORT_DIR / "phase2x_nishiwaki_sample_ccb.tsv"
OUT_MATCHED_TAXA = REPORT_DIR / "phase2x_nishiwaki_matched_curli_taxa.tsv"
OUT_TEST_SUMMARY = REPORT_DIR / "phase2x_nishiwaki_ccb_test_summary.tsv"
OUT_DECISION = REPORT_DIR / "phase2x_nishiwaki_final_decision_report.txt"
OUT_STATUS = REPORT_DIR / "phase2x_nishiwaki_ccb_status.txt"

WEIGHTS = {"high": 1.0, "medium": 0.5, "low": 0.25, "unknown": 0.5}


def norm(value: object) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def normalize_taxon(value: object) -> str:
    text = norm(value).replace("_", " ")
    text = re.sub(r"\b[kpcofgs]__", " ", text, flags=re.I)
    text = re.sub(r"\[[^\]]+\]", " ", text)
    text = re.sub(r"\b(strain|subsp|subspecies|serovar|serotype|complex|group)\b", " ", text, flags=re.I)
    text = re.sub(r"[^A-Za-z0-9 ]+", " ", text)
    return re.sub(r"\s+", " ", text).strip().lower()


def binomial(value: object) -> str:
    parts = normalize_taxon(value).split()
    if len(parts) >= 2:
        return " ".join(parts[:2])
    return " ".join(parts)


def detect_country(row: pd.Series) -> str:
    country_cols = [c for c in row.index if re.search(r"country|geo_loc_name", c, re.I)]
    values = [norm(row.get(c, "")) for c in country_cols if norm(row.get(c, ""))]
    return values[0] if values else ""


def select_nishiwaki_rows(df: pd.DataFrame) -> pd.DataFrame:
    if "study_name" not in df.columns:
        raise RuntimeError("ASAP-MAC metadata missing study_name")
    exact = df[df["study_name"].astype(str) == "NishiwakiH_2024"].copy()
    if not exact.empty:
        return exact
    return df[df["study_name"].astype(str).str.contains(r"Nishiwaki|NishiwakiH|NishiwakiH_2024", case=False, na=False)].copy()


def is_fecal(value: object) -> bool:
    return bool(re.search(r"feces|faeces|fecal|stool", norm(value), re.I))


def map_phenotype(row: pd.Series) -> str:
    control = norm(row.get("control", ""))
    disease = norm(row.get("disease", ""))
    target = norm(row.get("target_condition", ""))
    combined = f"{control} | {disease} | {target}"
    if re.search(r"study control|control|healthy", control, re.I):
        return "HC"
    if re.search(r"study case|case", control, re.I):
        return "PD"
    if re.search(r"parkinson'?s disease|parkinson disease|\bPD\b", combined, re.I):
        return "PD"
    return "UNMAPPED"


def build_manifest() -> pd.DataFrame:
    df = pd.read_csv(METADATA, sep="\t", dtype=str)
    nish = select_nishiwaki_rows(df)
    if nish.empty:
        return pd.DataFrame(
            columns=[
                "uuid",
                "study_name",
                "phenotype_CCB",
                "body_site",
                "target_condition",
                "control",
                "disease",
                "country",
                "age",
                "sex",
            ]
        )
    if "body_site" not in nish.columns:
        nish["body_site"] = ""
    nish = nish[nish["body_site"].map(is_fecal)].copy()
    nish["phenotype_CCB"] = nish.apply(map_phenotype, axis=1)
    nish["uuid"] = nish.get("uuid", "").map(norm)
    nish = nish[nish["phenotype_CCB"].isin(["PD", "HC"]) & nish["uuid"].ne("")].copy()
    nish["country"] = nish.apply(detect_country, axis=1)

    for col in ["target_condition", "control", "disease", "age", "sex"]:
        if col not in nish.columns:
            nish[col] = ""

    manifest = nish[
        [
            "uuid",
            "study_name",
            "phenotype_CCB",
            "body_site",
            "target_condition",
            "control",
            "disease",
            "country",
            "age",
            "sex",
        ]
    ].drop_duplicates("uuid")
    return manifest


def write_status(status: str, decision: str, notes: str = "", test: pd.DataFrame | None = None) -> None:
    OUT_STATUS.write_text(status + "\n", encoding="utf-8")
    lines = [
        "PHASE2X_NISHIWAKIH_2024_CCB_ANALYSIS",
        f"status={status}",
        f"decision={decision}",
        "cohort=NishiwakiH_2024",
        f"metadata_file={METADATA.relative_to(ROOT)}",
        f"hf_parquet_paths_file={HF_PARQUET_PATHS.relative_to(ROOT)}",
        f"processed_abundance_resource={RELATIVE_ABUNDANCE_PARQUET}",
        "analysis_type=processed_relative_abundance_only",
        "causality_claimed=no",
        "diagnostic_biomarker_claimed=no",
    ]
    if test is not None and not test.empty:
        for key, value in test.iloc[0].to_dict().items():
            lines.append(f"{key}={value}")
    if notes is not None and not (hasattr(notes, 'empty') and notes.empty):
        lines.append(f"notes={notes}")
    lines.append(
        "interpretation=NishiwakiH_2024 provides an independent processed Japanese fecal metagenomic cohort for CCB testing."
    )
    OUT_DECISION.write_text("\n".join(lines) + "\n", encoding="utf-8")


def query_processed_abundance(manifest: pd.DataFrame) -> pd.DataFrame:
    try:
        import duckdb  # type: ignore
    except Exception as exc:
        raise RuntimeError("duckdb_not_available") from exc

    con = duckdb.connect(database=":memory:")
    try:
        con.execute("LOAD httpfs")
    except Exception as exc:
        raise RuntimeError("duckdb_httpfs_not_available") from exc
    uuid_df = manifest[["uuid"]].drop_duplicates()
    con.register("uuid_filter", uuid_df)
    query = f"""
        SELECT *
        FROM read_parquet('{RELATIVE_ABUNDANCE_PARQUET}')
        WHERE uuid IN (SELECT uuid FROM uuid_filter)
    """
    return con.execute(query).fetchdf()


def species_label(row: pd.Series) -> str:
    for column in ["clade_name_species", "clade_name_terminal", "clade_name"]:
        value = norm(row.get(column, ""))
        if value:
            label = value.split("|")[-1]
            label = re.sub(r"^[a-z]__", "", label)
            return re.sub(r"\s+", " ", label.replace("_", " ")).strip()
    return ""


def species_like(abundance: pd.DataFrame) -> pd.DataFrame:
    if "clade_name_species" in abundance.columns:
        mask = abundance["clade_name_species"].map(norm).ne("")
    else:
        mask = abundance["clade_name"].astype(str).str.contains(r"\|s__|^s__", regex=True, na=False)
    out = abundance[mask].copy()
    out["nishiwaki_taxon"] = out.apply(species_label, axis=1)
    out["taxon_norm"] = out["nishiwaki_taxon"].map(normalize_taxon)
    out["taxon_binomial"] = out["nishiwaki_taxon"].map(binomial)
    return out


def candidate_file() -> Path:
    if CURATED_CANDIDATES.exists():
        return CURATED_CANDIDATES
    return FALLBACK_CANDIDATES


def load_candidates() -> pd.DataFrame:
    path = candidate_file()
    if not path.exists():
        raise RuntimeError("curli_candidate_file_missing")
    df = pd.read_csv(path, sep="\t", dtype=str)
    taxon_col = next((c for c in df.columns if re.search(r"taxon|species|organism|name", c, re.I)), None)
    if taxon_col is None:
        raise RuntimeError("curli_candidate_taxon_column_missing")
    conf_col = next((c for c in df.columns if re.search(r"confidence|evidence|tier", c, re.I)), None)
    weight_col = next((c for c in df.columns if re.search(r"weight", c, re.I)), None)
    out = pd.DataFrame({"candidate_taxon": df[taxon_col].map(norm)})
    out["candidate_norm"] = out["candidate_taxon"].map(normalize_taxon)
    out["candidate_binomial"] = out["candidate_taxon"].map(binomial)
    if conf_col:
        out["confidence"] = df[conf_col].astype(str).str.lower()
    else:
        out["confidence"] = "unknown"
    out["confidence"] = out["confidence"].map(lambda x: x if x in WEIGHTS else "unknown")
    if weight_col:
        out["weight"] = pd.to_numeric(df[weight_col], errors="coerce").fillna(out["confidence"].map(WEIGHTS))
    else:
        out["weight"] = out["confidence"].map(WEIGHTS)
    out = out[out["candidate_binomial"].str.contains(" ", regex=False, na=False)].copy()
    return out.drop_duplicates(["candidate_binomial", "candidate_taxon"])


def match_curli_taxa(species: pd.DataFrame, candidates: pd.DataFrame) -> pd.DataFrame:
    candidate_records = candidates.to_dict("records")
    records: dict[str, dict[str, Any]] = {}
    for clade_name, group in species.groupby("clade_name"):
        first = group.iloc[0]
        taxon_norm = first["taxon_norm"]
        taxon_binomial = first["taxon_binomial"]
        best: dict[str, Any] | None = None
        best_level = ""
        for cand in candidate_records:
            if taxon_norm == cand["candidate_norm"]:
                level = "exact_species_match"
            elif taxon_binomial and taxon_binomial == cand["candidate_binomial"]:
                level = "binomial_species_match"
            elif cand["candidate_binomial"] and (
                cand["candidate_binomial"] in taxon_norm or taxon_binomial in cand["candidate_norm"]
            ):
                level = "species_text_match"
            else:
                continue
            if best is None or float(cand["weight"]) > float(best["weight"]):
                best = cand
                best_level = level
        if best is not None:
            rel = pd.to_numeric(group["relative_abundance"], errors="coerce").fillna(0.0)
            records[str(clade_name)] = {
                "nishiwaki_taxon": first["nishiwaki_taxon"],
                "matched_candidate_taxon": best["candidate_taxon"],
                "match_level": best_level,
                "confidence": best["confidence"],
                "weight": float(best["weight"]),
                "n_samples_present": int((rel > 0).sum()),
                "mean_relative_abundance": float(rel.mean()),
                "max_relative_abundance": float(rel.max()),
                "clade_name": clade_name,
            }
    matched = pd.DataFrame(records.values())
    if matched.empty:
        return pd.DataFrame(
            columns=[
                "nishiwaki_taxon",
                "matched_candidate_taxon",
                "match_level",
                "confidence",
                "weight",
                "n_samples_present",
                "mean_relative_abundance",
                "max_relative_abundance",
                "clade_name",
            ]
        )
    return matched.sort_values(["match_level", "nishiwaki_taxon"])


def compute_ccb(manifest: pd.DataFrame, species: pd.DataFrame, matched: pd.DataFrame) -> pd.DataFrame:
    # CCB_s = sum_i abundance_is * weight_i
    weights = matched.set_index("clade_name")["weight"].astype(float).to_dict()
    subset = species[species["clade_name"].isin(weights)].copy()
    subset["relative_abundance"] = pd.to_numeric(subset["relative_abundance"], errors="coerce").fillna(0.0)
    subset["weighted_abundance"] = subset["relative_abundance"] * subset["clade_name"].map(weights).astype(float)
    ccb = subset.groupby("uuid")["weighted_abundance"].sum().rename("CCB")
    present = subset[subset["relative_abundance"] > 0].groupby("uuid")["clade_name"].nunique().rename("n_matched_curli_taxa_present")
    out = manifest[["uuid", "phenotype_CCB"]].copy()
    out = out.merge(ccb, on="uuid", how="left").merge(present, on="uuid", how="left")
    out["CCB"] = out["CCB"].fillna(0.0).astype(float)
    out["log1p_CCB"] = out["CCB"].map(math.log1p)
    out["n_matched_curli_taxa_present"] = out["n_matched_curli_taxa_present"].fillna(0).astype(int)
    return out[["uuid", "phenotype_CCB", "CCB", "log1p_CCB", "n_matched_curli_taxa_present"]]


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


def make_test_summary(sample_ccb: pd.DataFrame, matched: pd.DataFrame) -> pd.DataFrame:
    pd_vals = sample_ccb.loc[sample_ccb["phenotype_CCB"] == "PD", "CCB"]
    hc_vals = sample_ccb.loc[sample_ccb["phenotype_CCB"] == "HC", "CCB"]
    u_stat, p_value = mann_whitney(pd_vals, hc_vals)
    return pd.DataFrame(
        [
            {
                "n_samples": int(len(sample_ccb)),
                "n_PD": int((sample_ccb["phenotype_CCB"] == "PD").sum()),
                "n_HC": int((sample_ccb["phenotype_CCB"] == "HC").sum()),
                "n_matched_curli_taxa": int(len(matched)),
                "PD_mean_CCB": float(pd_vals.mean()) if len(pd_vals) else float("nan"),
                "HC_mean_CCB": float(hc_vals.mean()) if len(hc_vals) else float("nan"),
                "PD_median_CCB": float(pd_vals.median()) if len(pd_vals) else float("nan"),
                "HC_median_CCB": float(hc_vals.median()) if len(hc_vals) else float("nan"),
                "mannwhitney_U": u_stat,
                "mannwhitney_p": p_value,
                "cliffs_delta": cliffs_delta(pd_vals, hc_vals),
                "effect_direction_mean": "PD>HC" if pd_vals.mean() > hc_vals.mean() else "PD<=HC",
                "effect_direction_median": "PD>HC" if pd_vals.median() > hc_vals.median() else "PD<=HC",
            }
        ]
    )


def decide(summary: pd.Series) -> str:
    if int(summary["n_PD"]) < 2 or int(summary["n_HC"]) < 2 or int(summary["n_matched_curli_taxa"]) == 0:
        return "DATA_NOT_READY"
    mean_support = summary["effect_direction_mean"] == "PD>HC"
    median_support = summary["effect_direction_median"] == "PD>HC"
    p_value = pd.to_numeric(pd.Series([summary["mannwhitney_p"]]), errors="coerce").iloc[0]
    if mean_support and median_support and pd.notna(p_value) and float(p_value) < 0.05:
        return "STRONG_SUPPORT"
    if mean_support and median_support:
        return "DIRECTIONALLY_CONSISTENT_NOT_SIGNIFICANT"
    return "MIXED_OR_NO_SUPPORT"


def main() -> None:
    manifest = build_manifest()
    manifest.to_csv(OUT_MANIFEST, sep="\t", index=False)
    if manifest.empty or (manifest["phenotype_CCB"] == "PD").sum() < 2 or (manifest["phenotype_CCB"] == "HC").sum() < 2:
        write_status("DATA_NOT_READY", "DATA_NOT_READY", notes="insufficient NishiwakiH_2024 mapped PD/HC fecal UUIDs")
        print("STATUS=DATA_NOT_READY")
        return

    try:
        abundance = query_processed_abundance(manifest)
    except RuntimeError as exc:
        write_status("DATA_NOT_READY", "DATA_NOT_READY", notes=str(exc))
        print("STATUS=DATA_NOT_READY")
        return

    abundance.to_csv(OUT_LONG, sep="\t", index=False)
    species = species_like(abundance)
    species.to_csv(OUT_SPECIES_LIKE, sep="\t", index=False)

    candidates = load_candidates()
    matched = match_curli_taxa(species, candidates)
    matched_public = matched.drop(columns=["clade_name"], errors="ignore")
    matched_public.to_csv(OUT_MATCHED_TAXA, sep="\t", index=False)

    if matched.empty:
        pd.DataFrame().to_csv(OUT_SAMPLE_CCB, sep="\t", index=False)
        pd.DataFrame().to_csv(OUT_TEST_SUMMARY, sep="\t", index=False)
        write_status("DATA_NOT_READY", "DATA_NOT_READY", notes="no curated curli-carrier taxa matched NishiwakiH_2024 species-like abundance")
        print("STATUS=DATA_NOT_READY")
        return

    sample_ccb = compute_ccb(manifest, species, matched)
    sample_ccb.to_csv(OUT_SAMPLE_CCB, sep="\t", index=False)
    test_summary = make_test_summary(sample_ccb, matched)
    test_summary.to_csv(OUT_TEST_SUMMARY, sep="\t", index=False)
    decision = decide(test_summary.iloc[0])
    status = "PASS" if decision != "DATA_NOT_READY" else "DATA_NOT_READY"
    write_status(status, decision, notes=None)
    print(f"STATUS={status}")
    print(f"decision={decision}")


if __name__ == "__main__":
    main()
