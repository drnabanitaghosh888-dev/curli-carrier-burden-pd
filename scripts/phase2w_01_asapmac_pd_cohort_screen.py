#!/usr/bin/env python3
"""Screen ASAP-MAC metadata for additional eligible PD fecal cohorts.

This phase is metadata-only. It identifies candidate cohorts that may be
analyzable later using existing local processed relative-abundance UUID
workflows. It does not calculate Curli Carrier Burden.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
METADATA = ROOT / "data/phase2t/manual_sources/duru_asapmac/asapmac_sampleMetadata_from_source.tsv"
HF_PARQUET_PATHS = ROOT / "data/phase2t/manual_sources/duru_asapmac/hf_parquet_paths.txt"

REPORT_DIR = ROOT / "data/phase2w/reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)

OUT_TABLE = REPORT_DIR / "phase2w_asapmac_eligible_pd_cohorts.tsv"
OUT_SUMMARY = REPORT_DIR / "phase2w_asapmac_pd_cohort_screen_summary.txt"
OUT_STATUS = REPORT_DIR / "phase2w_asapmac_pd_cohort_screen_status.txt"

STATUS = "PASS"

MIN_PD = 20
MIN_HC = 20

ALREADY_ANALYZED_PATTERNS = {
    "DuruIC_2024": re.compile(r"\bDuruIC[_ -]?2024\b|Duru", re.I),
    "Wallen": re.compile(r"Wallen", re.I),
    "Mao": re.compile(r"\bMao|MaoL[_ -]?2021", re.I),
    "Romano": re.compile(r"Romano|Bedarf|Qian|JoS[_ -]?2022", re.I),
    "Integrated-US": re.compile(r"Integrated|Boktor", re.I),
}


def norm(value: object) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def joined_lower(row: pd.Series, columns: list[str]) -> str:
    return " | ".join(norm(row.get(column, "")) for column in columns).lower()


def unique_join(values: pd.Series, limit: int = 25) -> str:
    cleaned = sorted({norm(v) for v in values if norm(v)})
    if not cleaned:
        return "not_available"
    if len(cleaned) > limit:
        return "; ".join(cleaned[:limit]) + f"; ...(+{len(cleaned) - limit} more)"
    return "; ".join(cleaned)


def columns_matching(df: pd.DataFrame, patterns: list[str]) -> list[str]:
    regexes = [re.compile(pattern, re.I) for pattern in patterns]
    return [column for column in df.columns if any(regex.search(column) for regex in regexes)]


def is_fecal_row(row: pd.Series, body_columns: list[str]) -> bool:
    text = joined_lower(row, body_columns)
    return bool(re.search(r"\bfeces\b|\bfaeces\b|\bfecal\b|\bstool\b|human feces|gut metagenome", text))


def map_phenotype(row: pd.Series) -> str:
    control = norm(row.get("control", "")).lower()
    disease = norm(row.get("disease", "")).lower()
    target = norm(row.get("target_condition", "")).lower()
    combined = f"{control} | {disease} | {target}"

    if re.search(r"study control|\bcontrol\b|\bhealthy\b|\bhc\b", control) or re.search(r"\bhealthy\b|\bhc\b", disease):
        return "HC"
    if re.search(r"study case|\bcase\b", control):
        return "PD"
    if re.search(r"parkinson|parkinson'?s|\bpd\b", disease):
        return "PD"
    if re.search(r"parkinson|parkinson'?s|\bpd\b", combined) and not re.search(r"\bcontrol\b|\bhealthy\b", combined):
        return "PD"
    return "UNMAPPED"


def is_pd_related_study(group: pd.DataFrame) -> bool:
    fields = ["study_name", "target_condition", "disease", "control"]
    text = " | ".join(joined_lower(row, fields) for _, row in group.iterrows())
    return bool(re.search(r"parkinson|parkinson'?s|\bpd\b", text)) or bool((group["phenotype"] == "PD").any())


def already_analyzed_label(study_name: str) -> str:
    labels = [label for label, pattern in ALREADY_ANALYZED_PATTERNS.items() if pattern.search(study_name)]
    return ";".join(labels)


def make_row(group: pd.DataFrame, study_name: str, body_cols: list[str], country_cols: list[str]) -> dict[str, object]:
    n_pd = int((group["phenotype"] == "PD").sum())
    n_hc = int((group["phenotype"] == "HC").sum())
    analyzed = already_analyzed_label(study_name)
    has_uuid = bool("uuid" in group.columns and group["uuid"].map(norm).ne("").any())
    already = "yes" if analyzed else "no"
    if already == "yes":
        priority = "already_analyzed"
        action = f"do_not_prioritize_already_analyzed:{analyzed}"
    elif has_uuid:
        priority = "high"
        action = "prioritize_for_existing_processed_relative_abundance_uuid_workflow"
    else:
        priority = "blocked_no_uuid"
        action = "manual_review_needed_no_uuid_mapping"

    target_cols = [c for c in ["target_condition"] if c in group.columns]
    control_cols = [c for c in ["control"] if c in group.columns]
    body_series = [group[c] for c in body_cols if c in group.columns]
    country_series = [group[c] for c in country_cols if c in group.columns]

    return {
        "study_name": study_name,
        "n_total": int(len(group)),
        "n_PD": n_pd,
        "n_HC": n_hc,
        "body_site_values": unique_join(pd.concat(body_series, ignore_index=True)) if body_series else "not_available",
        "target_condition_values": unique_join(pd.concat([group[c] for c in target_cols], ignore_index=True)) if target_cols else "not_available",
        "control_values": unique_join(pd.concat([group[c] for c in control_cols], ignore_index=True)) if control_cols else "not_available",
        "country_values_if_available": unique_join(pd.concat(country_series, ignore_index=True)) if country_series else "not_available",
        "has_uuid": "yes" if has_uuid else "no",
        "already_analyzed": already,
        "priority": priority,
        "recommended_next_action": action,
    }


def nishiwaki_status(study_names: list[str]) -> str:
    hits = [s for s in study_names if re.search(r"Nishiwaki|NishiwakiH[_ -]?2024", s, re.I)]
    if hits:
        return "Nishiwaki present: " + "; ".join(sorted(set(hits)))
    return "Nishiwaki not found in eligible or candidate PD fecal studies."


def chen_taiwan_status(study_names: list[str], country_text: str) -> str:
    hits = [s for s in study_names if re.search(r"ChenS[_ -]?2022|Chen|Taiwan", s, re.I)]
    if re.search(r"Taiwan", country_text, re.I):
        hits.append("Taiwan country value detected")
    if hits:
        return "ChenS_2022/Taiwan-like present: " + "; ".join(sorted(set(hits)))
    return "ChenS_2022/Taiwan-like study not found in eligible or candidate PD fecal studies."


def main() -> None:
    if not METADATA.exists():
        OUT_TABLE.write_text(
            "\t".join(
                [
                    "study_name",
                    "n_total",
                    "n_PD",
                    "n_HC",
                    "body_site_values",
                    "target_condition_values",
                    "control_values",
                    "country_values_if_available",
                    "has_uuid",
                    "already_analyzed",
                    "priority",
                    "recommended_next_action",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        OUT_STATUS.write_text("FAIL\n", encoding="utf-8")
        OUT_SUMMARY.write_text(f"status=FAIL\nmissing_metadata_file={METADATA.relative_to(ROOT)}\n", encoding="utf-8")
        print("STATUS=FAIL")
        return

    df = pd.read_csv(METADATA, sep="\t", dtype=str)
    if "study_name" not in df.columns:
        raise RuntimeError("ASAP-MAC metadata is missing study_name")

    body_cols = columns_matching(df, [r"body_site", r"sample_type", r"environmental_medium", r"env_medium"])
    country_cols = columns_matching(df, [r"country", r"geo_loc_name"])
    if not body_cols:
        body_cols = ["body_site"] if "body_site" in df.columns else []

    df["phenotype"] = df.apply(map_phenotype, axis=1)
    df["is_fecal"] = df.apply(lambda row: is_fecal_row(row, body_cols), axis=1)

    fecal = df[df["is_fecal"] & df["phenotype"].isin(["PD", "HC"])].copy()
    candidate_study_names: list[str] = []
    rows: list[dict[str, object]] = []

    for study_name, group in fecal.groupby("study_name", dropna=False):
        study_name = norm(study_name) or "UNKNOWN_STUDY"
        if not is_pd_related_study(group):
            continue
        candidate_study_names.append(study_name)
        n_pd = int((group["phenotype"] == "PD").sum())
        n_hc = int((group["phenotype"] == "HC").sum())
        if n_pd >= MIN_PD and n_hc >= MIN_HC:
            rows.append(make_row(group, study_name, body_cols, country_cols))

    out = pd.DataFrame(
        rows,
        columns=[
            "study_name",
            "n_total",
            "n_PD",
            "n_HC",
            "body_site_values",
            "target_condition_values",
            "control_values",
            "country_values_if_available",
            "has_uuid",
            "already_analyzed",
            "priority",
            "recommended_next_action",
        ],
    )
    if not out.empty:
        out = out.sort_values(["already_analyzed", "priority", "n_PD", "n_HC", "study_name"], ascending=[True, True, False, False, True])
    out.to_csv(OUT_TABLE, sep="\t", index=False)

    all_studies_for_named_checks = sorted(set(candidate_study_names) | set(out["study_name"].astype(str).tolist() if not out.empty else []))
    country_text = " | ".join(out["country_values_if_available"].astype(str).tolist()) if not out.empty else ""
    n_not_analyzed = int((out["already_analyzed"] == "no").sum()) if not out.empty else 0

    summary = [
        "PHASE 2W ASAP-MAC PD COHORT SCREEN",
        f"status={STATUS}",
        f"metadata_file={METADATA.relative_to(ROOT)}",
        f"hf_parquet_paths_file={HF_PARQUET_PATHS.relative_to(ROOT)}",
        f"n_metadata_rows={len(df)}",
        f"n_fecal_pd_or_hc_rows={len(fecal)}",
        f"n_pd_related_fecal_candidate_studies={len(set(candidate_study_names))}",
        f"n_eligible_studies_min_PD_{MIN_PD}_min_HC_{MIN_HC}={len(out)}",
        f"n_not_yet_analyzed_eligible_studies={n_not_analyzed}",
        f"{nishiwaki_status(all_studies_for_named_checks)}",
        f"{chen_taiwan_status(all_studies_for_named_checks, country_text)}",
        "ccb_computed=no",
        "recommended_use=screening_only_then_manual_selection_for_processed_uuid_abundance_workflow",
    ]
    OUT_SUMMARY.write_text("\n".join(summary) + "\n", encoding="utf-8")
    OUT_STATUS.write_text(STATUS + "\n", encoding="utf-8")
    print(f"STATUS={STATUS}")
    print(f"n_eligible_studies={len(out)}")
    print(f"n_not_yet_analyzed_eligible_studies={n_not_analyzed}")


if __name__ == "__main__":
    main()
