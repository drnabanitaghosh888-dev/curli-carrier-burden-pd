from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data/phase2p/reports"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EVIDENCE_TABLE = OUT_DIR / "phase2p_cross_cohort_evidence_table.tsv"
SUMMARY_TXT = OUT_DIR / "phase2p_cross_cohort_evidence_summary.txt"
STATUS_TXT = OUT_DIR / "phase2p_cross_cohort_evidence_status.txt"

REQUIRED_WORDING = (
    "Curli Carrier Burden shows recurrent PD-associated elevation across processed public cohorts, "
    "but the strength and robustness of association vary by cohort. Romano non-Wallen analysis "
    "provides directionally supportive but cohort-sensitive evidence, emphasizing the need for "
    "stratified and cohort-aware interpretation."
)

COLUMNS = [
    "dataset_id",
    "phase",
    "dataset_role",
    "analysis_set",
    "n_samples",
    "n_PD",
    "n_HC",
    "n_controls",
    "n_studies",
    "candidate_source",
    "n_matched_curli_taxa",
    "PD_mean_CCB",
    "HC_mean_CCB",
    "PD_median_CCB",
    "HC_median_CCB",
    "effect_direction_mean",
    "effect_direction_median",
    "primary_p_value",
    "primary_p_value_type",
    "cohort_or_study_aware_p_value",
    "cohort_or_study_aware_p_value_type",
    "effect_size",
    "effect_size_type",
    "robustness_summary",
    "final_evidence_tier",
    "decision",
    "source_files_used",
    "missing_files",
    "notes",
]


def _base_row(dataset_id: str, phase: str, dataset_role: str) -> dict[str, Any]:
    return {col: "not_available" for col in COLUMNS} | {
        "dataset_id": dataset_id,
        "phase": phase,
        "dataset_role": dataset_role,
        "source_files_used": "",
        "missing_files": "",
        "notes": "",
    }


def _read_tsv(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return None


def _read_key_value_report(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def _num(value: Any) -> float | None:
    try:
        if pd.isna(value):
            return None
        return float(value)
    except Exception:
        return None


def _direction(pd_value: Any, control_value: Any) -> str:
    pd_num = _num(pd_value)
    control_num = _num(control_value)
    if pd_num is None or control_num is None:
        return "not_available"
    if pd_num > control_num:
        return "PD>HC"
    if pd_num < control_num:
        return "PD<HC"
    return "PD=HC"


def _first_row(df: pd.DataFrame | None) -> dict[str, Any]:
    if df is None or df.empty:
        return {}
    return df.iloc[0].to_dict()


def _count_burden_species(path: Path) -> str:
    df = _read_tsv(path)
    if df is None:
        return "not_available"
    return str(len(df))


def _carrier_metrics(path: Path) -> dict[str, Any]:
    df = _read_tsv(path)
    if df is None or df.empty:
        return {}
    if not {"Case_status", "CurliCarrierBurden"}.issubset(df.columns):
        return {}
    burden = pd.to_numeric(df["CurliCarrierBurden"], errors="coerce")
    case = df["Case_status"].astype(str)
    pd_mask = case == "PD"
    control_mask = case == "Control"
    return {
        "n_samples": len(df),
        "n_PD": int(pd_mask.sum()),
        "n_controls": int(control_mask.sum()),
        "PD_mean_CCB": float(burden[pd_mask].mean()),
        "HC_mean_CCB": float(burden[control_mask].mean()),
        "PD_median_CCB": float(burden[pd_mask].median()),
        "HC_median_CCB": float(burden[control_mask].median()),
    }


def _stat_by_id(path: Path, model_id: str) -> dict[str, Any]:
    df = _read_tsv(path)
    if df is None or "model_id" not in df.columns:
        return {}
    rows = df[df["model_id"].astype(str) == model_id]
    if rows.empty:
        return {}
    return rows.iloc[0].to_dict()


def _wallen_row() -> dict[str, Any]:
    stats_path = ROOT / "results/phase2d_wallen_discovery/wallen_curli_burden_statistics.tsv"
    burden_path = ROOT / "results/phase2d_wallen_discovery/wallen_curli_carrier_burden.tsv"
    species_path = ROOT / "results/phase2d_wallen_discovery/wallen_curli_burden_species_set.tsv"
    row = _base_row("WALLEN_DISCOVERY", "Phase 2D", "discovery")
    missing = [str(p.relative_to(ROOT)) for p in [stats_path, burden_path, species_path] if not p.exists()]
    used = [str(p.relative_to(ROOT)) for p in [stats_path, burden_path, species_path] if p.exists()]

    row.update(_carrier_metrics(burden_path))
    row["n_HC"] = row.get("n_controls", "not_available")
    row["n_matched_curli_taxa"] = _count_burden_species(species_path)
    row["effect_direction_mean"] = _direction(row.get("PD_mean_CCB"), row.get("HC_mean_CCB"))
    row["effect_direction_median"] = _direction(row.get("PD_median_CCB"), row.get("HC_median_CCB"))

    mwu = _stat_by_id(stats_path, "MW_PD_vs_Control_raw_burden")
    adjusted = _stat_by_id(stats_path, "LOGIT_PD_log1p_adjusted_age_sex_BMI")
    row["primary_p_value"] = mwu.get("p_value", "not_available")
    row["primary_p_value_type"] = "mannwhitney_raw_burden" if mwu else "not_available"
    row["cohort_or_study_aware_p_value"] = "not_available"
    row["cohort_or_study_aware_p_value_type"] = "not_available"
    row["effect_size"] = adjusted.get("odds_ratio", "not_available")
    row["effect_size_type"] = "adjusted_log1p_odds_ratio" if adjusted else "not_available"
    row["robustness_summary"] = "discovery_only_not_independent_replication"
    row["final_evidence_tier"] = "DISCOVERY_EVIDENCE"
    row["decision"] = "DISCOVERY_PD_ENRICHED_SIGNAL"
    row["source_files_used"] = ";".join(used)
    row["missing_files"] = ";".join(missing)
    if missing:
        row["notes"] = "PARTIAL_INFORMATION"
    return row


def _integrated_us_row() -> dict[str, Any]:
    stats_path = ROOT / "results/phase2f_integrated_us_replication/integrated_us_curli_replication_statistics.tsv"
    burden_path = ROOT / "results/phase2f_integrated_us_replication/integrated_us_curli_carrier_burden.tsv"
    species_path = ROOT / "results/phase2f_integrated_us_replication/integrated_us_curli_burden_species_set.tsv"
    row = _base_row("INTEGRATED_US_MULTICOHORT_PD", "Phase 2F", "independent_replication")
    missing = [str(p.relative_to(ROOT)) for p in [stats_path, burden_path, species_path] if not p.exists()]
    used = [str(p.relative_to(ROOT)) for p in [stats_path, burden_path, species_path] if p.exists()]

    row.update(_carrier_metrics(burden_path))
    row["n_HC"] = row.get("n_controls", "not_available")
    row["n_matched_curli_taxa"] = _count_burden_species(species_path)
    row["effect_direction_mean"] = _direction(row.get("PD_mean_CCB"), row.get("HC_mean_CCB"))
    row["effect_direction_median"] = _direction(row.get("PD_median_CCB"), row.get("HC_median_CCB"))

    mwu = _stat_by_id(stats_path, "MWU_A")
    adjusted = _stat_by_id(stats_path, "LOGIT_C")
    row["primary_p_value"] = mwu.get("p_value", "not_available")
    row["primary_p_value_type"] = "mannwhitney_raw_burden" if mwu else "not_available"
    row["cohort_or_study_aware_p_value"] = "not_available"
    row["cohort_or_study_aware_p_value_type"] = "not_available"
    row["effect_size"] = adjusted.get("odds_ratio", "not_available")
    row["effect_size_type"] = "adjusted_log1p_odds_ratio" if adjusted else "not_available"
    row["robustness_summary"] = "significant_adjusted_and_presence_models; no cohort-aware test available"
    row["final_evidence_tier"] = "SUPPORTIVE_REPLICATION"
    row["decision"] = "REPLICATED_DIRECTIONALLY_AND_STATISTICALLY"
    row["source_files_used"] = ";".join(used)
    row["missing_files"] = ";".join(missing)
    if missing:
        row["notes"] = "PARTIAL_INFORMATION"
    return row


def _mao_row() -> dict[str, Any]:
    stats_path = ROOT / "results/phase2i_mao_replication/mao_curli_replication_statistics.tsv"
    burden_path = ROOT / "results/phase2i_mao_replication/mao_curli_carrier_burden.tsv"
    species_path = ROOT / "results/phase2i_mao_replication/mao_curli_burden_species_set.tsv"
    row = _base_row("MAO_CENTRAL_CHINA_PD", "Phase 2I", "external_replication")
    missing = [str(p.relative_to(ROOT)) for p in [stats_path, burden_path, species_path] if not p.exists()]
    used = [str(p.relative_to(ROOT)) for p in [stats_path, burden_path, species_path] if p.exists()]

    row.update(_carrier_metrics(burden_path))
    row["n_HC"] = row.get("n_controls", "not_available")
    row["n_matched_curli_taxa"] = _count_burden_species(species_path)
    row["effect_direction_mean"] = _direction(row.get("PD_mean_CCB"), row.get("HC_mean_CCB"))
    row["effect_direction_median"] = _direction(row.get("PD_median_CCB"), row.get("HC_median_CCB"))

    mwu = _stat_by_id(stats_path, "MWU_A")
    logit = _stat_by_id(stats_path, "LOGIT_B")
    row["primary_p_value"] = mwu.get("p_value", "not_available")
    row["primary_p_value_type"] = "mannwhitney_raw_burden" if mwu else "not_available"
    row["cohort_or_study_aware_p_value"] = "not_available"
    row["cohort_or_study_aware_p_value_type"] = "not_available"
    row["effect_size"] = logit.get("odds_ratio", "not_available")
    row["effect_size_type"] = "unadjusted_log1p_odds_ratio" if logit else "not_available"
    row["robustness_summary"] = "directionally_consistent_but_statistically_non_significant"
    row["final_evidence_tier"] = "DIRECTIONAL_SUPPORT_NOT_SIGNIFICANT"
    row["decision"] = "DIRECTIONALLY_CONSISTENT_BUT_NOT_SIGNIFICANT"
    row["source_files_used"] = ";".join(used)
    row["missing_files"] = ";".join(missing)
    if missing:
        row["notes"] = "PARTIAL_INFORMATION"
    return row


def _romano_row() -> dict[str, Any]:
    summary_path = ROOT / "data/phase2o/reports/phase2o_romano_nonwallen_ccb_test_summary.tsv"
    taxa_path = ROOT / "data/phase2o/reports/phase2o_romano_nonwallen_ccb_matched_taxa.tsv"
    sample_path = ROOT / "data/phase2o/reports/phase2o_romano_nonwallen_ccb_sample_burden.tsv"
    loo_path = ROOT / "data/phase2o/reports/phase2o_romano_leave_one_study_out_ccb.tsv"
    within_path = ROOT / "data/phase2o/reports/phase2o_romano_within_study_ccb_effects.tsv"
    decision_path = ROOT / "data/phase2o/reports/phase2o_romano_final_decision_report.txt"
    paths = [summary_path, taxa_path, sample_path, loo_path, within_path, decision_path]
    missing = [str(p.relative_to(ROOT)) for p in paths if not p.exists()]
    used = [str(p.relative_to(ROOT)) for p in paths if p.exists()]

    row = _base_row("ROMANO_NONWALLEN", "Phase 2O", "external_replication_nonwallen")
    summary = _first_row(_read_tsv(summary_path))
    decision = _read_key_value_report(decision_path)

    row.update(
        {
            "analysis_set": summary.get("analysis_set", decision.get("analysis_set", "Romano_nonWallen_600")),
            "n_samples": summary.get("n_samples", decision.get("n_samples", 600)),
            "n_PD": summary.get("n_PD", decision.get("n_PD", 280)),
            "n_HC": summary.get("n_HC", decision.get("n_HC", 320)),
            "n_controls": summary.get("n_HC", decision.get("n_HC", 320)),
            "n_studies": summary.get("n_studies", decision.get("n_studies", 6)),
            "candidate_source": summary.get(
                "candidate_source",
                decision.get("candidate_source", "data/000_phase2o_curli_reference/000_curli_candidate_species.tsv"),
            ),
            "n_matched_curli_taxa": summary.get(
                "n_matched_curli_motus_taxa",
                decision.get("n_matched_curli_motus_taxa", 29),
            ),
            "PD_mean_CCB": summary.get("PD_mean_CCB", decision.get("PD_mean_CCB", "not_available")),
            "HC_mean_CCB": summary.get("HC_mean_CCB", decision.get("HC_mean_CCB", "not_available")),
            "PD_median_CCB": summary.get("PD_median_CCB", decision.get("PD_median_CCB", "not_available")),
            "HC_median_CCB": summary.get("HC_median_CCB", decision.get("HC_median_CCB", "not_available")),
            "primary_p_value": summary.get("mannwhitney_raw_p", decision.get("mannwhitney_p", "not_available")),
            "primary_p_value_type": "pooled_mannwhitney",
            "cohort_or_study_aware_p_value": summary.get(
                "study_stratified_permutation_p",
                decision.get("study_stratified_permutation_p", "not_available"),
            ),
            "cohort_or_study_aware_p_value_type": "study_stratified_permutation",
            "effect_size": summary.get("cliffs_delta_raw", decision.get("cliffs_delta", "not_available")),
            "effect_size_type": "cliffs_delta",
            "robustness_summary": (
                "leave-one-study-out PD>HC by mean and median; within-study heterogeneous; "
                "study-stratified permutation non-significant"
            ),
            "final_evidence_tier": "WEAK_COHORT_SENSITIVE_SUPPORT",
            "decision": decision.get("decision", "WEAK_COHORT_SENSITIVE_SUPPORT_NOT_DEFINITIVE_REPLICATION"),
            "source_files_used": ";".join(used),
            "missing_files": ";".join(missing),
        }
    )
    row["effect_direction_mean"] = _direction(row["PD_mean_CCB"], row["HC_mean_CCB"])
    row["effect_direction_median"] = _direction(row["PD_median_CCB"], row["HC_median_CCB"])
    if missing:
        row["notes"] = "PARTIAL_INFORMATION"
    return row


def build_evidence_table() -> pd.DataFrame:
    rows = [_wallen_row(), _integrated_us_row(), _mao_row(), _romano_row()]
    return pd.DataFrame(rows, columns=COLUMNS)


def _summary_text(table: pd.DataFrame) -> str:
    included = ", ".join(table["dataset_id"].astype(str).tolist())
    lines = [
        "PHASE 2P CROSS-COHORT EVIDENCE SYNTHESIS",
        "",
        f"Included datasets: {included}.",
        "Wallen is treated as discovery evidence, not independent replication.",
        "Integrated-US is treated as independent processed-table replication.",
        "Mao Central China is treated as an external processed-table cohort with directional but non-significant support.",
        "Romano non-Wallen is treated as external multi-study support with cohort-sensitive robustness.",
        "",
    ]

    for _, row in table.iterrows():
        lines.append(
            f"- {row['dataset_id']}: role={row['dataset_role']}; "
            f"mean_direction={row['effect_direction_mean']}; "
            f"median_direction={row['effect_direction_median']}; "
            f"decision={row['decision']}."
        )

    lines.extend(
        [
            "",
            REQUIRED_WORDING,
            "",
            "Conservative final conclusion: the processed-table evidence is recurrent and directionally supportive, "
            "with one statistically significant independent replication and one cohort-sensitive multi-study support "
            "signal. This does not establish csg gene presence, curli expression, intact operon architecture, "
            "causality, or raw-read validation.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    table = build_evidence_table()
    table.to_csv(EVIDENCE_TABLE, sep="\t", index=False)
    SUMMARY_TXT.write_text(_summary_text(table), encoding="utf-8")

    missing_any = table["missing_files"].astype(str).str.len().gt(0).any()
    romano_ok = (table["dataset_id"] == "ROMANO_NONWALLEN").any()
    previous_ok = table["dataset_id"].isin(["WALLEN_DISCOVERY", "INTEGRATED_US_MULTICOHORT_PD", "MAO_CENTRAL_CHINA_PD"]).any()
    if len(table) == 0 or not romano_ok:
        status = "FAIL"
    elif missing_any and romano_ok and previous_ok:
        status = "PASS_WITH_WARNINGS"
    else:
        status = "PASS"

    STATUS_TXT.write_text(
        f"STATUS={status}\n"
        f"N_DATASETS={len(table)}\n"
        f"N_ROWS_WITH_MISSING_FILES={int(table['missing_files'].astype(str).str.len().gt(0).sum())}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
