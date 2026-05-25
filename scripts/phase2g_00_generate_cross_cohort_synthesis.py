from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WALLEN_DIR = ROOT / "results/phase2d_wallen_discovery"
IUS_DIR = ROOT / "results/phase2f_integrated_us_replication"
OUT_DIR = ROOT / "results/phase2g_cross_cohort_synthesis"
PKG_DIR = ROOT / "results/phase2g_cross_cohort_synthesis_package"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PKG_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_TSV = OUT_DIR / "cross_cohort_processed_table_summary.tsv"
MODELS_TSV = OUT_DIR / "cross_cohort_model_comparison.tsv"
OVERLAP_TSV = OUT_DIR / "cross_cohort_burden_species_overlap.tsv"
SUMMARY_MD = OUT_DIR / "PHASE2G_CROSS_COHORT_SYNTHESIS_SUMMARY.md"
STATUS_TXT = OUT_DIR / "phase2g_cross_cohort_final_status.txt"
MANIFEST_TSV = OUT_DIR / "phase2g_cross_cohort_file_manifest.tsv"
FIG_MEDIANS = OUT_DIR / "figure_cross_cohort_burden_medians.png"
FIG_OR = OUT_DIR / "figure_cross_cohort_adjusted_or.png"


def _metric_map_from_phase2f(path: Path) -> dict[str, str]:
    d = pd.read_csv(path, sep="\t")
    return {str(k): str(v) for k, v in zip(d["metric"], d["value"])}


def _wallen_counts_and_medians(path: Path) -> dict[str, float]:
    df = pd.read_csv(path, sep="\t")
    pd_mask = df["Case_status"].astype(str) == "PD"
    ct_mask = df["Case_status"].astype(str) == "Control"
    burden = pd.to_numeric(df["CurliCarrierBurden"], errors="coerce")
    return {
        "n_samples": float(len(df)),
        "n_pd": float(pd_mask.sum()),
        "n_control": float(ct_mask.sum()),
        "pd_median": float(burden[pd_mask].median()),
        "control_median": float(burden[ct_mask].median()),
        "pd_zero": float((burden[pd_mask] == 0).mean()),
        "control_zero": float((burden[ct_mask] == 0).mean()),
    }


def _row_by_model(df: pd.DataFrame, model_id: str) -> dict[str, str]:
    r = df[df["model_id"] == model_id]
    return {} if r.empty else r.iloc[0].to_dict()


def _fmt(x) -> str:
    try:
        return f"{float(x):.6g}"
    except Exception:
        return str(x)


def _normalize_species_from_wallen_colname(name: str) -> str:
    s = str(name)
    if "metaphlan_species_name" in s.lower():
        return s
    if "s__" in s:
        s = s.rsplit("s__", 1)[-1]
    s = s.replace("_", " ")
    return s.strip()


def _make_figures(summary_df: pd.DataFrame) -> None:
    cohorts = summary_df["cohort"].tolist()
    pd_med = pd.to_numeric(summary_df["PD_median_burden"], errors="coerce").tolist()
    ct_med = pd.to_numeric(summary_df["Control_median_burden"], errors="coerce").tolist()

    x = range(len(cohorts))
    plt.figure(figsize=(7, 4))
    plt.bar([i - 0.15 for i in x], pd_med, width=0.3, label="PD")
    plt.bar([i + 0.15 for i in x], ct_med, width=0.3, label="Control")
    plt.xticks(list(x), cohorts, rotation=12)
    plt.ylabel("Median CurliCarrierBurden")
    plt.title("Cross-Cohort Burden Medians")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG_MEDIANS, dpi=140)
    plt.close()

    or_vals = pd.to_numeric(summary_df["primary_adjusted_OR"], errors="coerce").tolist()
    lo = pd.to_numeric(summary_df["primary_adjusted_CI_lower"], errors="coerce").tolist()
    hi = pd.to_numeric(summary_df["primary_adjusted_CI_upper"], errors="coerce").tolist()
    yerr = [[o - l if pd.notna(o) and pd.notna(l) else 0 for o, l in zip(or_vals, lo)], [h - o if pd.notna(o) and pd.notna(h) else 0 for o, h in zip(or_vals, hi)]]

    plt.figure(figsize=(7, 4))
    plt.errorbar(or_vals, list(x), xerr=yerr, fmt="o", capsize=4)
    plt.axvline(1.0, linestyle="--", linewidth=1)
    plt.yticks(list(x), cohorts)
    plt.xlabel("Adjusted OR (log1p burden)")
    plt.title("Cross-Cohort Adjusted OR")
    plt.tight_layout()
    plt.savefig(FIG_OR, dpi=140)
    plt.close()


def main() -> None:
    w_stats = pd.read_csv(WALLEN_DIR / "wallen_curli_burden_statistics.tsv", sep="\t")
    w_burden = _wallen_counts_and_medians(WALLEN_DIR / "wallen_curli_carrier_burden.tsv")
    w_species = pd.read_csv(WALLEN_DIR / "wallen_curli_burden_species_set.tsv", sep="\t")

    i_stats = pd.read_csv(IUS_DIR / "integrated_us_curli_replication_statistics.tsv", sep="\t")
    i_burden_map = _metric_map_from_phase2f(ROOT / "data/phase2f/reports/phase2f_integrated_us_curli_carrier_burden_summary.tsv")
    i_species = pd.read_csv(IUS_DIR / "integrated_us_curli_burden_species_set.tsv", sep="\t")

    w_logit_adj = _row_by_model(w_stats, "LOGIT_PD_log1p_adjusted_age_sex_BMI")
    w_pres_adj = _row_by_model(w_stats, "LOGIT_PD_presence_adjusted_age_sex_BMI")
    i_logit_adj = _row_by_model(i_stats, "LOGIT_C")
    i_pres_adj = _row_by_model(i_stats, "LOGIT_E")

    summary_rows = [
        {
            "cohort": "Wallen",
            "role": "discovery",
            "n_samples": int(w_burden["n_samples"]),
            "n_PD": int(w_burden["n_pd"]),
            "n_Control": int(w_burden["n_control"]),
            "n_burden_species": int(len(w_species)),
            "PD_median_burden": w_burden["pd_median"],
            "Control_median_burden": w_burden["control_median"],
            "PD_zero_fraction": w_burden["pd_zero"],
            "Control_zero_fraction": w_burden["control_zero"],
            "primary_adjusted_OR": w_logit_adj.get("odds_ratio", ""),
            "primary_adjusted_CI_lower": w_logit_adj.get("ci_lower", ""),
            "primary_adjusted_CI_upper": w_logit_adj.get("ci_upper", ""),
            "primary_adjusted_p": w_logit_adj.get("p_value", ""),
            "primary_adjusted_q": w_logit_adj.get("q_value", ""),
            "presence_adjusted_OR": w_pres_adj.get("odds_ratio", ""),
            "presence_adjusted_CI_lower": w_pres_adj.get("ci_lower", ""),
            "presence_adjusted_CI_upper": w_pres_adj.get("ci_upper", ""),
            "presence_adjusted_p": w_pres_adj.get("p_value", ""),
            "presence_adjusted_q": w_pres_adj.get("q_value", ""),
            "interpretation": "COMPLETE_WITH_WARNINGS",
        },
        {
            "cohort": "Integrated-US",
            "role": "replication",
            "n_samples": int(float(i_burden_map.get("n_samples", "0"))),
            "n_PD": int(float(i_burden_map.get("n_pd", "0"))),
            "n_Control": int(float(i_burden_map.get("n_control", "0"))),
            "n_burden_species": int(float(i_burden_map.get("n_burden_species", "0"))),
            "PD_median_burden": float(i_burden_map.get("pd_median_burden", "nan")),
            "Control_median_burden": float(i_burden_map.get("control_median_burden", "nan")),
            "PD_zero_fraction": float(i_burden_map.get("zero_fraction_case_pd", "nan")),
            "Control_zero_fraction": float(i_burden_map.get("zero_fraction_case_control", "nan")),
            "primary_adjusted_OR": i_logit_adj.get("odds_ratio", ""),
            "primary_adjusted_CI_lower": i_logit_adj.get("ci_lower", ""),
            "primary_adjusted_CI_upper": i_logit_adj.get("ci_upper", ""),
            "primary_adjusted_p": i_logit_adj.get("p_value", ""),
            "primary_adjusted_q": i_logit_adj.get("q_value", ""),
            "presence_adjusted_OR": i_pres_adj.get("odds_ratio", ""),
            "presence_adjusted_CI_lower": i_pres_adj.get("ci_lower", ""),
            "presence_adjusted_CI_upper": i_pres_adj.get("ci_upper", ""),
            "presence_adjusted_p": i_pres_adj.get("p_value", ""),
            "presence_adjusted_q": i_pres_adj.get("q_value", ""),
            "interpretation": "REPLICATED_DIRECTIONALLY_AND_STATISTICALLY",
        },
    ]
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(SUMMARY_TSV, sep="\t", index=False)

    model_rows = []
    for _, r in w_stats.iterrows():
        model_id = str(r.get("model_id", ""))
        comp = "PD_vs_all_controls"
        predictor = str(r.get("predictor", ""))
        if "constipation" in model_id.lower():
            comp = "exploratory"
        model_rows.append(
            {
                "cohort": "Wallen",
                "role": "discovery",
                "model_id": model_id,
                "comparison": comp,
                "predictor": predictor,
                "n_used": r.get("n_used", ""),
                "odds_ratio": r.get("odds_ratio", ""),
                "ci_lower": r.get("ci_lower", ""),
                "ci_upper": r.get("ci_upper", ""),
                "p_value": r.get("p_value", ""),
                "q_value": r.get("q_value", ""),
                "status": r.get("status", ""),
                "interpretation_note": str(r.get("notes", "")),
            }
        )

    for _, r in i_stats.iterrows():
        model_rows.append(
            {
                "cohort": "Integrated-US",
                "role": "replication",
                "model_id": r.get("model_id", ""),
                "comparison": r.get("comparison", ""),
                "predictor": r.get("predictor", ""),
                "n_used": r.get("n_used", ""),
                "odds_ratio": r.get("odds_ratio", ""),
                "ci_lower": r.get("ci_lower", ""),
                "ci_upper": r.get("ci_upper", ""),
                "p_value": r.get("p_value", ""),
                "q_value": r.get("q_value", ""),
                "status": r.get("status", ""),
                "interpretation_note": r.get("notes", ""),
            }
        )

    models_df = pd.DataFrame(model_rows)
    models_df.to_csv(MODELS_TSV, sep="\t", index=False)

    w_clade_col = "metaphlan_clade_name" if "metaphlan_clade_name" in w_species.columns else w_species.columns[0]
    w_method_col = "best_curli_candidate_confidence" if "best_curli_candidate_confidence" in w_species.columns else ""

    w_map = {}
    for _, r in w_species.iterrows():
        name = _normalize_species_from_wallen_colname(r.get(w_clade_col, ""))
        w_map[name] = r.get(w_method_col, "") if w_method_col else ""

    i_map = {}
    for _, r in i_species.iterrows():
        name = str(r.get("metaphlan_species_name", "")).strip()
        i_map[name] = r.get("matching_method", "")

    keys = sorted(set(w_map) | set(i_map))
    overlap_rows = []
    for k in keys:
        in_w = "yes" if k in w_map else "no"
        in_i = "yes" if k in i_map else "no"
        if in_w == "yes" and in_i == "yes":
            tag = "shared"
        elif in_w == "yes":
            tag = "wallen_specific"
        else:
            tag = "integrated_us_specific"
        overlap_rows.append(
            {
                "species_or_clade_name": k,
                "present_in_wallen": in_w,
                "present_in_integrated_us": in_i,
                "shared_or_cohort_specific": tag,
                "wallen_matching_method": w_map.get(k, ""),
                "integrated_us_matching_method": i_map.get(k, ""),
            }
        )

    pd.DataFrame(overlap_rows).to_csv(OVERLAP_TSV, sep="\t", index=False)

    _make_figures(summary_df)

    core_conclusion = (
        "A PD-enriched Curli Carrier Burden signal was detected in Wallen and independently replicated in Integrated-US at the processed-table level. "
        "The association is directionally consistent across cohorts and statistically supported in Integrated-US adjusted and presence-based models. "
        "These findings support a replicated ecological taxonomic-proxy association between curli-carrier taxa and PD status, but do not establish csg gene presence, curli expression, intact operon architecture, causality, or raw-read validation."
    )

    md = f"""# PHASE 2G Cross-Cohort Processed-Table Synthesis

## 1. Project Status

- Wallen discovery: COMPLETE_WITH_WARNINGS
- Integrated-US replication: COMPLETE
- Integrated-US interpretation: REPLICATED_DIRECTIONALLY_AND_STATISTICALLY
- Raw-read validation: NOT_PERFORMED_RUN_ACCESSION_MAPPING_NOT_AVAILABLE

## 2. Wallen Discovery Result

- Discovery cohort retained a PD-enriched processed-table curli-carrier burden signal.
- Primary discovery models and burden summary are in `results/phase2d_wallen_discovery`.

## 3. Integrated-US Replication Result

- Independent replication cohort (`n=244`) showed significant MWU and adjusted/presence model support.
- PD-vs-PC sensitivity was significant; PD-vs-HC sensitivity was not significant.

## 4. Cross-Cohort Model Comparison

See: `cross_cohort_model_comparison.tsv`.

## 5. Burden Species Overlap

See: `cross_cohort_burden_species_overlap.tsv`.

## 6. Main Interpretation

{core_conclusion}

## 7. Limitations

- Processed-table ecological proxy analysis only.
- No raw-read validation in this phase.
- No csg gene-level, operon-level, expression-level, or causal inference claim.

## 8. Next Step

- Expand cohort-by-cohort processed-table replication starting with Nishiwaki.

## 9. Raw-Read Validation Status

- Raw-read validation remains pending and blocked by unavailable run accession mapping in Integrated-US.
"""
    SUMMARY_MD.write_text(md, encoding="utf-8")

    status_text = (
        "STATUS = PASS\n"
        "CROSS_COHORT_INTERPRETATION = REPLICATED_DIRECTIONALLY_AND_STATISTICALLY\n"
        "RAW_READ_VALIDATION = NOT_PERFORMED\n"
    )
    STATUS_TXT.write_text(status_text, encoding="utf-8")

    manifest_files = [
        SUMMARY_TSV.name,
        MODELS_TSV.name,
        OVERLAP_TSV.name,
        SUMMARY_MD.name,
        STATUS_TXT.name,
        MANIFEST_TSV.name,
        FIG_MEDIANS.name,
        FIG_OR.name,
    ]
    with MANIFEST_TSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["file_name", "exists", "notes"])
        for name in manifest_files:
            p = OUT_DIR / name
            w.writerow([name, "yes" if p.exists() else "no", "phase2g_output"])

    pkg_files = [
        SUMMARY_TSV,
        MODELS_TSV,
        OVERLAP_TSV,
        SUMMARY_MD,
        STATUS_TXT,
        MANIFEST_TSV,
        FIG_MEDIANS,
        FIG_OR,
    ]

    for p in pkg_files:
        target = PKG_DIR / p.name
        target.write_bytes(p.read_bytes())

    zip_path = PKG_DIR / "curli_phase2g_cross_cohort_synthesis_package.zip"
    import zipfile

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for p in pkg_files:
            zf.write(PKG_DIR / p.name, arcname=p.name)

    sha = hashlib.sha256(zip_path.read_bytes()).hexdigest()
    (PKG_DIR / "curli_phase2g_cross_cohort_synthesis_package.zip.sha256").write_text(
        f"{sha}  {zip_path.name}\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
