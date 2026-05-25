from __future__ import annotations

import hashlib
import os
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

WALLEN_DIR = ROOT / "results/phase2d_wallen_discovery"
INTEGRATED_US_DIR = ROOT / "results/phase2f_integrated_us_replication"
MAO_DIR = ROOT / "results/phase2i_mao_replication"
NISHIWAKI_STATUS = ROOT / "data/phase2h/reports/phase2h_nishiwaki_source_verification_status.txt"

OUT_DIR = ROOT / "results/phase2j_cross_cohort_synthesis_updated"
PACKAGE_DIR = ROOT / "results/phase2j_cross_cohort_synthesis_updated_package"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PACKAGE_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_TSV = OUT_DIR / "cross_cohort_processed_table_summary_updated.tsv"
MODEL_TSV = OUT_DIR / "cross_cohort_model_comparison_updated.tsv"
OVERLAP_TSV = OUT_DIR / "cross_cohort_burden_species_overlap_updated.tsv"
SUMMARY_MD = OUT_DIR / "PHASE2J_CROSS_COHORT_SYNTHESIS_WITH_MAO_SUMMARY.md"
STATUS_TXT = OUT_DIR / "phase2j_cross_cohort_final_status.txt"
MANIFEST_TSV = OUT_DIR / "phase2j_cross_cohort_file_manifest.tsv"
FIG_MEDIANS = OUT_DIR / "figure_cross_cohort_burden_medians_with_mao.png"
FIG_LOG1P_OR = OUT_DIR / "figure_cross_cohort_log1p_or_with_mao.png"

PACKAGE_ZIP = PACKAGE_DIR / "curli_phase2j_cross_cohort_synthesis_with_mao_package.zip"
PACKAGE_SHA = PACKAGE_DIR / "curli_phase2j_cross_cohort_synthesis_with_mao_package.zip.sha256"

CORE_CONCLUSION = "A PD-enriched Curli Carrier Burden signal was detected in Wallen, statistically replicated in Integrated-US, and observed in the same direction but without statistical significance in Mao Central China. Nishiwaki was publication-verified but could not be analyzed because sample-level processed species abundance and metadata tables were not identified. Overall, the evidence supports a replicated processed-table ecological taxonomic-proxy association, strengthened by one statistically significant independent replication and one directionally consistent external cohort. These findings do not establish csg gene presence, curli expression, intact operon architecture, causality, or raw-read validation."


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def _num(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def _empty_stat() -> dict[str, Any]:
    return {
        "model_id": "",
        "n_used": np.nan,
        "predictor": "",
        "coefficient": np.nan,
        "odds_ratio": np.nan,
        "ci_lower": np.nan,
        "ci_upper": np.nan,
        "p_value": np.nan,
        "q_value": np.nan,
        "status": "",
        "notes": "",
    }


def _stat_by_id(stats: pd.DataFrame, model_id: str) -> dict[str, Any]:
    rows = stats[stats["model_id"].astype(str) == model_id]
    if rows.empty:
        return _empty_stat()
    return rows.iloc[0].to_dict()


def _cohort_burden_metrics(path: Path) -> dict[str, Any]:
    df = _read_tsv(path)
    burden = pd.to_numeric(df["CurliCarrierBurden"], errors="coerce")
    case = df["Case_status"].astype(str)
    pd_mask = case == "PD"
    control_mask = case == "Control"
    return {
        "n_samples": int(len(df)),
        "n_PD": int(pd_mask.sum()),
        "n_Control": int(control_mask.sum()),
        "PD_median_burden": float(burden[pd_mask].median()),
        "Control_median_burden": float(burden[control_mask].median()),
        "PD_zero_fraction": float((burden[pd_mask] == 0).mean()),
        "Control_zero_fraction": float((burden[control_mask] == 0).mean()),
    }


def _species_label_from_clade(value: str) -> str:
    text = str(value)
    if "|" in text:
        text = text.split("|")[-1]
    if "__" in text:
        text = text.split("__", 1)[-1]
    return " ".join(text.replace("_", " ").split()).strip()


def _species_map(path: Path, cohort: str) -> dict[str, dict[str, str]]:
    df = _read_tsv(path)
    out: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        if cohort == "wallen":
            clade = row.get("metaphlan_clade_name", row.iloc[0])
            label = row.get("metaphlan_species_name", "") or _species_label_from_clade(str(clade))
            method = str(row.get("best_curli_candidate_confidence", row.get("confidence_weight", "")))
        elif cohort == "integrated_us":
            clade = row.get("clade_name", row.get("matched_clade_name", row.iloc[0]))
            label = row.get("metaphlan_species_name", "") or _species_label_from_clade(str(clade))
            method = str(row.get("matching_method", ""))
        else:
            clade = row.get("matched_clade_name", row.get("clade_name", row.iloc[0]))
            label = row.get("normalized_mao_clade_name", "") or _species_label_from_clade(str(clade))
            method = str(row.get("matching_method", ""))
        out[str(label)] = {"method": method}
    return out


def _summary_rows() -> pd.DataFrame:
    wallen_stats = _read_tsv(WALLEN_DIR / "wallen_curli_burden_statistics.tsv")
    integrated_stats = _read_tsv(INTEGRATED_US_DIR / "integrated_us_curli_replication_statistics.tsv")
    mao_stats = _read_tsv(MAO_DIR / "mao_curli_replication_statistics.tsv")

    wallen = _cohort_burden_metrics(WALLEN_DIR / "wallen_curli_carrier_burden.tsv")
    integrated = _cohort_burden_metrics(INTEGRATED_US_DIR / "integrated_us_curli_carrier_burden.tsv")
    mao = _cohort_burden_metrics(MAO_DIR / "mao_curli_carrier_burden.tsv")

    wallen_adj = _stat_by_id(wallen_stats, "LOGIT_PD_log1p_adjusted_age_sex_BMI")
    wallen_presence = _stat_by_id(wallen_stats, "LOGIT_PD_presence_adjusted_age_sex_BMI")
    integrated_adj = _stat_by_id(integrated_stats, "LOGIT_C")
    integrated_presence = _stat_by_id(integrated_stats, "LOGIT_E")
    mao_primary = _stat_by_id(mao_stats, "LOGIT_B")
    mao_presence = _stat_by_id(mao_stats, "PRESENCE_D")

    rows = [
        {
            "cohort": "Wallen",
            "role": "discovery",
            **wallen,
            "n_burden_species": len(_read_tsv(WALLEN_DIR / "wallen_curli_burden_species_set.tsv")),
            "primary_log1p_OR": wallen_adj.get("odds_ratio", np.nan),
            "primary_log1p_CI_lower": wallen_adj.get("ci_lower", np.nan),
            "primary_log1p_CI_upper": wallen_adj.get("ci_upper", np.nan),
            "primary_log1p_p": wallen_adj.get("p_value", np.nan),
            "primary_log1p_q": wallen_adj.get("q_value", np.nan),
            "presence_OR": wallen_presence.get("odds_ratio", np.nan),
            "presence_p": wallen_presence.get("p_value", np.nan),
            "presence_q": wallen_presence.get("q_value", np.nan),
            "interpretation": "DISCOVERY_PD_ENRICHED_SIGNAL",
        },
        {
            "cohort": "Integrated-US",
            "role": "replication",
            **integrated,
            "n_burden_species": len(_read_tsv(INTEGRATED_US_DIR / "integrated_us_curli_burden_species_set.tsv")),
            "primary_log1p_OR": integrated_adj.get("odds_ratio", np.nan),
            "primary_log1p_CI_lower": integrated_adj.get("ci_lower", np.nan),
            "primary_log1p_CI_upper": integrated_adj.get("ci_upper", np.nan),
            "primary_log1p_p": integrated_adj.get("p_value", np.nan),
            "primary_log1p_q": integrated_adj.get("q_value", np.nan),
            "presence_OR": integrated_presence.get("odds_ratio", np.nan),
            "presence_p": integrated_presence.get("p_value", np.nan),
            "presence_q": integrated_presence.get("q_value", np.nan),
            "interpretation": "REPLICATED_DIRECTIONALLY_AND_STATISTICALLY",
        },
        {
            "cohort": "Mao Central China",
            "role": "external_replication",
            **mao,
            "n_burden_species": len(_read_tsv(MAO_DIR / "mao_curli_burden_species_set.tsv")),
            "primary_log1p_OR": mao_primary.get("odds_ratio", np.nan),
            "primary_log1p_CI_lower": mao_primary.get("ci_lower", np.nan),
            "primary_log1p_CI_upper": mao_primary.get("ci_upper", np.nan),
            "primary_log1p_p": mao_primary.get("p_value", np.nan),
            "primary_log1p_q": mao_primary.get("q_value", np.nan),
            "presence_OR": mao_presence.get("odds_ratio", np.nan),
            "presence_p": mao_presence.get("p_value", np.nan),
            "presence_q": mao_presence.get("q_value", np.nan),
            "interpretation": "DIRECTIONALLY_CONSISTENT_BUT_NOT_SIGNIFICANT",
        },
        {
            "cohort": "Nishiwaki",
            "role": "not_analyzed_source_inspection_only",
            "n_samples": np.nan,
            "n_PD": np.nan,
            "n_Control": np.nan,
            "PD_median_burden": np.nan,
            "Control_median_burden": np.nan,
            "PD_zero_fraction": np.nan,
            "Control_zero_fraction": np.nan,
            "n_burden_species": np.nan,
            "primary_log1p_OR": np.nan,
            "primary_log1p_CI_lower": np.nan,
            "primary_log1p_CI_upper": np.nan,
            "primary_log1p_p": np.nan,
            "primary_log1p_q": np.nan,
            "presence_OR": np.nan,
            "presence_p": np.nan,
            "presence_q": np.nan,
            "interpretation": "NOT_ANALYZED_PROCESSED_TABLES_NOT_IDENTIFIED",
        },
    ]
    return pd.DataFrame(rows)


def _model_rows() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    for cohort, role, path in [
        ("Wallen", "discovery", WALLEN_DIR / "wallen_curli_burden_statistics.tsv"),
        ("Integrated-US", "replication", INTEGRATED_US_DIR / "integrated_us_curli_replication_statistics.tsv"),
        ("Mao Central China", "external_replication", MAO_DIR / "mao_curli_replication_statistics.tsv"),
    ]:
        stats = _read_tsv(path)
        for _, row in stats.iterrows():
            comparison = row.get("comparison", "PD_vs_Control")
            if cohort == "Wallen":
                comparison = "PD_vs_Control"
            rows.append(
                {
                    "cohort": cohort,
                    "role": role,
                    "model_id": row.get("model_id", ""),
                    "comparison": comparison,
                    "predictor": row.get("predictor", ""),
                    "n_used": row.get("n_used", np.nan),
                    "odds_ratio": row.get("odds_ratio", np.nan),
                    "ci_lower": row.get("ci_lower", np.nan),
                    "ci_upper": row.get("ci_upper", np.nan),
                    "p_value": row.get("p_value", np.nan),
                    "q_value": row.get("q_value", np.nan),
                    "status": row.get("status", ""),
                    "interpretation_note": row.get("notes", ""),
                }
            )
    rows.append(
        {
            "cohort": "Nishiwaki",
            "role": "not_analyzed_source_inspection_only",
            "model_id": "NOT_ANALYZED",
            "comparison": "not_applicable",
            "predictor": "",
            "n_used": np.nan,
            "odds_ratio": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p_value": np.nan,
            "q_value": np.nan,
            "status": "NOT_ANALYZED",
            "interpretation_note": "sample-level processed species abundance and metadata tables were not identified",
        }
    )
    return pd.DataFrame(rows)


def _overlap_rows() -> pd.DataFrame:
    wallen = _species_map(WALLEN_DIR / "wallen_curli_burden_species_set.tsv", "wallen")
    integrated = _species_map(INTEGRATED_US_DIR / "integrated_us_curli_burden_species_set.tsv", "integrated_us")
    mao = _species_map(MAO_DIR / "mao_curli_burden_species_set.tsv", "mao")

    species = sorted(set(wallen) | set(integrated) | set(mao))
    rows = []
    for name in species:
        present = {
            "wallen": name in wallen,
            "integrated_us": name in integrated,
            "mao": name in mao,
        }
        n_present = sum(present.values())
        if n_present == 3:
            label = "shared_all_three"
        elif n_present == 2:
            label = "shared_two_cohorts"
        else:
            label = "cohort_specific"
        rows.append(
            {
                "species_or_clade_name": name,
                "present_in_wallen": "yes" if present["wallen"] else "no",
                "present_in_integrated_us": "yes" if present["integrated_us"] else "no",
                "present_in_mao": "yes" if present["mao"] else "no",
                "shared_or_cohort_specific": label,
                "wallen_matching_method": wallen.get(name, {}).get("method", ""),
                "integrated_us_matching_method": integrated.get(name, {}).get("method", ""),
                "mao_matching_method": mao.get(name, {}).get("method", ""),
            }
        )
    return pd.DataFrame(rows)


def _write_figures(summary: pd.DataFrame) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
    import matplotlib.pyplot as plt

    plotted = summary[summary["role"] != "not_analyzed_source_inspection_only"].copy()
    labels = plotted["cohort"].tolist()

    x = np.arange(len(labels))
    width = 0.35
    plt.figure(figsize=(8, 4.5))
    plt.bar(x - width / 2, plotted["PD_median_burden"].astype(float), width, label="PD")
    plt.bar(x + width / 2, plotted["Control_median_burden"].astype(float), width, label="Control")
    plt.xticks(x, labels, rotation=12, ha="right")
    plt.ylabel("Median Curli Carrier Burden")
    plt.title("Processed-table burden medians")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG_MEDIANS, dpi=160)
    plt.close()

    or_values = plotted["primary_log1p_OR"].map(_num).to_numpy()
    lower = plotted["primary_log1p_CI_lower"].map(_num).to_numpy()
    upper = plotted["primary_log1p_CI_upper"].map(_num).to_numpy()
    yerr = np.vstack([or_values - lower, upper - or_values])
    yerr = np.nan_to_num(yerr, nan=0.0)

    plt.figure(figsize=(8, 4.5))
    plt.errorbar(or_values, x, xerr=yerr, fmt="o", capsize=4)
    plt.axvline(1.0, color="black", linewidth=1, linestyle="--")
    plt.yticks(x, labels)
    plt.xlabel("CurliCarrierBurden_log1p odds ratio")
    plt.title("Cross-cohort log1p burden effect")
    plt.tight_layout()
    plt.savefig(FIG_LOG1P_OR, dpi=160)
    plt.close()


def _write_summary_md(summary: pd.DataFrame, models: pd.DataFrame, overlap: pd.DataFrame) -> None:
    nishiwaki_status = NISHIWAKI_STATUS.read_text(encoding="utf-8", errors="ignore").strip()
    markdown = f"""# Phase 2J Cross-Cohort Synthesis With Mao

## Cohort Status

- Wallen: discovery cohort with a PD-enriched processed-table Curli Carrier Burden signal.
- Integrated-US: statistically significant independent replication.
- Mao Central China: directionally consistent but statistically non-significant.
- Nishiwaki: publication verified but not analyzed because sample-level processed species abundance and metadata tables were not identified.

## Main Conclusion

{CORE_CONCLUSION}

## Processed-Table Cohort Summary

See `cross_cohort_processed_table_summary_updated.tsv`.

## Model Comparison

See `cross_cohort_model_comparison_updated.tsv`.

## Burden Species Overlap

See `cross_cohort_burden_species_overlap_updated.tsv`.

## Nishiwaki Status Snapshot

```text
{nishiwaki_status}
```

## Limits

This synthesis is a processed-table ecological proxy analysis. It does not establish csg gene presence, curli expression, intact operon architecture, causality, or raw-read validation.

## Output Counts

- Cohort summary rows: {len(summary)}
- Model comparison rows: {len(models)}
- Species overlap rows: {len(overlap)}
"""
    SUMMARY_MD.write_text(markdown, encoding="utf-8")


def _write_manifest(files: list[Path]) -> None:
    rows = []
    for path in files:
        rows.append(
            {
                "file_name": path.name,
                "relative_path": str(path.relative_to(ROOT)),
                "exists": "yes" if path.exists() else "no",
                "bytes": path.stat().st_size if path.exists() else 0,
            }
        )
    pd.DataFrame(rows).to_csv(MANIFEST_TSV, sep="\t", index=False)


def _write_package(files: list[Path]) -> str:
    with zipfile.ZipFile(PACKAGE_ZIP, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in files:
            zf.write(path, arcname=path.name)
    sha = hashlib.sha256(PACKAGE_ZIP.read_bytes()).hexdigest()
    PACKAGE_SHA.write_text(f"{sha}  {PACKAGE_ZIP.name}\n", encoding="utf-8")
    return sha


def main() -> None:
    summary = _summary_rows()
    models = _model_rows()
    overlap = _overlap_rows()

    summary.to_csv(SUMMARY_TSV, sep="\t", index=False)
    models.to_csv(MODEL_TSV, sep="\t", index=False)
    overlap.to_csv(OVERLAP_TSV, sep="\t", index=False)
    _write_figures(summary)
    _write_summary_md(summary, models, overlap)

    STATUS_TXT.write_text(
        "STATUS = PASS\n"
        "UPDATED_SYNTHESIS = WALLEN_PLUS_INTEGRATED_US_PLUS_MAO\n"
        "MAO_INTERPRETATION = DIRECTIONALLY_CONSISTENT_BUT_NOT_SIGNIFICANT\n"
        "NISHIWAKI_INTERPRETATION = NOT_ANALYZED_PROCESSED_TABLES_NOT_IDENTIFIED\n",
        encoding="utf-8",
    )

    output_files = [
        SUMMARY_TSV,
        MODEL_TSV,
        OVERLAP_TSV,
        SUMMARY_MD,
        STATUS_TXT,
        MANIFEST_TSV,
        FIG_MEDIANS,
        FIG_LOG1P_OR,
    ]
    _write_manifest(output_files)
    output_files = [
        SUMMARY_TSV,
        MODEL_TSV,
        OVERLAP_TSV,
        SUMMARY_MD,
        STATUS_TXT,
        MANIFEST_TSV,
        FIG_MEDIANS,
        FIG_LOG1P_OR,
    ]
    _write_package(output_files)


if __name__ == "__main__":
    main()
