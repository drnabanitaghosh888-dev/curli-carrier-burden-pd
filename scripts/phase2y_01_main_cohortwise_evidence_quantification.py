#!/usr/bin/env python3
from pathlib import Path
import math
import re
import pandas as pd

ROOT = Path(".")
OUTDIR = ROOT / "data/phase2y/reports"
OUTDIR.mkdir(parents=True, exist_ok=True)

MAIN_COHORT_TERMS = [
    "Wallen",
    "Integrated",
    "INTEGRATED",
    "Mao",
    "MAO",
    "Romano",
    "Duru",
    "DuruIC_2024",
]

EXCLUDE_TERMS = [
    "Nishiwaki",
    "Sampson",
    "Stagaman",
    "FoxDEN",
]

INPUTS = [
    ROOT / "data/phase2p/reports/phase2p_cross_cohort_evidence_table.tsv",
    ROOT / "data/phase2u/reports/phase2u_duru_ccb_test_summary.tsv",
]

def read_tsv_if_exists(path):
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()

def safe_float(x):
    try:
        if pd.isna(x):
            return math.nan
        return float(x)
    except Exception:
        return math.nan

def safe_int(x):
    try:
        if pd.isna(x):
            return 0
        return int(float(x))
    except Exception:
        return 0

def first_existing(row, names, default=""):
    for n in names:
        if n in row.index and pd.notna(row[n]):
            return row[n]
    return default

def infer_cohort_name(row):
    candidates = [
        "cohort", "dataset", "dataset_id", "analysis_set", "study_name",
        "source", "analysis_name", "candidate_source"
    ]
    vals = []
    for c in candidates:
        if c in row.index and pd.notna(row[c]):
            vals.append(str(row[c]))
    joined = " ".join(vals)

    if re.search("Duru|DuruIC_2024", joined, re.I):
        return "DuruIC_2024"
    if re.search("Romano", joined, re.I):
        return "Romano_nonWallen"
    if re.search("Mao|MAO", joined, re.I):
        return "Mao_2021_Central_China"
    if re.search("Integrated|INTEGRATED", joined, re.I):
        return "Integrated_US_multicohort"
    if re.search("Wallen", joined, re.I):
        return "Wallen_2022"
    return vals[0] if vals else "UNKNOWN_COHORT"

def include_main_cohort(name):
    if any(re.search(t, name, re.I) for t in EXCLUDE_TERMS):
        return False
    return any(re.search(t, name, re.I) for t in MAIN_COHORT_TERMS)

def evidence_components(row):
    pd_mean = safe_float(row.get("PD_mean_CCB", math.nan))
    hc_mean = safe_float(row.get("HC_mean_CCB", math.nan))
    pd_median = safe_float(row.get("PD_median_CCB", math.nan))
    hc_median = safe_float(row.get("HC_median_CCB", math.nan))

    p = safe_float(first_existing(row, [
        "mannwhitney_p",
        "mannwhitney_raw_p",
        "mannwhitney_log1p_p",
        "p_value",
        "p",
    ], math.nan))

    delta = safe_float(first_existing(row, [
        "cliffs_delta",
        "cliffs_delta_raw",
        "cliffs_delta_log1p",
        "cliff_delta",
    ], math.nan))

    n = safe_int(first_existing(row, ["n_samples", "n_total", "N", "total_samples"], 0))
    n_pd = safe_int(first_existing(row, ["n_PD", "PD_n", "n_cases"], 0))
    n_hc = safe_int(first_existing(row, ["n_HC", "HC_n", "n_controls"], 0))
    n_taxa = safe_int(first_existing(row, [
        "n_matched_curli_taxa",
        "n_matched_curli_motus_taxa",
        "matched_curli_taxa",
        "n_curli_taxa",
    ], 0))

    mean_dir = "PD>HC" if pd_mean > hc_mean else "HC>=PD"
    median_dir = "PD>HC" if pd_median > hc_median else "HC>=PD"

    direction_points = 0
    if mean_dir == "PD>HC":
        direction_points += 2
    if median_dir == "PD>HC":
        direction_points += 2

    if not math.isnan(p):
        if p < 0.01:
            p_points = 2
            p_class = "strong_nominal"
        elif p < 0.05:
            p_points = 1.5
            p_class = "nominal"
        elif p < 0.10:
            p_points = 0.75
            p_class = "trend"
        else:
            p_points = 0
            p_class = "not_significant"
    else:
        p_points = 0
        p_class = "not_available"

    if not math.isnan(delta):
        ad = abs(delta)
        if ad >= 0.33:
            effect_points = 2
            effect_class = "moderate_or_larger"
        elif ad >= 0.147:
            effect_points = 1
            effect_class = "small"
        elif ad > 0:
            effect_points = 0.5
            effect_class = "very_small"
        else:
            effect_points = 0
            effect_class = "zero"
    else:
        effect_points = 0
        effect_class = "not_available"

    if n_taxa >= 10:
        compatibility_points = 1
        compatibility_class = "high"
    elif n_taxa >= 5:
        compatibility_points = 0.5
        compatibility_class = "moderate"
    elif n_taxa > 0:
        compatibility_points = 0.25
        compatibility_class = "low"
    else:
        compatibility_points = 0
        compatibility_class = "unknown_or_zero"

    if n >= 200:
        sample_points = 1
        sample_class = "large"
    elif n >= 100:
        sample_points = 0.75
        sample_class = "moderate"
    elif n >= 40:
        sample_points = 0.5
        sample_class = "small"
    elif n > 0:
        sample_points = 0.25
        sample_class = "very_small"
    else:
        sample_points = 0
        sample_class = "unknown"

    descriptive_score_0_10 = direction_points + p_points + effect_points + compatibility_points + sample_points

    if direction_points == 4 and not math.isnan(p) and p < 0.05 and not math.isnan(delta) and delta > 0:
        evidence_class = "positive_nominal_support"
    elif direction_points == 4 and (math.isnan(p) or p >= 0.05):
        evidence_class = "directionally_consistent_support"
    elif direction_points in [2]:
        evidence_class = "mixed_direction"
    elif direction_points == 0:
        evidence_class = "no_positive_direction"
    else:
        evidence_class = "unclassified"

    return {
        "n_samples": n,
        "n_PD": n_pd,
        "n_HC": n_hc,
        "n_matched_curli_taxa": n_taxa,
        "PD_mean_CCB": pd_mean,
        "HC_mean_CCB": hc_mean,
        "PD_median_CCB": pd_median,
        "HC_median_CCB": hc_median,
        "effect_direction_mean": mean_dir,
        "effect_direction_median": median_dir,
        "mannwhitney_p": p,
        "p_class": p_class,
        "cliffs_delta": delta,
        "effect_class": effect_class,
        "compatibility_class": compatibility_class,
        "sample_size_class": sample_class,
        "direction_points_0_4": direction_points,
        "p_points_0_2": p_points,
        "effect_points_0_2": effect_points,
        "compatibility_points_0_1": compatibility_points,
        "sample_points_0_1": sample_points,
        "descriptive_evidence_score_0_10": descriptive_score_0_10,
        "cohort_evidence_class": evidence_class,
    }

def collect_rows():
    rows = []

    for p in INPUTS:
        df = read_tsv_if_exists(p)
        if df.empty:
            continue

        for _, row in df.iterrows():
            cohort = infer_cohort_name(row)
            if not include_main_cohort(cohort):
                continue

            comp = evidence_components(row)
            comp["cohort"] = cohort
            comp["source_file"] = str(p)
            rows.append(comp)

    # Remove duplicate Duru rows if accidentally collected twice
    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out = out.drop_duplicates(subset=["cohort", "source_file"], keep="first")
    preferred_order = {
        "Wallen_2022": 1,
        "Integrated_US_multicohort": 2,
        "Mao_2021_Central_China": 3,
        "Romano_nonWallen": 4,
        "DuruIC_2024": 5,
    }
    out["display_order"] = out["cohort"].map(preferred_order).fillna(99)
    out = out.sort_values(["display_order", "cohort"]).drop(columns=["display_order"])
    return out

def main():
    table = collect_rows()

    table_path = OUTDIR / "phase2y_main_cohortwise_evidence_table.tsv"
    summary_path = OUTDIR / "phase2y_main_cohortwise_evidence_summary.txt"
    status_path = OUTDIR / "phase2y_main_cohortwise_evidence_status.txt"

    if table.empty:
        status = "NO_MAIN_COHORT_ROWS_FOUND"
        table.to_csv(table_path, sep="\t", index=False)
        summary_path.write_text(
            "PHASE2Y_MAIN_COHORTWISE_EVIDENCE_QUANTIFICATION\n"
            f"status={status}\n"
            "No main-manuscript cohort rows could be extracted from the available Phase 2P/2U reports.\n"
        )
        status_path.write_text(status + "\n")
        print(status)
        return

    table.to_csv(table_path, sep="\t", index=False)

    n_cohorts = len(table)
    n_positive_direction_both = int(((table["effect_direction_mean"] == "PD>HC") & (table["effect_direction_median"] == "PD>HC")).sum())
    n_nominal_positive = int((table["cohort_evidence_class"] == "positive_nominal_support").sum())
    n_directional = int((table["cohort_evidence_class"] == "directionally_consistent_support").sum())
    n_high_or_moderate_compat = int(table["compatibility_class"].isin(["high", "moderate"]).sum())

    score_mean = table["descriptive_evidence_score_0_10"].mean()
    score_median = table["descriptive_evidence_score_0_10"].median()

    summary = [
        "PHASE2Y_MAIN_COHORTWISE_EVIDENCE_QUANTIFICATION",
        "status=PASS",
        "analysis_type=cohort-wise descriptive quantification, not formal meta-analysis",
        f"n_main_manuscript_cohorts={n_cohorts}",
        f"n_cohorts_PD_gt_HC_by_both_mean_and_median={n_positive_direction_both}",
        f"n_cohorts_positive_nominal_support={n_nominal_positive}",
        f"n_cohorts_directionally_consistent_non_significant_support={n_directional}",
        f"n_cohorts_high_or_moderate_taxon_compatibility={n_high_or_moderate_compat}",
        f"mean_descriptive_evidence_score_0_10={score_mean:.3f}",
        f"median_descriptive_evidence_score_0_10={score_median:.3f}",
        "",
        "Included main-manuscript cohorts:",
    ]

    for _, r in table.iterrows():
        summary.append(
            f"- {r['cohort']}: class={r['cohort_evidence_class']}; "
            f"score={r['descriptive_evidence_score_0_10']:.2f}/10; "
            f"n={int(r['n_samples'])}; PD={int(r['n_PD'])}; HC={int(r['n_HC'])}; "
            f"matched_taxa={int(r['n_matched_curli_taxa'])}; "
            f"mean_direction={r['effect_direction_mean']}; median_direction={r['effect_direction_median']}; "
            f"p={r['mannwhitney_p']}; delta={r['cliffs_delta']}"
        )

    summary += [
        "",
        "Interpretation:",
        "This table quantifies cohort-wise evidence for the Curli Carrier Burden index without pooling effects across studies.",
        "The descriptive evidence score is intended for transparent comparison of direction, nominal support, effect size, taxon compatibility, and sample size.",
        "It should not be interpreted as a formal meta-analytic statistic or diagnostic performance metric.",
    ]

    summary_path.write_text("\n".join(summary) + "\n")
    status_path.write_text("PASS\n")

    print("PASS")
    print(f"[DONE] {table_path}")
    print(f"[DONE] {summary_path}")
    print(f"[DONE] {status_path}")

if __name__ == "__main__":
    main()
