from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IN_BURDEN = ROOT / "results/phase2i_mao_replication/mao_curli_carrier_burden.tsv"
OUT_DIR = ROOT / "results/phase2i_mao_replication"
REP = ROOT / "data/phase2i/reports"
OUT_DIR.mkdir(parents=True, exist_ok=True)
REP.mkdir(parents=True, exist_ok=True)

OUT = OUT_DIR / "mao_curli_replication_statistics.tsv"
PLAN_RUNTIME = REP / "phase2i_mao_statistics_model_plan_runtime.tsv"
STATUS = REP / "phase2i_mao_statistics_status.txt"

MODELS = [
    ("MWU_A", "PD_vs_SP_controls", "CurliCarrierBurden", "mann_whitney_u"),
    ("LOGIT_B", "PD_vs_SP_controls", "CurliCarrierBurden_log1p", "primary_unadjusted_logit"),
    ("PAIRED_C", "PD_vs_SP_matched_pairs", "CurliCarrierBurden", "paired_wilcoxon"),
    ("PRESENCE_D", "PD_vs_SP_controls", "CurliCarrierPresence", "presence_model_expected_skip_if_constant"),
    ("LOGIT_E", "PD_vs_SP_controls", "CurliCarrierBurden_log1p", "sequencing_adjusted_sensitivity"),
]

SEQUENCING_COVARIATES = [
    "source_Table2_Clean_Reads_Ratio_pct",
    "source_Table2_Clean_ReadsEliminating_host_reads",
    "source_Table2_Raw_Reads",
]

try:
    from scipy.stats import mannwhitneyu, wilcoxon  # type: ignore
except Exception:  # pragma: no cover
    mannwhitneyu = None
    wilcoxon = None

try:
    import statsmodels.formula.api as smf  # type: ignore
except Exception:  # pragma: no cover
    smf = None


def _empty_row(model_id: str, comparison: str, predictor: str) -> dict[str, Any]:
    return {
        "model_id": model_id,
        "comparison": comparison,
        "n_used": np.nan,
        "predictor": predictor,
        "coefficient": np.nan,
        "odds_ratio": np.nan,
        "ci_lower": np.nan,
        "ci_upper": np.nan,
        "p_value": np.nan,
        "q_value": np.nan,
        "status": "SKIPPED",
        "notes": "",
    }


def _bh_qvalues(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    p = pd.to_numeric(out["p_value"], errors="coerce")
    valid = p[p.notna()].index.tolist()
    if not valid:
        return out

    pvals = p.loc[valid].astype(float).to_numpy()
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    m = len(sorted_p)
    adjusted = np.empty(m, dtype=float)
    running = 1.0
    for i in range(m - 1, -1, -1):
        rank = i + 1
        q = min((sorted_p[i] * m) / rank, running, 1.0)
        adjusted[i] = q
        running = q

    qvals = np.empty(m, dtype=float)
    qvals[order] = adjusted
    out.loc[valid, "q_value"] = qvals
    return out


def _write_plan_runtime(execute: bool) -> None:
    with PLAN_RUNTIME.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["model_id", "comparison", "predictor", "role", "planned", "executed"])
        for model_id, comparison, predictor, role in MODELS:
            writer.writerow([model_id, comparison, predictor, role, "yes", "yes" if execute else "no"])


def _run_mwu(df: pd.DataFrame) -> dict[str, Any]:
    row = _empty_row("MWU_A", "PD_vs_SP_controls", "CurliCarrierBurden")
    if mannwhitneyu is None:
        row["status"] = "FAILED"
        row["notes"] = "scipy_not_available"
        return row

    d = df[["Case_status", "CurliCarrierBurden"]].copy()
    d["CurliCarrierBurden"] = pd.to_numeric(d["CurliCarrierBurden"], errors="coerce")
    d = d.dropna()
    pd_vals = d.loc[d["Case_status"].astype(str) == "PD", "CurliCarrierBurden"].astype(float)
    control_vals = d.loc[d["Case_status"].astype(str) == "Control", "CurliCarrierBurden"].astype(float)
    row["n_used"] = int(len(pd_vals) + len(control_vals))

    if len(pd_vals) < 2 or len(control_vals) < 2:
        row["status"] = "FAILED"
        row["notes"] = "insufficient_group_sizes"
        return row

    try:
        result = mannwhitneyu(pd_vals, control_vals, alternative="two-sided")
        row["coefficient"] = float(np.median(pd_vals) - np.median(control_vals))
        row["p_value"] = float(result.pvalue)
        row["status"] = "SUCCESS"
        row["notes"] = "median_difference_pd_minus_control"
    except Exception as exc:
        row["status"] = "FAILED"
        row["notes"] = f"mwu_error:{type(exc).__name__}"
    return row


def _run_logit(
    df: pd.DataFrame,
    model_id: str,
    comparison: str,
    predictor: str,
    covariates: list[str] | None = None,
) -> dict[str, Any]:
    row = _empty_row(model_id, comparison, predictor)
    if smf is None:
        row["status"] = "FAILED"
        row["notes"] = "statsmodels_not_available"
        return row

    covariates = covariates or []
    columns = ["PD_binary", predictor] + covariates
    missing = [c for c in columns if c not in df.columns]
    if missing:
        row["status"] = "SKIPPED"
        row["notes"] = "missing_columns:" + ",".join(missing)
        return row

    d = df[columns].copy()
    for col in ["PD_binary", predictor] + covariates:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d = d.dropna()
    row["n_used"] = int(len(d))

    if len(d) < 20:
        row["status"] = "SKIPPED"
        row["notes"] = "too_few_rows_after_missing_drop"
        return row
    if d["PD_binary"].nunique() < 2:
        row["status"] = "SKIPPED"
        row["notes"] = "outcome_has_no_variation"
        return row
    if d[predictor].nunique() < 2:
        row["status"] = "SKIPPED"
        row["notes"] = "predictor_has_no_variation"
        return row

    terms = [predictor] + covariates
    formula = "PD_binary ~ " + " + ".join(terms)
    try:
        fit = smf.logit(formula=formula, data=d).fit(disp=False)
        beta = float(fit.params[predictor])
        se = float(fit.bse[predictor])
        row["coefficient"] = beta
        row["odds_ratio"] = float(math.exp(beta))
        row["ci_lower"] = float(math.exp(beta - 1.96 * se))
        row["ci_upper"] = float(math.exp(beta + 1.96 * se))
        row["p_value"] = float(fit.pvalues[predictor])
        row["status"] = "SUCCESS"
        row["notes"] = "logistic_regression"
    except Exception as exc:
        row["status"] = "FAILED"
        row["notes"] = f"logit_error:{type(exc).__name__}"
    return row


def _pair_id(sample_name: str) -> str:
    match = re.match(r"^(?:PD|SP)_(\d+)$", str(sample_name))
    return match.group(1) if match else ""


def _run_paired_wilcoxon(df: pd.DataFrame) -> dict[str, Any]:
    row = _empty_row("PAIRED_C", "PD_vs_SP_matched_pairs", "CurliCarrierBurden")
    if wilcoxon is None:
        row["status"] = "FAILED"
        row["notes"] = "scipy_not_available"
        return row

    d = df[["sample_name", "donor_group", "CurliCarrierBurden"]].copy()
    d["pair_id"] = d["sample_name"].map(_pair_id)
    d["CurliCarrierBurden"] = pd.to_numeric(d["CurliCarrierBurden"], errors="coerce")
    d = d[(d["pair_id"] != "") & d["donor_group"].isin(["PD", "SP"])].dropna()

    pivot = d.pivot_table(index="pair_id", columns="donor_group", values="CurliCarrierBurden", aggfunc="first")
    if not {"PD", "SP"}.issubset(pivot.columns):
        row["status"] = "SKIPPED"
        row["notes"] = "could_not_infer_complete_pd_sp_pairs"
        return row

    pairs = pivot.dropna(subset=["PD", "SP"])
    row["n_used"] = int(len(pairs) * 2)
    if len(pairs) < 2:
        row["status"] = "SKIPPED"
        row["notes"] = "insufficient_complete_pd_sp_pairs"
        return row

    diffs = pairs["PD"].astype(float) - pairs["SP"].astype(float)
    row["coefficient"] = float(np.median(diffs))
    try:
        result = wilcoxon(diffs, alternative="two-sided", zero_method="wilcox")
        row["p_value"] = float(result.pvalue)
        row["status"] = "SUCCESS"
        row["notes"] = "paired_pd_minus_sp_wilcoxon"
    except ValueError as exc:
        row["status"] = "SKIPPED"
        row["notes"] = f"wilcoxon_not_informative:{type(exc).__name__}"
    except Exception as exc:
        row["status"] = "FAILED"
        row["notes"] = f"wilcoxon_error:{type(exc).__name__}"
    return row


def _suitable_sequencing_covariates(df: pd.DataFrame) -> list[str]:
    usable: list[str] = []
    for col in SEQUENCING_COVARIATES:
        if col not in df.columns:
            continue
        vals = pd.to_numeric(df[col], errors="coerce")
        table = pd.DataFrame({"Case_status": df["Case_status"].astype(str), col: vals}).dropna()
        has_both_groups = set(table["Case_status"]) >= {"PD", "Control"}
        enough_values = table.groupby("Case_status")[col].size().min() >= 10 if has_both_groups else False
        has_variation = table[col].nunique() >= 2
        if has_both_groups and enough_values and has_variation:
            usable.append(col)
    return usable


def _presence_row(df: pd.DataFrame) -> dict[str, Any]:
    row = _empty_row("PRESENCE_D", "PD_vs_SP_controls", "CurliCarrierPresence")
    if "CurliCarrierPresence" not in df.columns:
        row["notes"] = "presence_column_missing"
        return row
    vals = pd.to_numeric(df["CurliCarrierPresence"], errors="coerce").dropna()
    row["n_used"] = int(len(vals))
    if vals.nunique() < 2:
        row["status"] = "SKIPPED"
        row["notes"] = "CurliCarrierPresence has no variation; all samples have burden > 0"
        return row
    return _run_logit(df, "PRESENCE_D", "PD_vs_SP_controls", "CurliCarrierPresence")


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase2I Mao replication statistics")
    parser.add_argument("--execute", action="store_true", help="Compute and write statistics")
    args = parser.parse_args()

    _write_plan_runtime(execute=args.execute)
    if not args.execute:
        STATUS.write_text("STATUS = DRY_RUN_ONLY\n", encoding="utf-8")
        return

    df = pd.read_csv(IN_BURDEN, sep="\t")
    rows: list[dict[str, Any]] = [
        _run_mwu(df),
        _run_logit(df, "LOGIT_B", "PD_vs_SP_controls", "CurliCarrierBurden_log1p"),
        _run_paired_wilcoxon(df),
        _presence_row(df),
    ]

    covariates = _suitable_sequencing_covariates(df)
    if covariates:
        rows.append(
            _run_logit(
                df,
                "LOGIT_E",
                "PD_vs_SP_controls",
                "CurliCarrierBurden_log1p",
                covariates=covariates,
            )
        )
    else:
        row = _empty_row("LOGIT_E", "PD_vs_SP_controls", "CurliCarrierBurden_log1p")
        row["status"] = "SKIPPED"
        row["notes"] = "no_suitable_numeric_sequencing_covariates"
        rows.append(row)

    out = _bh_qvalues(pd.DataFrame(rows))
    out = out[
        [
            "model_id",
            "comparison",
            "n_used",
            "predictor",
            "coefficient",
            "odds_ratio",
            "ci_lower",
            "ci_upper",
            "p_value",
            "q_value",
            "status",
            "notes",
        ]
    ]
    out.to_csv(OUT, sep="\t", index=False)

    successes = set(out.loc[out["status"] == "SUCCESS", "model_id"])
    n_failed = int((out["status"] == "FAILED").sum())
    n_skipped = int((out["status"] == "SKIPPED").sum())
    if {"MWU_A", "LOGIT_B"}.issubset(successes):
        final = "PASS" if n_failed == 0 and n_skipped == 0 else "PASS_WITH_WARNINGS"
    elif successes:
        final = "FAIL"
    else:
        final = "FAIL"

    STATUS.write_text(
        f"STATUS = {final}\n"
        f"N_SUCCESS = {int((out['status'] == 'SUCCESS').sum())}\n"
        f"N_FAILED = {n_failed}\n"
        f"N_SKIPPED = {n_skipped}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
