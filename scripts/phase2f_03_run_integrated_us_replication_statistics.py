from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IN_BURDEN = ROOT / "results/phase2f_integrated_us_replication/integrated_us_curli_carrier_burden.tsv"
OUT_DIR = ROOT / "results/phase2f_integrated_us_replication"
REP = ROOT / "data/phase2f/reports"
OUT_DIR.mkdir(parents=True, exist_ok=True)
REP.mkdir(parents=True, exist_ok=True)

OUT = OUT_DIR / "integrated_us_curli_replication_statistics.tsv"
PLAN_RUNTIME = REP / "phase2f_integrated_us_statistics_model_plan_runtime.tsv"
STATUS = REP / "phase2f_integrated_us_statistics_status.txt"

MODELS = [
    ("MWU_A", "PD_vs_all_controls", "CurliCarrierBurden", "nonparametric"),
    ("LOGIT_B", "PD_vs_all_controls", "CurliCarrierBurden_log1p", "primary_unadjusted"),
    ("LOGIT_C", "PD_vs_all_controls", "CurliCarrierBurden_log1p", "primary_adjusted"),
    ("LOGIT_D", "PD_vs_all_controls", "CurliCarrierPresence", "secondary_unadjusted"),
    ("LOGIT_E", "PD_vs_all_controls", "CurliCarrierPresence", "secondary_adjusted"),
    ("LOGIT_F", "PD_vs_PC_only", "CurliCarrierBurden_log1p", "sensitivity"),
    ("LOGIT_G", "PD_vs_HC_only", "CurliCarrierBurden_log1p", "sensitivity"),
]


try:
    from scipy.stats import mannwhitneyu  # type: ignore
except Exception:  # pragma: no cover
    mannwhitneyu = None

try:
    import statsmodels.formula.api as smf  # type: ignore
except Exception:  # pragma: no cover
    smf = None


def _na() -> float:
    return float("nan")


def _result_row(model_id: str, comparison: str, predictor: str) -> dict[str, Any]:
    return {
        "model_id": model_id,
        "comparison": comparison,
        "n_used": _na(),
        "predictor": predictor,
        "coefficient": _na(),
        "odds_ratio": _na(),
        "ci_lower": _na(),
        "ci_upper": _na(),
        "p_value": _na(),
        "q_value": _na(),
        "status": "SKIPPED",
        "notes": "",
    }


def _apply_bh_qvalues(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    p = pd.to_numeric(out["p_value"], errors="coerce")
    idx = p[p.notna()].index.tolist()
    if not idx:
        return out

    pvals = p.loc[idx].astype(float).values
    m = len(pvals)
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    adj = np.empty(m, dtype=float)

    prev = 1.0
    for i in range(m - 1, -1, -1):
        rank = i + 1
        q = (sorted_p[i] * m) / rank
        q = min(q, prev, 1.0)
        adj[i] = q
        prev = q

    qvals = np.empty(m, dtype=float)
    qvals[order] = adj
    out.loc[idx, "q_value"] = qvals
    return out


def _subset_for_comparison(df: pd.DataFrame, comparison: str) -> pd.DataFrame:
    if comparison == "PD_vs_all_controls":
        return df[df["Case_status"].isin(["PD", "Control"])].copy()
    if comparison == "PD_vs_PC_only":
        return df[df["donor_group"].isin(["PD", "PC"])].copy()
    if comparison == "PD_vs_HC_only":
        return df[df["donor_group"].isin(["PD", "HC"])].copy()
    return df.copy()


def _run_mwu(df: pd.DataFrame) -> dict[str, Any]:
    row = _result_row("MWU_A", "PD_vs_all_controls", "CurliCarrierBurden")
    if mannwhitneyu is None:
        row["status"] = "FAILED"
        row["notes"] = "scipy_not_available"
        return row

    d = _subset_for_comparison(df, "PD_vs_all_controls")
    d["CurliCarrierBurden"] = pd.to_numeric(d["CurliCarrierBurden"], errors="coerce")
    d = d.dropna(subset=["Case_status", "CurliCarrierBurden"])

    pd_vals = d.loc[d["Case_status"] == "PD", "CurliCarrierBurden"].astype(float)
    ct_vals = d.loc[d["Case_status"] == "Control", "CurliCarrierBurden"].astype(float)

    row["n_used"] = float(len(pd_vals) + len(ct_vals))
    if len(pd_vals) < 2 or len(ct_vals) < 2:
        row["status"] = "FAILED"
        row["notes"] = "insufficient_group_sizes_for_mwu"
        return row

    try:
        res = mannwhitneyu(pd_vals, ct_vals, alternative="two-sided")
        row["coefficient"] = float(np.median(pd_vals) - np.median(ct_vals))
        row["p_value"] = float(res.pvalue)
        row["status"] = "SUCCESS"
        row["notes"] = "median_difference_pd_minus_control"
    except Exception as exc:  # pragma: no cover
        row["status"] = "FAILED"
        row["notes"] = f"mwu_error:{type(exc).__name__}"
    return row


def _run_logit(
    df: pd.DataFrame,
    model_id: str,
    comparison: str,
    predictor: str,
    adjusted: bool,
) -> dict[str, Any]:
    row = _result_row(model_id, comparison, predictor)
    if smf is None:
        row["status"] = "FAILED"
        row["notes"] = "statsmodels_not_available"
        return row

    d = _subset_for_comparison(df, comparison)

    need = ["PD_binary", predictor]
    formula = f"PD_binary ~ {predictor}"
    if adjusted:
        need += ["host_age", "sex", "host_body_mass_index"]
        formula = f"PD_binary ~ {predictor} + host_age + C(sex) + host_body_mass_index"

    for c in ["PD_binary", predictor, "host_age", "host_body_mass_index"]:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")

    d = d.dropna(subset=[c for c in need if c in d.columns]).copy()

    if adjusted and d["sex"].nunique(dropna=True) < 2:
        row["status"] = "SKIPPED"
        row["notes"] = "insufficient_sex_variation_for_adjusted_model"
        row["n_used"] = float(len(d))
        return row

    if d["PD_binary"].nunique(dropna=True) < 2:
        row["status"] = "SKIPPED"
        row["notes"] = "insufficient_outcome_variation"
        row["n_used"] = float(len(d))
        return row

    row["n_used"] = float(len(d))
    if len(d) < 20:
        row["status"] = "SKIPPED"
        row["notes"] = "too_few_rows_after_dropping_missing"
        return row

    try:
        fit = smf.logit(formula=formula, data=d).fit(disp=False)
        if predictor not in fit.params.index:
            row["status"] = "FAILED"
            row["notes"] = "predictor_not_in_model_params"
            return row

        beta = float(fit.params[predictor])
        se = float(fit.bse[predictor])
        pval = float(fit.pvalues[predictor])

        row["coefficient"] = beta
        row["odds_ratio"] = float(math.exp(beta))
        row["ci_lower"] = float(math.exp(beta - 1.96 * se))
        row["ci_upper"] = float(math.exp(beta + 1.96 * se))
        row["p_value"] = pval
        row["status"] = "SUCCESS"
        row["notes"] = ""
    except Exception as exc:
        row["status"] = "FAILED"
        row["notes"] = f"logit_error:{type(exc).__name__}"

    return row


def _write_plan_runtime(execute: bool) -> None:
    with PLAN_RUNTIME.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["model_id", "comparison", "predictor", "role", "planned", "executed"])
        for mid, comp, pred, role in MODELS:
            w.writerow([mid, comp, pred, role, "yes", "yes" if execute else "no"])


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase2F integrated-us replication statistics")
    parser.add_argument("--execute", action="store_true", help="Run statistics and write outputs")
    args = parser.parse_args()

    _write_plan_runtime(execute=args.execute)

    if not args.execute:
        STATUS.write_text("STATUS = DRY_RUN_ONLY\n", encoding="utf-8")
        return

    df = pd.read_csv(IN_BURDEN, sep="\t")

    rows: list[dict[str, Any]] = []
    rows.append(_run_mwu(df))
    rows.append(_run_logit(df, "LOGIT_B", "PD_vs_all_controls", "CurliCarrierBurden_log1p", adjusted=False))
    rows.append(_run_logit(df, "LOGIT_C", "PD_vs_all_controls", "CurliCarrierBurden_log1p", adjusted=True))
    rows.append(_run_logit(df, "LOGIT_D", "PD_vs_all_controls", "CurliCarrierPresence", adjusted=False))
    rows.append(_run_logit(df, "LOGIT_E", "PD_vs_all_controls", "CurliCarrierPresence", adjusted=True))
    rows.append(_run_logit(df, "LOGIT_F", "PD_vs_PC_only", "CurliCarrierBurden_log1p", adjusted=True))
    rows.append(_run_logit(df, "LOGIT_G", "PD_vs_HC_only", "CurliCarrierBurden_log1p", adjusted=True))

    out = pd.DataFrame(rows)
    out = _apply_bh_qvalues(out)

    success = int((out["status"] == "SUCCESS").sum())
    failed = int((out["status"] == "FAILED").sum())
    skipped = int((out["status"] == "SKIPPED").sum())

    if (out["status"] == "NOT_RUN_IN_SCAFFOLD").any():
        final_status = "FAIL"
        note = "not_run_rows_detected"
    elif success == len(MODELS):
        final_status = "PASS"
        note = "all_models_success"
    elif success >= 2 and set(out.loc[out["status"] == "SUCCESS", "model_id"]).issuperset({"MWU_A", "LOGIT_B"}):
        final_status = "PASS_WITH_WARNINGS"
        note = "core_models_success_with_partial_failures"
    else:
        final_status = "FAIL"
        note = "insufficient_successful_models"

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

    STATUS.write_text(
        f"STATUS = {final_status}\n"
        f"N_SUCCESS = {success}\n"
        f"N_FAILED = {failed}\n"
        f"N_SKIPPED = {skipped}\n"
        f"NOTES = {note}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
