#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run Wallen discovery statistics (processed-table proxy analysis).")
    p.add_argument("--execute", action="store_true", help="Run statistics and write result tables.")
    return p.parse_args()


def _load_cfg(root: Path) -> dict:
    with (root / "config/phase2d_config.yml").open("r", encoding="utf-8") as h:
        return yaml.safe_load(h)


def _bh_adjust(pvals: list[float]) -> list[float]:
    if not pvals:
        return []
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [1.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i_rank = n - rank + 1
        val = min(prev, pvals[idx] * n / i_rank)
        adj[idx] = val
        prev = val
    return adj


def _prepare_analysis_frame(burden_df: pd.DataFrame, metadata_df: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    # If burden already has key phenotype columns, use it directly to avoid duplicate *_x/*_y columns.
    use_cols = [sample_col, "Case_status", "PD_binary", "CurliCarrierBurden"]
    if all(c in burden_df.columns for c in use_cols):
        df = burden_df.copy()
        for c in metadata_df.columns:
            if c not in df.columns:
                df = df.merge(metadata_df[[sample_col, c]], on=sample_col, how="left")
    else:
        df = metadata_df.merge(burden_df, on=sample_col, how="inner")
    return df


def _add_burden_features(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["CurliCarrierBurden"] = pd.to_numeric(out["CurliCarrierBurden"], errors="coerce")
    out["CurliCarrierBurden_log1p"] = np.log1p(out["CurliCarrierBurden"].clip(lower=0))
    out["CurliCarrierPresence"] = (out["CurliCarrierBurden"] > 0).astype(int)
    return out


def _safe_ci(model_result, term: str) -> tuple[float, float]:
    ci = model_result.conf_int()
    if term not in ci.index:
        return (math.nan, math.nan)
    return float(ci.loc[term, 0]), float(ci.loc[term, 1])


def _fit_logit(df: pd.DataFrame, outcome: str, predictor: str, covariates: list[str], model_id: str) -> dict:
    import statsmodels.formula.api as smf

    cols = [outcome, predictor] + covariates
    d = df[cols].copy()
    if "Sex" in d.columns:
        d["Sex"] = d["Sex"].astype(str)
    d = d.dropna()
    if outcome == "PD_binary":
        d = d[d[outcome].isin([0, 1])]
    if d.empty:
        return {
            "model_id": model_id,
            "n_used": 0,
            "predictor": predictor,
            "coefficient": math.nan,
            "odds_ratio": math.nan,
            "ci_lower": math.nan,
            "ci_upper": math.nan,
            "p_value": math.nan,
            "q_value": math.nan,
            "status": "SKIPPED",
            "notes": "no complete-case rows available",
        }

    terms = [predictor] + covariates
    formula_terms = [predictor] + [("C(Sex)" if c == "Sex" else c) for c in covariates]
    formula = f"{outcome} ~ " + " + ".join(formula_terms)
    model = smf.logit(formula=formula, data=d).fit(disp=False, maxiter=200)
    term_name = predictor
    coef = float(model.params.get(term_name, math.nan))
    pval = float(model.pvalues.get(term_name, math.nan))
    ci_l, ci_u = _safe_ci(model, term_name)
    return {
        "model_id": model_id,
        "n_used": int(d.shape[0]),
        "predictor": predictor,
        "coefficient": coef,
        "odds_ratio": float(np.exp(coef)) if not math.isnan(coef) else math.nan,
        "ci_lower": float(np.exp(ci_l)) if not math.isnan(ci_l) else math.nan,
        "ci_upper": float(np.exp(ci_u)) if not math.isnan(ci_u) else math.nan,
        "p_value": pval,
        "q_value": math.nan,
        "status": "PASS",
        "notes": f"logistic regression with terms: {', '.join(terms)}",
    }


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg = _load_cfg(root)
    in_burden = root / cfg["output_dirs"]["results"] / "wallen_curli_carrier_burden.tsv"
    in_meta = root / cfg["output_dirs"]["results"] / "wallen_analysis_metadata.tsv"
    out_results = root / cfg["output_dirs"]["results"]
    out_reports = root / cfg["output_dirs"]["reports"]
    out_results.mkdir(parents=True, exist_ok=True)
    out_reports.mkdir(parents=True, exist_ok=True)
    out_stats = out_results / "wallen_curli_burden_statistics.tsv"
    out_status = out_reports / "phase2d_statistics_status.txt"
    out_runtime_plan = out_reports / "phase2d_statistics_model_plan_runtime.tsv"

    planned = [
        {
            "model_id": "MW_PD_vs_Control_raw_burden",
            "predictor": "CurliCarrierBurden",
            "model_type": "mann_whitney_u",
            "status": "PLANNED",
            "notes": "non-parametric group comparison on raw burden",
        },
        {
            "model_id": "LOGIT_PD_log1p_unadjusted",
            "predictor": "CurliCarrierBurden_log1p",
            "model_type": "logistic_regression",
            "status": "PLANNED",
            "notes": "primary regression predictor is log1p burden",
        },
        {
            "model_id": "LOGIT_PD_log1p_adjusted_age_sex_BMI",
            "predictor": "CurliCarrierBurden_log1p",
            "model_type": "logistic_regression",
            "status": "PLANNED",
            "notes": "adjusted primary model",
        },
        {
            "model_id": "LOGIT_PD_presence_unadjusted",
            "predictor": "CurliCarrierPresence",
            "model_type": "logistic_regression",
            "status": "PLANNED",
            "notes": "secondary binary exposure model",
        },
        {
            "model_id": "LOGIT_PD_presence_adjusted_age_sex_BMI",
            "predictor": "CurliCarrierPresence",
            "model_type": "logistic_regression",
            "status": "PLANNED",
            "notes": "secondary adjusted binary exposure model",
        },
        {
            "model_id": "LOGIT_constipation_log1p_exploratory",
            "predictor": "CurliCarrierBurden_log1p",
            "model_type": "logistic_regression_exploratory",
            "status": "PLANNED_IF_AVAILABLE",
            "notes": "run only if constipation variable can be mapped to binary",
        },
    ]
    pd.DataFrame(planned).to_csv(out_runtime_plan, sep="\t", index=False)

    if not args.execute:
        out_status.write_text("DRY_RUN_ONLY\n", encoding="utf-8")
        print("[INFO] Dry-run only. Use --execute to run statistics.")
        print(f"[DONE] {out_runtime_plan}")
        print(f"[DONE] {out_status}")
        return 0

    missing = [str(p) for p in [in_burden, in_meta] if not p.exists()]
    if missing:
        print("[ERROR] Missing required input(s):")
        for m in missing:
            print(m)
        return 1

    burden_df = pd.read_csv(in_burden, sep="\t")
    meta_df = pd.read_csv(in_meta, sep="\t")
    sample_col = cfg["sample_id_column"]
    df = _prepare_analysis_frame(burden_df, meta_df, sample_col)
    df = _add_burden_features(df)
    case_col = cfg["case_control_column"]
    case = cfg["case_label"]
    ctrl = cfg["control_label"]

    if "PD_binary" not in df.columns:
        df["PD_binary"] = df[case_col].map({case: 1, ctrl: 0})
    df["PD_binary"] = pd.to_numeric(df["PD_binary"], errors="coerce")

    pd_vals = pd.to_numeric(df.loc[df[case_col] == case, "CurliCarrierBurden"], errors="coerce").dropna()
    ct_vals = pd.to_numeric(df.loc[df[case_col] == ctrl, "CurliCarrierBurden"], errors="coerce").dropna()
    delta_median = float(pd_vals.median() - ct_vals.median()) if len(pd_vals) and len(ct_vals) else math.nan

    # Mann-Whitney U for raw burden
    from scipy.stats import mannwhitneyu

    try:
        mw = mannwhitneyu(pd_vals, ct_vals, alternative="two-sided")
        mw_p = float(mw.pvalue)
        mw_status = "PASS"
        mw_notes = "raw burden non-parametric comparison"
    except Exception as exc:  # pragma: no cover
        mw_p = math.nan
        mw_status = "FAIL"
        mw_notes = f"mann-whitney failed: {exc}"

    rows = []
    rows.append(
        {
            "model_id": "MW_PD_vs_Control_raw_burden",
            "n_used": int(len(pd_vals) + len(ct_vals)),
            "predictor": "CurliCarrierBurden",
            "coefficient": delta_median,
            "odds_ratio": math.nan,
            "ci_lower": math.nan,
            "ci_upper": math.nan,
            "p_value": mw_p,
            "q_value": math.nan,
            "status": mw_status,
            "notes": mw_notes,
        }
    )

    # Logistic models with log1p burden and presence
    model_specs = [
        ("LOGIT_PD_log1p_unadjusted", "PD_binary", "CurliCarrierBurden_log1p", []),
        ("LOGIT_PD_log1p_adjusted_age_sex_BMI", "PD_binary", "CurliCarrierBurden_log1p", ["Age_at_collection", "Sex", "BMI"]),
        ("LOGIT_PD_presence_unadjusted", "PD_binary", "CurliCarrierPresence", []),
        ("LOGIT_PD_presence_adjusted_age_sex_BMI", "PD_binary", "CurliCarrierPresence", ["Age_at_collection", "Sex", "BMI"]),
    ]
    for mid, outcome, predictor, covs in model_specs:
        try:
            rows.append(_fit_logit(df, outcome, predictor, covs, mid))
        except Exception as exc:  # pragma: no cover
            rows.append(
                {
                    "model_id": mid,
                    "n_used": 0,
                    "predictor": predictor,
                    "coefficient": math.nan,
                    "odds_ratio": math.nan,
                    "ci_lower": math.nan,
                    "ci_upper": math.nan,
                    "p_value": math.nan,
                    "q_value": math.nan,
                    "status": "FAIL",
                    "notes": f"logistic fit failed: {exc}",
                }
            )

    # Exploratory constipation model if a suitable variable exists
    constipation_col = "Day_of_stool_collection_constipation"
    if constipation_col in df.columns:
        dd = df.copy()
        cc = pd.to_numeric(dd[constipation_col], errors="coerce")
        if cc.notna().sum() > 0 and set(cc.dropna().unique()).issubset({0, 1}):
            dd["constipation_binary"] = cc
            try:
                rows.append(
                    _fit_logit(
                        dd,
                        "constipation_binary",
                        "CurliCarrierBurden_log1p",
                        ["Age_at_collection", "Sex", "BMI"],
                        "LOGIT_constipation_log1p_exploratory",
                    )
                )
            except Exception as exc:  # pragma: no cover
                rows.append(
                    {
                        "model_id": "LOGIT_constipation_log1p_exploratory",
                        "n_used": 0,
                        "predictor": "CurliCarrierBurden_log1p",
                        "coefficient": math.nan,
                        "odds_ratio": math.nan,
                        "ci_lower": math.nan,
                        "ci_upper": math.nan,
                        "p_value": math.nan,
                        "q_value": math.nan,
                        "status": "FAIL",
                        "notes": f"exploratory constipation model failed: {exc}",
                    }
                )
        else:
            rows.append(
                {
                    "model_id": "LOGIT_constipation_log1p_exploratory",
                    "n_used": 0,
                    "predictor": "CurliCarrierBurden_log1p",
                    "coefficient": math.nan,
                    "odds_ratio": math.nan,
                    "ci_lower": math.nan,
                    "ci_upper": math.nan,
                    "p_value": math.nan,
                    "q_value": math.nan,
                    "status": "SKIPPED",
                    "notes": "constipation variable not binary/interpretable",
                }
            )
    else:
        rows.append(
            {
                "model_id": "LOGIT_constipation_log1p_exploratory",
                "n_used": 0,
                "predictor": "CurliCarrierBurden_log1p",
                "coefficient": math.nan,
                "odds_ratio": math.nan,
                "ci_lower": math.nan,
                "ci_upper": math.nan,
                "p_value": math.nan,
                "q_value": math.nan,
                "status": "SKIPPED",
                "notes": "constipation variable unavailable",
            }
        )

    out_df = pd.DataFrame(
        rows,
        columns=[
            "model_id",
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
        ],
    )
    valid_idx = out_df.index[out_df["p_value"].apply(lambda x: pd.notna(x) and not isinstance(x, str))]
    valid_p = [float(out_df.loc[i, "p_value"]) for i in valid_idx]
    if valid_p:
        qvals = _bh_adjust(valid_p)
        for i, q in zip(valid_idx, qvals):
            out_df.loc[i, "q_value"] = q

    out_df.to_csv(out_stats, sep="\t", index=False)
    status = "PASS"
    if (out_df["status"] == "FAIL").any():
        status = "PASS_WITH_WARNINGS"
    out_status.write_text(
        status
        + "\nprocessed-table proxy statistics written; log1p burden is the primary regression predictor and presence is secondary; raw-read validation still required.\n",
        encoding="utf-8",
    )
    print(f"[DONE] {out_stats}")
    print(f"[DONE] {out_runtime_plan}")
    print(f"[DONE] {out_status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
