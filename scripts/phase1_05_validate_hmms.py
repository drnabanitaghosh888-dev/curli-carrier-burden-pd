#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


GENES = ["csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"]
SUMMARY_COLUMNS = [
    "gene",
    "pos_total",
    "pos_passing",
    "pos_best_match_correct",
    "pos_best_match_wrong",
    "negative_hits_passing_threshold",
    "tp_rate",
    "mismatch_rate",
    "median_best_evalue_correct",
    "median_best_evalue_wrong",
    "coverage_available_fraction",
    "ambiguity_fraction",
    "classification",
    "warning",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Phase 1 HMM validation and reclassification.")
    parser.add_argument("--config", default="config/phase1_config.yml")
    parser.add_argument("--outdir", default="data/phase1/validation")
    parser.add_argument(
        "--reclassify-existing",
        action="store_true",
        help="Rebuild summary from existing hmm_validation_results.tsv without rerunning hmmsearch.",
    )
    return parser.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        obj = yaml.safe_load(handle)
    if not isinstance(obj, dict):
        raise ValueError(f"Invalid YAML mapping: {path}")
    return obj


def _safe_float(x: Any) -> float | None:
    try:
        if x is None:
            return None
        sx = str(x).strip()
        if sx == "" or sx.lower() in {"na", "nan", "none"}:
            return None
        return float(sx)
    except Exception:
        return None


def _coverage_value_and_status(row: dict[str, Any], legacy_all_zero_mode: bool = False) -> tuple[float | None, str]:
    # Coverage can come from coverage_proxy in existing results. It may be unavailable
    # or malformed in historical files; unavailable coverage must not be treated as 0.
    status_raw = str(row.get("coverage_status", "")).strip().lower()
    cov = _safe_float(row.get("coverage_proxy"))
    if cov is None or (isinstance(cov, float) and math.isnan(cov)):
        return None, "not_available"
    if legacy_all_zero_mode and cov == 0.0:
        return None, "not_available"
    # Legacy safeguard: historical reclassification artifacts stored 0.0 coverage
    # for all rows when coverage was not actually computed. Treat as unavailable
    # unless an explicit "available" status exists.
    if cov == 0.0 and status_raw != "available":
        return None, "not_available"
    return cov, "available"


def _median_or_na(values: list[float]) -> str:
    if not values:
        return "NA"
    s = pd.Series(values, dtype="float64").median()
    return f"{float(s):.3e}"


def _log10_margin(best_evalue: float, second_evalue: float | None) -> float | None:
    if second_evalue is None:
        return None
    b = max(best_evalue, 1e-300)
    s = max(second_evalue, 1e-300)
    return math.log10(s) - math.log10(b)


def _classification(tp_rate: float, mismatch_rate: float, neg_hits: int) -> str:
    if tp_rate >= 0.90 and mismatch_rate <= 0.05 and neg_hits == 0:
        return "reliable"
    if tp_rate >= 0.70 and neg_hits == 0:
        return "borderline"
    return "weak"


def build_summary_df(rows: list[dict[str, Any]], evalue_thr: float, coverage_thr: float) -> pd.DataFrame:
    # Positive rows are grouped by expected_gene (required behavior).
    pos_rows = [r for r in rows if str(r.get("control_set", "")).strip().lower() == "positive"]
    neg_rows = [r for r in rows if str(r.get("control_set", "")).strip().lower() == "negative"]

    # Legacy artifact guard: if all parseable positive coverages are exactly 0.0,
    # coverage was likely not computed and should be treated as unavailable.
    pos_cov_values = []
    for r in pos_rows:
        cv = _safe_float(r.get("coverage_proxy"))
        if cv is not None and not math.isnan(cv):
            pos_cov_values.append(cv)
    legacy_all_zero_mode = bool(pos_cov_values) and all(cv == 0.0 for cv in pos_cov_values)

    # Precompute ambiguity by target based on best vs second-best e-value among positive hits.
    pos_df = pd.DataFrame(pos_rows)
    ambiguity_by_target: dict[str, bool] = {}
    if not pos_df.empty and "target_name" in pos_df.columns and "full_evalue" in pos_df.columns:
        for target, grp in pos_df.groupby("target_name", dropna=False):
            evalues = []
            for _, rr in grp.iterrows():
                ev = _safe_float(rr.get("full_evalue"))
                if ev is not None:
                    evalues.append((ev, rr))
            evalues.sort(key=lambda x: x[0])
            if len(evalues) < 2:
                ambiguity_by_target[str(target)] = False
                continue
            best_ev = evalues[0][0]
            second_ev = evalues[1][0]
            margin = _log10_margin(best_ev, second_ev)
            second_strong = second_ev <= evalue_thr
            ambiguity_by_target[str(target)] = bool(
                second_strong and margin is not None and margin < 5.0
            )

    summary_rows: list[dict[str, Any]] = []
    for gene in GENES:
        gene_pos = [
            r
            for r in pos_rows
            if str(r.get("expected_gene", "")).strip() == gene
            and str(r.get("status", "")).startswith("OK")
        ]

        pos_total = len(gene_pos)
        pos_correct = 0
        pos_wrong = 0
        pos_passing = 0
        cov_available = 0
        ambiguity_count = 0
        correct_evals: list[float] = []
        wrong_evals: list[float] = []
        pair_warning_ab = False
        pair_warning_ce = False

        for r in gene_pos:
            predicted = str(r.get("predicted_hmm", "")).strip()
            expected = str(r.get("expected_gene", "")).strip()
            ev = _safe_float(r.get("full_evalue"))
            evalue_pass = (ev is not None) and (ev <= evalue_thr)
            correct_pass = predicted == expected

            cov, cov_status = _coverage_value_and_status(r, legacy_all_zero_mode=legacy_all_zero_mode)
            if cov_status == "available":
                cov_available += 1
            coverage_pass_or_unavailable = (cov is None) or (cov >= coverage_thr)
            positive_pass = bool(evalue_pass and correct_pass and coverage_pass_or_unavailable)

            if correct_pass:
                pos_correct += 1
                if ev is not None:
                    correct_evals.append(ev)
            else:
                pos_wrong += 1
                if ev is not None:
                    wrong_evals.append(ev)

            if positive_pass:
                pos_passing += 1

            tname = str(r.get("target_name", ""))
            if tname and ambiguity_by_target.get(tname, False):
                ambiguity_count += 1

            # Gene-specific pair-aware cross-reactivity warnings
            if gene in {"csgA", "csgB"} and {gene, predicted} == {"csgA", "csgB"} and predicted != gene:
                pair_warning_ab = True
            if gene in {"csgC", "csgE"} and {gene, predicted} == {"csgC", "csgE"} and predicted != gene:
                pair_warning_ce = True

        neg_gene = [
            r
            for r in neg_rows
            if str(r.get("predicted_hmm", "")).strip() == gene
            and str(r.get("status", "")).startswith("OK")
        ]
        neg_hits = 0
        for r in neg_gene:
            ev = _safe_float(r.get("full_evalue"))
            if ev is not None and ev <= evalue_thr:
                neg_hits += 1

        tp_rate = (pos_passing / pos_total) if pos_total else 0.0
        mismatch_rate = (pos_wrong / pos_total) if pos_total else 0.0
        cov_frac = (cov_available / pos_total) if pos_total else 0.0
        amb_frac = (ambiguity_count / pos_total) if pos_total else 0.0
        klass = _classification(tp_rate, mismatch_rate, neg_hits)

        warnings = []
        if pair_warning_ab:
            warnings.append("csgA_csgB_crossreactivity")
        if pair_warning_ce:
            warnings.append("csgC_csgE_crossreactivity")
        if cov_frac < 1.0:
            warnings.append("coverage_not_fully_available")

        summary_rows.append(
            {
                "gene": gene,
                "pos_total": int(pos_total),
                "pos_passing": int(pos_passing),
                "pos_best_match_correct": int(pos_correct),
                "pos_best_match_wrong": int(pos_wrong),
                "negative_hits_passing_threshold": int(neg_hits),
                "tp_rate": float(tp_rate),
                "mismatch_rate": float(mismatch_rate),
                "median_best_evalue_correct": _median_or_na(correct_evals),
                "median_best_evalue_wrong": _median_or_na(wrong_evals),
                "coverage_available_fraction": float(cov_frac),
                "ambiguity_fraction": float(amb_frac),
                "classification": klass,
                "warning": ";".join(warnings),
            }
        )

    df = pd.DataFrame(summary_rows)
    df = df[SUMMARY_COLUMNS]
    return df


def reclassify_existing(results_path: Path, summary_path: Path, evalue_thr: float, coverage_thr: float) -> int:
    if not results_path.exists():
        raise FileNotFoundError(f"Missing results file: {results_path}")
    # Read raw rows, ignore stale pass_thresholds from previous runs.
    rows = pd.read_csv(results_path, sep="\t", dtype=str).fillna("").to_dict(orient="records")
    df = build_summary_df(rows, evalue_thr=evalue_thr, coverage_thr=coverage_thr)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(summary_path, sep="\t", index=False)
    print(f"[DONE] reclassified summary written: {summary_path}")
    return 0


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg = load_yaml(root / args.config)
    outdir = (root / args.outdir).resolve()
    results_path = outdir / "hmm_validation_results.tsv"
    summary_path = outdir / "hmm_validation_summary.tsv"
    evalue_thr = float(cfg.get("validation", {}).get("evalue_threshold", 1e-20))
    coverage_thr = float(cfg.get("validation", {}).get("coverage_threshold", 0.50))

    # Coverage note:
    # - unavailable coverage is never treated as zero.
    # - preliminary validation can use strong e-value + correct gene identity
    #   when negative controls are clean.
    # - later metagenomic screening can enforce stricter coverage thresholds.

    if args.reclassify_existing:
        return reclassify_existing(results_path, summary_path, evalue_thr, coverage_thr)

    raise RuntimeError(
        "This patched script is intended for --reclassify-existing mode in this workflow step. "
        "Full hmmsearch execution is intentionally not run here."
    )


if __name__ == "__main__":
    raise SystemExit(main())
    # Legacy artifact guard: if all parseable positive coverages are exactly 0.0,
    # coverage was likely not computed and should be treated as unavailable.
    pos_cov_values = []
    for r in pos_rows:
        cv = _safe_float(r.get("coverage_proxy"))
        if cv is not None and not math.isnan(cv):
            pos_cov_values.append(cv)
    legacy_all_zero_mode = bool(pos_cov_values) and all(cv == 0.0 for cv in pos_cov_values)
