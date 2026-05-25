from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/phase1_05_validate_hmms.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("phase1_05_validate_hmms", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def test_groups_by_expected_gene_not_predicted_gene():
    m = _load_module()
    rows = [
        {
            "gene": "csgA",
            "control_set": "positive",
            "status": "OK",
            "target_name": "t1",
            "predicted_hmm": "csgB",
            "expected_gene": "csgA",
            "full_evalue": "1e-50",
            "coverage_proxy": "NA",
            "pass_thresholds": "false",
        },
        {
            "gene": "csgB",
            "control_set": "positive",
            "status": "OK",
            "target_name": "t2",
            "predicted_hmm": "csgA",
            "expected_gene": "csgB",
            "full_evalue": "1e-40",
            "coverage_proxy": "NA",
            "pass_thresholds": "false",
        },
    ]
    df = m.build_summary_df(rows, evalue_thr=1e-20, coverage_thr=0.5)
    by_gene = {r["gene"]: r for r in df.to_dict(orient="records")}
    assert by_gene["csgA"]["pos_total"] == 1
    assert by_gene["csgB"]["pos_total"] == 1


def test_stale_pass_thresholds_is_ignored_and_recomputed():
    m = _load_module()
    rows = [
        {
            "control_set": "positive",
            "status": "OK",
            "target_name": "t1",
            "predicted_hmm": "csgD",
            "expected_gene": "csgD",
            "full_evalue": "1e-80",
            "coverage_proxy": "NA",
            "pass_thresholds": "false",  # stale
        }
    ]
    df = m.build_summary_df(rows, evalue_thr=1e-20, coverage_thr=0.5)
    row = df[df["gene"] == "csgD"].iloc[0]
    assert int(row["pos_total"]) == 1
    assert int(row["pos_passing"]) == 1
    assert float(row["tp_rate"]) == 1.0


def test_unavailable_coverage_does_not_fail_strong_correct_hit():
    m = _load_module()
    rows = [
        {
            "control_set": "positive",
            "status": "OK",
            "target_name": "t1",
            "predicted_hmm": "csgG",
            "expected_gene": "csgG",
            "full_evalue": "1e-60",
            "coverage_proxy": "not_a_number",
            "pass_thresholds": "false",
        }
    ]
    df = m.build_summary_df(rows, evalue_thr=1e-20, coverage_thr=0.5)
    row = df[df["gene"] == "csgG"].iloc[0]
    assert int(row["pos_passing"]) == 1
    assert float(row["coverage_available_fraction"]) == 0.0


def test_dataframe_columns_exact_and_warnings_gene_specific():
    m = _load_module()
    rows = [
        {
            "control_set": "positive",
            "status": "OK",
            "target_name": "tA",
            "predicted_hmm": "csgB",
            "expected_gene": "csgA",
            "full_evalue": "1e-45",
            "coverage_proxy": "NA",
            "pass_thresholds": "false",
        },
        {
            "control_set": "positive",
            "status": "OK",
            "target_name": "tB",
            "predicted_hmm": "csgC",
            "expected_gene": "csgC",
            "full_evalue": "1e-60",
            "coverage_proxy": "NA",
            "pass_thresholds": "false",
        },
    ]
    df = m.build_summary_df(rows, evalue_thr=1e-20, coverage_thr=0.5)
    assert list(df.columns) == m.SUMMARY_COLUMNS
    by_gene = {r["gene"]: r for r in df.to_dict(orient="records")}
    assert "csgA_csgB_crossreactivity" in str(by_gene["csgA"]["warning"])
    assert "csgA_csgB_crossreactivity" not in str(by_gene["csgD"]["warning"])

