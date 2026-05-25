from __future__ import annotations

import py_compile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/phase2i_06_run_mao_replication_statistics.py"


def test_mao_statistics_script_exists_and_compiles():
    assert SCRIPT.exists()
    py_compile.compile(str(SCRIPT), doraise=True)


def test_statistics_logic_is_implemented():
    text = SCRIPT.read_text(encoding="utf-8", errors="ignore").lower()
    assert "curlicarrierburden_log1p" in text
    assert "mannwhitneyu" in text
    assert "smf.logit" in text
    assert "wilcoxon" in text
    assert "_bh_qvalues" in text
    assert "curlicarrierpresence has no variation" in text
    assert "pd_vs_sp_matched_pairs" in text
    assert "source_table2_clean_reads_ratio_pct" in text


def test_dry_run_and_execute_gate_are_present():
    text = SCRIPT.read_text(encoding="utf-8", errors="ignore").lower()
    assert "if not args.execute" in text
    assert "status = dry_run_only" in text
    assert "mao_curli_replication_statistics.tsv" in text
    assert "phase2i_mao_statistics_model_plan_runtime.tsv" in text
    assert "phase2i_mao_statistics_status.txt" in text


def test_no_scaffold_placeholder_success_rows():
    text = SCRIPT.read_text(encoding="utf-8", errors="ignore").lower()
    assert "not_run_in_scaffold" not in text


def test_no_download_rawread_heavy_or_phase3_command_patterns():
    text = SCRIPT.read_text(encoding="utf-8", errors="ignore").lower()
    forbidden = [
        "wget ",
        "curl ",
        "prefetch ",
        "fasterq-dump ",
        "hmmsearch ",
        "hmmbuild ",
        "kraken ",
        "humann ",
        "metaphlan ",
        "bowtie2 ",
        "blast ",
        "diamond ",
        "mafft ",
        "python scripts/phase3",
    ]
    for token in forbidden:
        assert token not in text
