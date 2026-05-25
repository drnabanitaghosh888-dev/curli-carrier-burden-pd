from __future__ import annotations

import py_compile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/phase2j_00_update_cross_cohort_synthesis_with_mao.py"


def test_phase2j_script_exists_and_compiles():
    assert SCRIPT.exists()
    py_compile.compile(str(SCRIPT), doraise=True)


def test_outputs_and_required_interpretations_declared():
    text = SCRIPT.read_text(encoding="utf-8", errors="ignore")
    lower = text.lower()
    assert "cross_cohort_processed_table_summary_updated.tsv" in lower
    assert "cross_cohort_model_comparison_updated.tsv" in lower
    assert "cross_cohort_burden_species_overlap_updated.tsv" in lower
    assert "mao central china" in lower
    assert "DIRECTIONALLY_CONSISTENT_BUT_NOT_SIGNIFICANT" in text
    assert "REPLICATED_DIRECTIONALLY_AND_STATISTICALLY" in text
    assert "NOT_ANALYZED_PROCESSED_TABLES_NOT_IDENTIFIED" in text


def test_required_conclusion_has_no_positive_raw_or_gene_claims():
    text = SCRIPT.read_text(encoding="utf-8", errors="ignore")
    assert "processed-table ecological taxonomic-proxy association" in text
    assert "do not establish csg gene presence" in text
    assert "curli expression" in text
    assert "intact operon architecture" in text
    assert "raw-read validation" in text


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
