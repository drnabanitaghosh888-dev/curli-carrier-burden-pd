from __future__ import annotations

import py_compile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/phase2g_00_generate_cross_cohort_synthesis.py"


def test_phase2g_script_exists_and_compiles():
    assert SCRIPT.exists()
    py_compile.compile(str(SCRIPT), doraise=True)


def test_output_paths_declared():
    t = SCRIPT.read_text(encoding="utf-8", errors="ignore").lower()
    must_have = [
        "cross_cohort_processed_table_summary.tsv",
        "cross_cohort_model_comparison.tsv",
        "cross_cohort_burden_species_overlap.tsv",
        "phase2g_cross_cohort_synthesis_summary.md",
        "phase2g_cross_cohort_final_status.txt",
        "phase2g_cross_cohort_file_manifest.tsv",
    ]
    for s in must_have:
        assert s in t


def test_summary_conclusion_has_processed_table_and_no_overclaims():
    t = SCRIPT.read_text(encoding="utf-8", errors="ignore").lower()
    assert "processed-table" in t
    assert "do not establish csg gene presence" in t
    assert "do not establish" in t
    assert "raw-read validation" in t


def test_no_download_rawread_heavy_or_phase3_commands_in_phase2g_script():
    t = SCRIPT.read_text(encoding="utf-8", errors="ignore").lower()
    forbidden = [
        "wget",
        "prefetch",
        "fasterq-dump",
        "fastq",
        "sra",
        "qiita",
        "phase 3",
        "hmmsearch",
        "hmmbuild",
        "kraken",
        "humann",
        "metaphlan cli",
        "bowtie2",
        "blast",
        "diamond",
        "mafft",
    ]
    for token in forbidden:
        assert token not in t
