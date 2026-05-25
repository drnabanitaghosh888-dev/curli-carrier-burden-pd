from __future__ import annotations

import py_compile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_phase2i_scripts_exist_and_compile():
    s1 = ROOT / "scripts/phase2i_01_extract_mao_standardized_tables.py"
    s2 = ROOT / "scripts/phase2i_02_validate_mao_standardized_tables.py"
    assert s1.exists()
    assert s2.exists()
    py_compile.compile(str(s1), doraise=True)
    py_compile.compile(str(s2), doraise=True)


def test_source_path_uses_manual_verified_and_not_supplementary_exact():
    t = (ROOT / "scripts/phase2i_01_extract_mao_standardized_tables.py").read_text(encoding="utf-8", errors="ignore").lower()
    assert "supplementary_manual_verified" in t
    assert "supplementary_exact" not in t


def test_execute_gate_and_pd_sp_mapping_and_outputs_present():
    t = (ROOT / "scripts/phase2i_01_extract_mao_standardized_tables.py").read_text(encoding="utf-8", errors="ignore").lower()
    must = [
        "if not args.execute",
        "pd_",
        "sp_",
        "case_status",
        "pd_binary",
        "donor_group",
        "sample_metadata.tsv",
        "species_abundance.tsv",
        "genus_abundance.tsv",
        "phase2i_mao_standardized_extraction_plan.tsv",
        "phase2i_mao_standardized_extraction_status.txt",
    ]
    for m in must:
        assert m in t


def test_no_download_rawread_heavy_or_phase3_commands_present():
    texts = "\n".join(
        [
            (ROOT / "scripts/phase2i_01_extract_mao_standardized_tables.py").read_text(encoding="utf-8", errors="ignore").lower(),
            (ROOT / "scripts/phase2i_02_validate_mao_standardized_tables.py").read_text(encoding="utf-8", errors="ignore").lower(),
        ]
    )
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
        "metaphlan",
        "bowtie2",
        "blast",
        "diamond",
        "mafft",
    ]
    for token in forbidden:
        assert token not in texts
