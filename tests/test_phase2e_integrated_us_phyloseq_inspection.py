from __future__ import annotations

import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_script_exists_and_compiles():
    script = ROOT / "scripts/phase2e_04_inspect_integrated_us_phyloseq_object.py"
    assert script.exists()
    py_compile.compile(str(script), doraise=True)


def test_script_declares_required_outputs_and_no_forbidden_actions():
    script = ROOT / "scripts/phase2e_04_inspect_integrated_us_phyloseq_object.py"
    text = script.read_text(encoding="utf-8").lower()

    required_outputs = [
        "phase2e_integrated_us_phyloseq_object_inventory.tsv",
        "phase2e_integrated_us_phyloseq_sample_metadata_inventory.tsv",
        "phase2e_integrated_us_phyloseq_taxa_inventory.tsv",
        "phase2e_integrated_us_metadata_column_inventory.tsv",
        "phase2e_integrated_us_extraction_plan.tsv",
        "phase2e_integrated_us_d5_status.txt",
    ]
    for out in required_outputs:
        assert out in text

    forbidden = [
        "prefetch(",
        "fasterq-dump",
        "wget ",
        "curl ",
        "kraken",
        "humann",
        "bowtie2",
        "blast",
        "diamond",
        "mafft",
    ]
    for token in forbidden:
        assert token not in text
