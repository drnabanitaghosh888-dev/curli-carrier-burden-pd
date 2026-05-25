from __future__ import annotations

import py_compile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_component_inspection_outputs_declared_and_script_compiles():
    script = ROOT / "scripts/phase2e_04_inspect_integrated_us_phyloseq_object.py"
    assert script.exists()
    py_compile.compile(str(script), doraise=True)
    text = script.read_text(encoding="utf-8").lower()

    required_outputs = [
        "phase2e_integrated_us_phyloseq_component_inventory.tsv",
        "phase2e_integrated_us_candidate_component_selection.tsv",
        "phase2e_integrated_us_label_mapping_candidates.tsv",
        "phase2e_integrated_us_taxonomy_rank_readiness.tsv",
        "phase2e_integrated_us_d5b_status.txt",
        "phase2e_integrated_us_d5b_report.md",
    ]
    for out in required_outputs:
        assert out in text


def test_no_forbidden_actions_and_ready_guard_present():
    script = ROOT / "scripts/phase2e_04_inspect_integrated_us_phyloseq_object.py"
    text = script.read_text(encoding="utf-8").lower()

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

    # Guard for list-only unresolved state should be explicit
    assert "partial_ready_needs_list_component_inspection" in text
