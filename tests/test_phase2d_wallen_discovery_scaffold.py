from __future__ import annotations

import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_phase2d_config_and_metadata_exist():
    required = [
        ROOT / "config/phase2d_config.yml",
        ROOT / "metadata/phase2d_wallen_analysis_variables.tsv",
        ROOT / "metadata/phase2d_curli_taxon_matching_rules.tsv",
        ROOT / "metadata/phase2d_model_plan.tsv",
    ]
    missing = [str(p) for p in required if not p.exists()]
    assert not missing, f"Missing Phase 2D config/metadata files: {missing}"


def test_phase2d_docs_exist_and_state_proxy_limits():
    docs = [
        ROOT / "docs/PHASE2D_WALLEN_DISCOVERY_ANALYSIS_PLAN.md",
        ROOT / "docs/PHASE2D_CURLI_CARRIER_BURDEN_INTERPRETATION.md",
        ROOT / "docs/PHASE2D_STATISTICAL_MODEL_PLAN.md",
    ]
    missing = [str(p) for p in docs if not p.exists()]
    assert not missing, f"Missing Phase 2D docs: {missing}"
    text = (ROOT / "docs/PHASE2D_CURLI_CARRIER_BURDEN_INTERPRETATION.md").read_text(encoding="utf-8").lower()
    assert "processed-table proxy" in text
    assert "not proof of csg gene presence" in text
    assert "phase 3 remains blocked" in text


def test_phase2d_scripts_exist_and_compile():
    scripts = [
        ROOT / "scripts/phase2d_00_validate_wallen_inputs.py",
        ROOT / "scripts/phase2d_01_prepare_wallen_analysis_metadata.py",
        ROOT / "scripts/phase2d_02_match_curli_candidate_taxa_to_wallen_species.py",
        ROOT / "scripts/phase2d_03_compute_curli_carrier_burden.py",
        ROOT / "scripts/phase2d_04_run_wallen_discovery_statistics.py",
        ROOT / "scripts/phase2d_05_generate_wallen_discovery_summary.py",
    ]
    missing = [str(p) for p in scripts if not p.exists()]
    assert not missing, f"Missing Phase 2D scripts: {missing}"
    for s in scripts:
        py_compile.compile(str(s), doraise=True)


def test_execute_guards_and_no_forbidden_download_or_heavy_tool_calls():
    guarded = [
        ROOT / "scripts/phase2d_01_prepare_wallen_analysis_metadata.py",
        ROOT / "scripts/phase2d_02_match_curli_candidate_taxa_to_wallen_species.py",
        ROOT / "scripts/phase2d_03_compute_curli_carrier_burden.py",
        ROOT / "scripts/phase2d_04_run_wallen_discovery_statistics.py",
    ]
    for p in guarded:
        text = p.read_text(encoding="utf-8")
        assert "--execute" in text
        assert "if not args.execute" in text

    forbidden = [
        "prefetch",
        "fasterq-dump",
        "sra",
        "fastq",
        "hmmsearch",
        "hmmbuild",
        "kraken",
        "humann",
        "bowtie2",
        "blast",
        "diamond",
        "mafft",
    ]
    for p in guarded + [ROOT / "scripts/phase2d_00_validate_wallen_inputs.py", ROOT / "scripts/phase2d_05_generate_wallen_discovery_summary.py"]:
        txt = p.read_text(encoding="utf-8").lower()
        for token in forbidden:
            assert token not in txt, f"Forbidden token '{token}' found in {p}"


def test_phase2d_does_not_require_run_accession_mapping():
    cfg = (ROOT / "config/phase2d_config.yml").read_text(encoding="utf-8")
    assert "input_files:" in cfg
    assert "run_accession_mapping:" not in cfg
    validator = (ROOT / "scripts/phase2d_00_validate_wallen_inputs.py").read_text(encoding="utf-8").lower()
    assert "not required for phase 2d" in validator
