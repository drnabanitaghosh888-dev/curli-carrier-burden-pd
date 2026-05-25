from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/phase2u_01_duru_ccb_analysis.py"


def test_phase2u_script_exists_and_declares_expected_outputs():
    assert SCRIPT.exists()
    text = SCRIPT.read_text(encoding="utf-8")
    assert "phase2u_duru_ccb_status.txt" in text
    assert "phase2u_duru_final_decision_report.txt" in text
    assert "phase2u_duru_ccb_sample_burden.tsv" in text
    assert "phase2u_duru_ccb_matched_taxa.tsv" in text
    assert "phase2u_duru_ccb_test_summary.tsv" in text


def test_phase2u_uses_local_processed_duru_inputs_only():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "duruIC_2024_uuid_phenotype.tsv" in text
    assert "duruIC_2024_relative_abundance_species_like.tsv" in text
    assert "metadata/curli_candidate_taxa_from_phase1.tsv" in text
    assert "requests" not in text
    assert "urllib" not in text
    assert "subprocess" not in text


def test_phase2u_contains_ccb_logic_without_sequence_workflows():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "CurliCarrierBurden" in text
    assert "CurliCarrierBurden_log1p" in text
    assert "CurliCarrierPresence" in text
    assert "mannwhitneyu" in text
    assert "cliffs_delta" in text
    blocked_terms = [
        "HMMER",
        "BLAST",
        "DIAMOND",
        "MAFFT",
        "MetaPhlAn",
        "HUMAnN",
        "fasterq",
        "prefetch",
    ]
    for term in blocked_terms:
        assert term not in text


def test_phase2u_analysis_script_does_not_use_process_execution():
    script_text = Path("scripts/phase2u_01_duru_ccb_analysis.py").read_text(encoding="utf-8")
    forbidden_patterns = [
        "sub" + "process",
        "os.system",
        "pytest.main",
        "Popen",
        "run(",
    ]
    for pattern in forbidden_patterns:
        assert pattern not in script_text
