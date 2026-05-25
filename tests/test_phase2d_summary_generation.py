from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/phase2d_05_generate_wallen_discovery_summary.py"


def test_summary_script_exists_and_uses_expected_inputs_outputs():
    assert SCRIPT.exists()
    text = SCRIPT.read_text(encoding="utf-8")
    assert "wallen_analysis_metadata.tsv" in text
    assert "wallen_curli_candidate_species_matches.tsv" in text
    assert "wallen_curli_burden_species_set.tsv" in text
    assert "wallen_curli_carrier_burden.tsv" in text
    assert "wallen_curli_burden_statistics.tsv" in text
    assert "PHASE2D_WALLEN_DISCOVERY_SUMMARY.md" in text
    assert "phase2d_wallen_discovery_file_manifest.tsv" in text
    assert "phase2d_final_status.txt" in text


def test_summary_contains_required_scientific_language_and_nonfabrication():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "processed-table proxy only" in text
    assert "attenuate after age, sex, and BMI adjustment" in text
    assert "does not fabricate missing results" in text
    assert "raw-read csg validation are required" in text
    assert "run-accession mapping becomes available" in text


def test_no_phase3_claims_and_no_forbidden_commands():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    assert "phase 3 raw-read validation" in text
    forbidden = [
        "prefetch",
        "fasterq-dump",
        "fastq",
        "hmmsearch",
        "hmmbuild",
        "kraken",
        "humann",
        "metaphlan ",
        "bowtie2",
        "blast",
        "diamond",
        "mafft",
    ]
    for token in forbidden:
        assert token not in text
