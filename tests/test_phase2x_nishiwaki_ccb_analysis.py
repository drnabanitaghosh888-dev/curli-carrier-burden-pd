from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/phase2x_01_nishiwaki_ccb_analysis.py"


def test_phase2x_script_exists_and_targets_nishiwaki():
    assert SCRIPT.exists()
    text = SCRIPT.read_text(encoding="utf-8")
    assert "NishiwakiH_2024" in text
    assert "Nishiwaki" in text
    assert "phase2x_nishiwaki_ccb_status.txt" in text
    assert "phase2x_nishiwaki_final_decision_report.txt" in text


def test_phase2x_contains_ccb_formula_and_outputs():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "CCB_s = sum_i abundance_is * weight_i" in text
    assert "log1p_CCB" in text
    assert "phase2x_nishiwaki_sample_ccb.tsv" in text
    assert "phase2x_nishiwaki_matched_curli_taxa.tsv" in text
    assert "phase2x_nishiwaki_ccb_test_summary.tsv" in text


def test_phase2x_uses_processed_relative_abundance_uuid_resource():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "relative_abundance_uuid.parquet" in text
    assert "https://huggingface.co/datasets/waldronlab/metagenomics_mac/resolve/main/" in text
    assert "duckdb" in text
    assert "LOAD httpfs" in text


def test_phase2x_has_no_raw_read_or_external_command_execution():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    assert "subprocess" not in text
    assert "requests" not in text
    assert "urllib" not in text
    blocked_terms = [
        "prefetch",
        "fasterq",
        "fastq",
        " sra",
        "bowtie2",
        "hmmer",
        "blast",
        "diamond",
        "mafft",
        "metaphlan",
        "humann",
    ]
    for term in blocked_terms:
        assert term not in text


def test_phase2x_analysis_script_does_not_use_process_execution():
    script_text = Path("scripts/phase2x_01_nishiwaki_ccb_analysis.py").read_text(encoding="utf-8").lower()
    blocked_patterns = [
        r"\bsubprocess\b",
        r"\bos\.system\b",
        r"\bpytest\.main\b",
        r"\bpopen\b",
        r"\bprefetch\b",
        r"\bfasterq\b",
        r"fastq-dump",
        r"\bbowtie2\b",
        r"\bhmmer\b",
        r"\bblast\b",
        r"\bdiamond\b",
        r"\bmafft\b",
        r"\bmetaphlan\b",
        r"\bhumann\b",
        r"\bwget\b",
        r"\bcurl\b",
        r"requests\.",
        r"urllib",
    ]
    for pattern in blocked_patterns:
        assert re.search(pattern, script_text) is None
