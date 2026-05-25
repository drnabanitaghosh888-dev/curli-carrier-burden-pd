from pathlib import Path
import re

SCRIPT = Path("scripts/phase2z_01_sampson_ccb_analysis.py")


def test_phase2z_script_exists_and_targets_sampson():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "SampsonTR_2025" in text
    assert "relative_abundance_uuid.parquet" in text


def test_phase2z_script_contains_locked_ccb_logic():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "CCB" in text
    assert "weighted_abundance" in text
    assert "mannwhitneyu" in text
    assert "cliffs_delta" in text


def test_phase2z_script_uses_sampson_phenotype_mapping():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "Study Control" in text or "control" in text
    assert "Case" in text or "case" in text
    assert "Parkinson Disease" in text or "parkinson" in text.lower()


def test_phase2z_script_does_not_use_raw_read_workflows():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    blocked_patterns = [
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
    ]
    for pattern in blocked_patterns:
        assert re.search(pattern, text) is None


def test_phase2z_test_does_not_execute_analysis_script():
    script_text = Path("scripts/phase2z_01_sampson_ccb_analysis.py").read_text(encoding="utf-8").lower()
    blocked_patterns = [
        r"\bsubprocess\b",
        r"\bos\.system\b",
        r"\bpopen\b",
        r"\bpytest\.main\b",
    ]
    for pattern in blocked_patterns:
        assert re.search(pattern, script_text) is None
