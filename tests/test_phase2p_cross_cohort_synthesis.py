from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/phase2p_01_cross_cohort_evidence_synthesis.py"
OUT_DIR = ROOT / "data/phase2p/reports"
TABLE = OUT_DIR / "phase2p_cross_cohort_evidence_table.tsv"
SUMMARY = OUT_DIR / "phase2p_cross_cohort_evidence_summary.txt"
STATUS = OUT_DIR / "phase2p_cross_cohort_evidence_status.txt"


def test_phase2p_cross_cohort_synthesis_smoke():
    assert SCRIPT.exists()
    subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)

    assert TABLE.exists()
    assert SUMMARY.exists()
    assert STATUS.exists()

    table = pd.read_csv(TABLE, sep="\t")
    romano = table[table["dataset_id"] == "ROMANO_NONWALLEN"]
    assert not romano.empty
    assert (
        romano.iloc[0]["decision"]
        == "WEAK_COHORT_SENSITIVE_SUPPORT_NOT_DEFINITIVE_REPLICATION"
    )

    summary = SUMMARY.read_text(encoding="utf-8")
    assert "cohort-sensitive" in summary
