import subprocess
import sys
from pathlib import Path
import re

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/phase2w_01_asapmac_pd_cohort_screen.py"
REPORT_DIR = ROOT / "data/phase2w/reports"
OUT_TABLE = REPORT_DIR / "phase2w_asapmac_eligible_pd_cohorts.tsv"
OUT_STATUS = REPORT_DIR / "phase2w_asapmac_pd_cohort_screen_status.txt"
OUT_SUMMARY = REPORT_DIR / "phase2w_asapmac_pd_cohort_screen_summary.txt"


def test_phase2w_asapmac_pd_cohort_screen_outputs():
    assert SCRIPT.exists()
    subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)

    assert OUT_TABLE.exists()
    assert OUT_STATUS.exists()
    assert OUT_SUMMARY.exists()
    assert OUT_STATUS.read_text(encoding="utf-8").strip() == "PASS"

    table = pd.read_csv(OUT_TABLE, sep="\t")
    assert not table.empty
    assert (table["n_PD"].astype(int) >= 20).all()
    assert (table["n_HC"].astype(int) >= 20).all()

    summary = OUT_SUMMARY.read_text(encoding="utf-8")
    assert "Nishiwaki" in summary
    assert "ccb_computed=no" in summary


def test_phase2w_script_has_no_raw_read_commands():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    blocked_patterns = [
        r"\\bprefetch\\b",
        r"\\bfasterq\\b",
        r"fastq-dump",
        r"\\bbowtie2\\b",
        r"\\bhmmer\\b",
        r"\\bblast\\b",
        r"\\bdiamond\\b",
        r"\\bmafft\\b",
        r"\\bmetaphlan\\b",
        r"\\bhumann\\b",
        r"\\bwget\\b",
        r"\\bcurl\\b",
        r"requests\\.",
        r"urllib",
    ]
    for pattern in blocked_patterns:
        assert re.search(pattern, text) is None
