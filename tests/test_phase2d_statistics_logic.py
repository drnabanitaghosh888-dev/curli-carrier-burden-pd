from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/phase2d_04_run_wallen_discovery_statistics.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("phase2d_04", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_log1p_and_presence_features_created():
    mod = _load_module()
    df = pd.DataFrame({"CurliCarrierBurden": [0.0, 0.1, 1.0]})
    out = mod._add_burden_features(df)
    assert "CurliCarrierBurden_log1p" in out.columns
    assert "CurliCarrierPresence" in out.columns
    assert np.isclose(out.loc[0, "CurliCarrierBurden_log1p"], 0.0)
    assert out["CurliCarrierPresence"].tolist() == [0, 1, 1]


def test_script_uses_log1p_and_presence_models_and_execute_gate():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "CurliCarrierBurden_log1p" in text
    assert "CurliCarrierPresence" in text
    assert "--execute" in text
    assert "if not args.execute" in text
    assert "out_runtime_plan" in text


def test_no_fastq_sra_or_raw_read_commands_present():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    forbidden = [
        "prefetch",
        "fasterq-dump",
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
    for token in forbidden:
        assert token not in text

