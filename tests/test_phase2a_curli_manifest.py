from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]

EXPECTED_GENES = ["csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"]
VALID_WEIGHTS = {"high", "moderate", "low"}
REQUIRED_COLUMNS = [
    "gene_symbol",
    "functional_role",
    "biological_weight",
    "core_score",
    "operon_completeness",
]


def _load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        obj = yaml.safe_load(handle)
    assert isinstance(obj, dict)
    return obj


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), [dict(r) for r in reader]


def test_phase2a_files_exist():
    required = [
        ROOT / "config/curli_genes.yml",
        ROOT / "config/phase2_reference_config.yml",
        ROOT / "data/reference_manifest/curli_gene_manifest.tsv",
        ROOT / "scripts/build_curli_reference_manifest.py",
        ROOT / "scripts/validate_curli_reference_manifest.py",
        ROOT / "tests/test_phase2a_curli_manifest.py",
    ]
    missing = [str(p) for p in required if not p.exists()]
    assert not missing, f"Missing Phase 2A files: {missing}"


def test_manifest_schema_and_genes():
    manifest = ROOT / "data/reference_manifest/curli_gene_manifest.tsv"
    cols, rows = _read_tsv(manifest)
    missing_cols = [c for c in REQUIRED_COLUMNS if c not in cols]
    assert not missing_cols, f"Manifest missing columns: {missing_cols}"
    assert len(rows) == 7, "Manifest must contain exactly seven curli genes"
    genes = [r["gene_symbol"] for r in rows]
    assert set(genes) == set(EXPECTED_GENES), f"Unexpected genes in manifest: {genes}"


def test_manifest_values_valid():
    manifest = ROOT / "data/reference_manifest/curli_gene_manifest.tsv"
    _, rows = _read_tsv(manifest)
    bool_like = {"yes", "no", "true", "false", "1", "0"}
    for row in rows:
        assert row["functional_role"].strip(), f"Empty role for {row['gene_symbol']}"
        assert row["biological_weight"].strip().lower() in VALID_WEIGHTS
        assert row["core_score"].strip().lower() in bool_like
        assert row["operon_completeness"].strip().lower() in bool_like


def test_configs_define_expected_gene_set():
    genes_cfg = _load_yaml(ROOT / "config/curli_genes.yml")
    phase2_cfg = _load_yaml(ROOT / "config/phase2_reference_config.yml")
    genes = list((genes_cfg.get("curli_genes") or {}).keys())
    expected = list((phase2_cfg.get("phase2_reference") or {}).get("expected_genes", []))
    assert set(genes) == set(EXPECTED_GENES)
    assert set(expected) == set(EXPECTED_GENES)


def test_phase2_config_has_safety_defaults():
    cfg = _load_yaml(ROOT / "config/phase2_reference_config.yml")
    safety = (cfg.get("phase2_reference") or {}).get("safety", {})
    assert safety.get("default_dry_run") is True
    assert safety.get("require_explicit_heavy_flag") is True
    assert safety.get("heavy_execution_flag_name") == "--allow-heavy"


def test_scripts_dry_run_mode_is_lightweight():
    build = subprocess.run(
        [
            sys.executable,
            "scripts/build_curli_reference_manifest.py",
            "--dry-run",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "[DRY-RUN]" in build.stdout

    validate = subprocess.run(
        [
            sys.executable,
            "scripts/validate_curli_reference_manifest.py",
            "--dry-run",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "[DRY-RUN]" in validate.stdout
    assert "heavy tools blocked" in validate.stdout

