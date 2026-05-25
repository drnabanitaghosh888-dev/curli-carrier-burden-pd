from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "data/reference_manifest/curli_reference_plan.tsv"
EXPECTED_GENES = {"csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"}
REQUIRED_COLUMNS = {
    "gene_symbol",
    "reference_group",
    "target_taxa",
    "intended_reference_source",
    "expected_output_fasta",
    "expected_output_alignment",
    "expected_output_hmm",
    "requires_download",
    "requires_alignment",
    "requires_hmmbuild",
    "heavy_step",
}
BOOL_LIKE = {"yes", "no", "true", "false", "1", "0"}


def _read_plan() -> tuple[list[str], list[dict[str, str]]]:
    with PLAN.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return list(r.fieldnames or []), [dict(x) for x in r]


def _inside_project(p: Path) -> bool:
    root = ROOT.resolve()
    target = p.resolve()
    return target == root or root in target.parents


def test_phase2b_files_exist():
    required = [
        ROOT / "config/phase2_reference_config.yml",
        ROOT / "data/reference_manifest/curli_reference_plan.tsv",
        ROOT / "scripts/plan_curli_reference_build.py",
        ROOT / "scripts/validate_curli_reference_plan.py",
        ROOT / "tests/test_phase2b_reference_plan.py",
    ]
    missing = [str(p) for p in required if not p.exists()]
    assert not missing, f"Missing files: {missing}"


def test_plan_has_required_columns_and_all_genes():
    cols, rows = _read_plan()
    missing_cols = sorted(REQUIRED_COLUMNS - set(cols))
    assert not missing_cols, f"Missing columns: {missing_cols}"
    genes = {row["gene_symbol"] for row in rows}
    assert genes == EXPECTED_GENES


def test_path_fields_stay_inside_project_and_bool_fields_valid():
    _, rows = _read_plan()
    for row in rows:
        for field in ["expected_output_fasta", "expected_output_alignment", "expected_output_hmm"]:
            p = ROOT / row[field]
            assert _inside_project(p), f"Path escapes project: {p}"
        for field in ["heavy_step", "requires_download", "requires_alignment", "requires_hmmbuild"]:
            assert row[field].strip().lower() in BOOL_LIKE


def test_dry_run_does_not_create_heavy_outputs_and_reports_actions():
    cp = subprocess.run(
        [sys.executable, "scripts/plan_curli_reference_build.py", "--dry-run"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "[DRY-RUN]" in cp.stdout
    assert "heavy steps are declared only" in cp.stdout

    cp2 = subprocess.run(
        [sys.executable, "scripts/validate_curli_reference_plan.py", "--dry-run"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "[DRY-RUN]" in cp2.stdout
    assert "no heavy execution" in cp2.stdout

    _, rows = _read_plan()
    # Phase 2B should not have produced heavy outputs.
    for row in rows:
        for field in ["expected_output_fasta", "expected_output_alignment", "expected_output_hmm"]:
            out = ROOT / row[field]
            assert not out.exists(), f"Heavy output exists unexpectedly during Phase 2B: {out}"


def test_plan_generation_is_deterministic():
    # Regenerate in a temp path twice and compare bytes.
    out1 = ROOT / "data/reference_manifest/.tmp_plan_1.tsv"
    out2 = ROOT / "data/reference_manifest/.tmp_plan_2.tsv"
    try:
        subprocess.run(
            [sys.executable, "scripts/plan_curli_reference_build.py", "--output", str(out1.relative_to(ROOT))],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        subprocess.run(
            [sys.executable, "scripts/plan_curli_reference_build.py", "--output", str(out2.relative_to(ROOT))],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        b1 = out1.read_bytes()
        b2 = out2.read_bytes()
        assert b1 == b2, "Plan generation is not deterministic"
    finally:
        if out1.exists():
            out1.unlink()
        if out2.exists():
            out2.unlink()


def test_config_declares_expected_reference_scope():
    cfg = yaml.safe_load((ROOT / "config/phase2_reference_config.yml").read_text())
    p2 = cfg["phase2_reference"]
    expected = set(p2["expected_genes"])
    assert expected == EXPECTED_GENES
    target_taxa = p2["reference_plan"]["target_taxa"]
    assert "Escherichia coli" in target_taxa
    assert "Salmonella enterica" in target_taxa
    assert "Enterobacteriaceae" in target_taxa

