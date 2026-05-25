from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return list(r.fieldnames or []), [dict(x) for x in r]


def test_phase2c_w1_files_exist():
    required = [
        ROOT / "metadata/phase2c_wallen_manual_source_verification.tsv",
        ROOT / "metadata/phase2c_wallen_file_renaming_plan.tsv",
        ROOT / "docs/PHASE2C_WALLEN_MANUAL_VERIFICATION_RUNBOOK.md",
        ROOT / "scripts/phase2c_01_check_wallen_local_files.py",
        ROOT / "tests/test_phase2c_wallen_verification_scaffold.py",
    ]
    missing = [str(p) for p in required if not p.exists()]
    assert not missing, f"Missing Phase 2C-W1 scaffold files: {missing}"


def test_manual_source_verification_schema_and_defaults():
    cols, rows = _read_tsv(ROOT / "metadata/phase2c_wallen_manual_source_verification.tsv")
    req_cols = [
        "item_id",
        "resource_needed",
        "expected_processed_or_metadata_file",
        "source_to_check",
        "verified_source_url",
        "verified_file_name_at_source",
        "local_target_path",
        "verification_status",
        "required_for_phase2d",
        "required_for_phase3",
        "notes",
    ]
    for c in req_cols:
        assert c in cols
    assert len(rows) >= 10
    ids = {r["item_id"] for r in rows}
    expected_ids = {
        "sample_metadata",
        "species_abundance",
        "genus_abundance",
        "gene_family_abundance",
        "pathway_abundance",
        "run_accession_mapping",
        "clinical_metadata",
        "raw_accession_metadata",
        "README_or_data_dictionary",
        "license_or_access_terms",
    }
    assert expected_ids.issubset(ids)
    for r in rows:
        assert r["verification_status"] == "pending_manual_verification"


def test_renaming_plan_has_standardized_targets():
    _, rows = _read_tsv(ROOT / "metadata/phase2c_wallen_file_renaming_plan.tsv")
    targets = {r["standardized_local_file_name"] for r in rows}
    required = {
        "sample_metadata.tsv",
        "species_abundance.tsv",
        "genus_abundance.tsv",
        "gene_family_abundance.tsv",
        "pathway_abundance.tsv",
        "run_accession_mapping.tsv",
    }
    assert required.issubset(targets)

