from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return list(r.fieldnames or []), [dict(x) for x in r]


def test_phase2c_files_exist():
    required = [
        ROOT / "metadata/phase2c_processed_table_verification_matrix.tsv",
        ROOT / "metadata/phase2c_wallen_processed_table_inventory.tsv",
        ROOT / "metadata/phase2c_replication_cohort_manual_checklist.tsv",
        ROOT / "metadata/phase2c_status_vocabulary.tsv",
        ROOT / "docs/PHASE2C_MANUAL_PROCESSED_TABLE_VERIFICATION.md",
        ROOT / "docs/PHASE2C_COHORT_DECISION_RULES.md",
        ROOT / "docs/PHASE2C_METADATA_PROVENANCE_NOTES.md",
        ROOT / "scripts/phase2c_00_validate_verification_matrix.py",
        ROOT / "tests/test_phase2c_verification_matrix.py",
    ]
    missing = [str(p) for p in required if not p.exists()]
    assert not missing, f"Missing Phase 2C files: {missing}"


def test_matrix_columns_and_wallen_and_replication_count():
    cols, rows = _read_tsv(ROOT / "metadata/phase2c_processed_table_verification_matrix.tsv")
    for c in [
        "dataset_id",
        "cohort_role",
        "sequencing_type",
        "processed_species_table_status",
        "sample_metadata_status",
        "pd_control_label_status",
        "run_accession_mapping_status",
        "ready_for_phase2d_processed_analysis",
        "ready_for_phase3_raw_read_subset_selection",
    ]:
        assert c in cols
    ids = {r["dataset_id"] for r in rows}
    assert "WALLEN_PRJNA834801" in ids
    reps = [r for r in rows if "replication_candidate" in r["cohort_role"]]
    assert len(reps) >= 3


def test_no_16s_in_active_matrix_and_all_shotgun():
    _, rows = _read_tsv(ROOT / "metadata/phase2c_processed_table_verification_matrix.tsv")
    ids = {r["dataset_id"] for r in rows}
    assert "SIXTEEN_S_ONLY_PD_COHORTS" not in ids
    for r in rows:
        assert r["sequencing_type"] == "shotgun_metagenomics"


def test_status_vocab_exists_and_readiness_rules_not_violated():
    _, vocab_rows = _read_tsv(ROOT / "metadata/phase2c_status_vocabulary.tsv")
    vocab = {r["status"] for r in vocab_rows}
    assert "yes" in vocab
    assert "wallen_processed_tables_verified_phase2d_ready" in vocab
    _, rows = _read_tsv(ROOT / "metadata/phase2c_processed_table_verification_matrix.tsv")
    status_fields = [
        "processed_species_table_status",
        "processed_genus_table_status",
        "processed_gene_family_table_status",
        "processed_pathway_table_status",
        "sample_metadata_status",
        "pd_control_label_status",
        "clinical_metadata_status",
        "vagal_body_first_metadata_status",
        "run_accession_mapping_status",
        "manual_verification_status",
        "ready_for_phase2d_processed_analysis",
        "ready_for_phase3_raw_read_subset_selection",
    ]
    for r in rows:
        for f in status_fields:
            assert r[f] in vocab
        if r["ready_for_phase2d_processed_analysis"] == "yes":
            assert r["processed_species_table_status"] == "available"
            assert r["sample_metadata_status"] == "available"
            assert r["pd_control_label_status"] == "available"
        if r["ready_for_phase3_raw_read_subset_selection"] == "yes":
            assert r["run_accession_mapping_status"] == "available"
            assert r["sample_metadata_status"] == "available"
            assert r["pd_control_label_status"] == "available"


def test_wallen_phase2d_yes_phase3_pending_allowed_with_missing_run_mapping():
    _, rows = _read_tsv(ROOT / "metadata/phase2c_processed_table_verification_matrix.tsv")
    wallen = [r for r in rows if r["dataset_id"] == "WALLEN_PRJNA834801"]
    assert len(wallen) == 1
    r = wallen[0]
    assert r["manual_verification_status"] == "wallen_processed_tables_verified_phase2d_ready"
    assert r["ready_for_phase2d_processed_analysis"] == "yes"
    assert r["ready_for_phase3_raw_read_subset_selection"] in {"pending", "no"}
    assert r["run_accession_mapping_status"] == "not_available"
