from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return list(r.fieldnames or []), [dict(x) for x in r]


def test_phase2b_files_exist():
    required = [
        ROOT / "config/phase2b_config.yml",
        ROOT / "metadata/phase2b_processed_table_sources.tsv",
        ROOT / "metadata/phase2b_required_table_schema.tsv",
        ROOT / "metadata/phase2b_cohort_verification_queue.tsv",
        ROOT / "metadata/phase2b_excluded_or_context_only_datasets.tsv",
        ROOT / "scripts/phase2b_00_validate_processed_table_plan.py",
        ROOT / "scripts/phase2b_01_build_curli_candidate_taxa.py",
        ROOT / "scripts/phase2b_02_check_local_processed_tables.py",
        ROOT / "scripts/phase2b_03_prepare_curli_carrier_burden_scaffold.py",
        ROOT / "scripts/manual_downloads/download_wallen_processed_tables_manual.sh",
        ROOT / "scripts/manual_downloads/download_replication_processed_tables_manual.sh",
        ROOT / "docs/PHASE2B_PROCESSED_TABLE_PLAN.md",
        ROOT / "docs/PHASE2B_MULTICOOHORT_STRATEGY.md",
        ROOT / "docs/PHASE2B_TABLE_SCHEMA.md",
        ROOT / "tests/test_phase2b_processed_table_plan.py",
    ]
    missing = [str(p) for p in required if not p.exists()]
    assert not missing, f"Missing Phase 2B files: {missing}"


def test_sources_columns_and_wallen_priority():
    cols, rows = _read_tsv(ROOT / "metadata/phase2b_processed_table_sources.tsv")
    required_cols = [
        "dataset_id",
        "cohort_role",
        "sequencing_type",
        "priority",
        "manual_verification_status",
        "download_status",
        "notes",
    ]
    for c in required_cols:
        assert c in cols
    by_id = {r["dataset_id"]: r for r in rows}
    assert "WALLEN_PRJNA834801" in by_id
    assert by_id["WALLEN_PRJNA834801"]["priority"] == "highest"


def test_replication_candidates_count():
    _, rows = _read_tsv(ROOT / "metadata/phase2b_processed_table_sources.tsv")
    reps = [r for r in rows if "replication_candidate" in r["cohort_role"]]
    assert len(reps) >= 2


def test_active_registry_contains_no_16s_rows():
    _, rows = _read_tsv(ROOT / "metadata/phase2b_processed_table_sources.tsv")
    has_16s = any("16s" in (r["sequencing_type"] or "").lower() for r in rows)
    assert not has_16s, "Active Phase 2B registry must not contain 16S rows."


def test_16s_dataset_only_in_excluded_table_and_forbidden_from_core_use():
    _, rows_active = _read_tsv(ROOT / "metadata/phase2b_processed_table_sources.tsv")
    active_ids = {r["dataset_id"] for r in rows_active}
    assert "SIXTEEN_S_ONLY_PD_COHORTS" not in active_ids

    cols_ex, rows_ex = _read_tsv(ROOT / "metadata/phase2b_excluded_or_context_only_datasets.tsv")
    assert "dataset_id" in cols_ex
    by_id = {r["dataset_id"]: r for r in rows_ex}
    assert "SIXTEEN_S_ONLY_PD_COHORTS" in by_id
    r = by_id["SIXTEEN_S_ONLY_PD_COHORTS"]
    assert r["study_type"] == "16S_amplicon"
    forbid = r["not_allowed_use"]
    for token in [
        "curli_gene_detection",
        "csg_operon_analysis",
        "raw_csg_screening",
        "strain_reconstruction",
        "primary_association_testing",
    ]:
        assert token in forbid


def test_manual_scripts_guarded_and_no_unguarded_heavy_commands():
    for p in [
        ROOT / "scripts/manual_downloads/download_wallen_processed_tables_manual.sh",
        ROOT / "scripts/manual_downloads/download_replication_processed_tables_manual.sh",
    ]:
        text = p.read_text(encoding="utf-8")
        assert 'RUN_MANUAL_DOWNLOAD' in text
        assert 'Manual download guard active' in text
        # ensure no active wget/curl/prefetch/fasterq commands (commented is acceptable)
        active_lines = [
            ln.strip()
            for ln in text.splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
        forbidden = ["wget ", "curl ", "prefetch ", "fasterq-dump ", "kraken", "metaphlan", "humann", "bowtie2"]
        for line in active_lines:
            assert not any(tok in line.lower() for tok in forbidden), f"forbidden active token in {p}: {line}"


def test_phase2b_reports_dir_exists_or_is_creatable():
    report_dir = ROOT / "data/phase2b/reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    assert report_dir.exists()
