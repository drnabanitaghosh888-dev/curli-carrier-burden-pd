from __future__ import annotations

import csv
import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return list(r.fieldnames or []), [dict(x) for x in r]


def test_phase2e_files_exist():
    required = [
        ROOT / "metadata/phase2e_replication_cohort_verification_matrix.tsv",
        ROOT / "metadata/phase2e_replication_expected_files.tsv",
        ROOT / "metadata/phase2e_replication_readiness_rules.tsv",
        ROOT / "metadata/phase2e_replication_manual_source_tracking.tsv",
        ROOT / "metadata/phase2e_dataset_id_aliases.tsv",
        ROOT / "docs/PHASE2E_REPLICATION_PROCESSED_TABLE_VERIFICATION_PLAN.md",
        ROOT / "docs/PHASE2E_REPLICATION_DECISION_RULES.md",
        ROOT / "docs/PHASE2E_INTERPRETATION_AFTER_WALLEN_DISCOVERY.md",
        ROOT / "scripts/phase2e_00_validate_replication_verification_matrix.py",
        ROOT / "scripts/phase2e_01_check_local_replication_processed_tables.py",
        ROOT / "tests/test_phase2e_replication_scaffold.py",
    ]
    missing = [str(p) for p in required if not p.exists()]
    assert not missing, f"Missing Phase 2E scaffold files: {missing}"


def test_matrix_columns_and_cohorts_and_16s_exclusion():
    cols, rows = _read_tsv(ROOT / "metadata/phase2e_replication_cohort_verification_matrix.tsv")
    for c in [
        "dataset_id",
        "sequencing_type",
        "processed_species_table_status",
        "sample_metadata_status",
        "pd_control_label_status",
        "ready_for_replication_processed_analysis",
        "ready_for_raw_read_validation",
    ]:
        assert c in cols

    ids = {r["dataset_id"] for r in rows}
    for cohort in [
        "PALACIOS_PRODROMAL_PD",
        "NISHIWAKI_MULTICOUNTRY_PD",
        "MAO_CENTRAL_CHINA_PD",
        "INTEGRATED_US_MULTICOHORT_PD",
    ]:
        assert cohort in ids
    assert "INTEGRATED_US_MULTICOOHORT_PD" not in ids
    assert "SIXTEEN_S_ONLY_PD_COHORTS" not in ids

    for r in rows:
        assert r["sequencing_type"] == "shotgun_metagenomics"


def test_expected_paths_keep_multicohort_slug():
    _, rows = _read_tsv(ROOT / "metadata/phase2e_replication_expected_files.tsv")
    us_rows = [r for r in rows if r["dataset_id"] == "INTEGRATED_US_MULTICOHORT_PD"]
    assert us_rows, "Missing canonical Integrated-US rows"
    for r in us_rows:
        assert "integrated_us_multicohort_pd/" in r["local_expected_path"]
        assert "integrated_us_multicoohort_pd/" not in r["local_expected_path"]


def test_scripts_compile_and_no_forbidden_commands():
    scripts = [
        ROOT / "scripts/phase2e_00_validate_replication_verification_matrix.py",
        ROOT / "scripts/phase2e_01_check_local_replication_processed_tables.py",
    ]
    for s in scripts:
        py_compile.compile(str(s), doraise=True)

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
    for s in scripts:
        text = s.read_text(encoding="utf-8").lower()
        for token in forbidden:
            assert token not in text, f"Forbidden token '{token}' in {s}"
