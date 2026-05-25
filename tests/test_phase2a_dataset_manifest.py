from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "metadata/pd_metagenome_sources.tsv"
PRIORITY = ROOT / "metadata/phase2a_dataset_priority.tsv"

REQUIRED_COLUMNS = [
    "dataset_id",
    "study_name",
    "disease",
    "sample_type",
    "sequencing_type",
    "data_level_available",
    "accession_or_source",
    "raw_data_available",
    "processed_tables_available",
    "clinical_metadata_available",
    "vagal_body_first_metadata_available",
    "estimated_sample_count",
    "geography",
    "access_status",
    "priority",
    "intended_phase2_use",
    "download_status",
    "notes",
]

VALID_DOWNLOAD_STATUS = {
    "not_started",
    "manual_required",
    "downloaded",
    "unavailable",
    "candidate_pending_manual_verification",
}


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), [dict(r) for r in reader]


def test_phase2a_files_exist():
    required = [
        ROOT / "metadata/pd_metagenome_sources.tsv",
        ROOT / "metadata/phase2a_dataset_priority.tsv",
        ROOT / "scripts/phase2a_00_validate_dataset_manifest.py",
        ROOT / "scripts/manual_downloads/download_wallen_prjna834801_manual.sh",
        ROOT / "scripts/manual_downloads/download_additional_pd_cohorts_manual.sh",
        ROOT / "docs/PHASE2A_DATASET_ACQUISITION_PLAN.md",
        ROOT / "tests/test_phase2a_dataset_manifest.py",
    ]
    missing = [str(p) for p in required if not p.exists()]
    assert not missing, f"Missing required Phase 2A files: {missing}"


def test_manifest_required_columns_and_primary_dataset():
    cols, rows = _read_tsv(MANIFEST)
    missing_cols = [c for c in REQUIRED_COLUMNS if c not in cols]
    assert not missing_cols, f"Missing manifest columns: {missing_cols}"
    ids = {r["dataset_id"] for r in rows}
    assert "WALLEN_PRJNA834801" in ids
    wallen = [r for r in rows if r["dataset_id"] == "WALLEN_PRJNA834801"][0]
    assert wallen["accession_or_source"] == "PRJNA834801"
    assert "shotgun" in wallen["sequencing_type"].lower()


def test_16s_not_marked_for_raw_csg_screening():
    _, rows = _read_tsv(MANIFEST)
    for r in rows:
        if "16s" in r["sequencing_type"].lower():
            use = r["intended_phase2_use"].lower()
            assert not ("raw" in use and "csg" in use), (
                f"16S dataset wrongly marked for raw csg screening: {r['dataset_id']}"
            )


def test_download_status_values_valid():
    _, rows = _read_tsv(MANIFEST)
    bad = [r["dataset_id"] for r in rows if r["download_status"] not in VALID_DOWNLOAD_STATUS]
    assert not bad, f"Invalid download_status values for: {bad}"


def test_manual_download_scripts_are_guarded():
    for p in [
        ROOT / "scripts/manual_downloads/download_wallen_prjna834801_manual.sh",
        ROOT / "scripts/manual_downloads/download_additional_pd_cohorts_manual.sh",
    ]:
        text = p.read_text(encoding="utf-8")
        assert "RUN_MANUAL_DOWNLOAD" in text
        assert "Exiting without download" in text


def test_validator_script_runs_lightweight_and_writes_reports():
    cp = subprocess.run(
        [sys.executable, "scripts/phase2a_00_validate_dataset_manifest.py"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    # PASS_WITH_WARNINGS is acceptable in this phase because candidate rows are expected.
    assert cp.returncode in {0, 1}
    assert (ROOT / "data/phase2a/reports/phase2a_manifest_validation.tsv").exists()
    assert (ROOT / "data/phase2a/reports/phase2a_status.txt").exists()

