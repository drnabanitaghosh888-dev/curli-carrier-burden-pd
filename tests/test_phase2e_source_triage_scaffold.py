from __future__ import annotations

import csv
import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return list(r.fieldnames or []), [dict(x) for x in r]


def test_triage_and_download_plan_files_exist():
    required = [
        ROOT / "metadata/phase2e_replication_source_triage.tsv",
        ROOT / "metadata/phase2e_replication_manual_download_plan.tsv",
        ROOT / "metadata/phase2e_dataset_id_aliases.tsv",
        ROOT / "docs/PHASE2E_MANUAL_SOURCE_VERIFICATION_RUNBOOK.md",
        ROOT / "scripts/phase2e_02_summarize_replication_source_triage.py",
        ROOT / "tests/test_phase2e_source_triage_scaffold.py",
    ]
    missing = [str(p) for p in required if not p.exists()]
    assert not missing, f"Missing Phase 2E-B files: {missing}"


def test_active_cohorts_present_and_16s_absent_and_not_ready():
    _, triage = _read_tsv(ROOT / "metadata/phase2e_replication_source_triage.tsv")
    ids = {r["dataset_id"] for r in triage}
    expected = {
        "PALACIOS_PRODROMAL_PD",
        "NISHIWAKI_MULTICOUNTRY_PD",
        "MAO_CENTRAL_CHINA_PD",
        "INTEGRATED_US_MULTICOHORT_PD",
    }
    assert expected.issubset(ids)
    assert "INTEGRATED_US_MULTICOOHORT_PD" not in ids
    assert "SIXTEEN_S_ONLY_PD_COHORTS" not in ids

    for r in triage:
        if r["dataset_id"] == "INTEGRATED_US_MULTICOHORT_PD":
            assert r["ready_for_manual_download"] in {
                "no",
                "completed_or_not_needed_after_download",
            }
            assert r["ready_for_processed_replication_analysis"] != "yes"
            assert r["ready_for_raw_read_validation"] != "yes"
        else:
            assert r["ready_for_manual_download"] == "no"
            assert r["ready_for_processed_replication_analysis"] == "pending"
            assert r["ready_for_raw_read_validation"] == "pending"


def test_integrated_us_post_download_status_artifacts():
    status_file = ROOT / "data/phase2e/reports/phase2e_integrated_us_source_download_status.txt"
    summary_file = (
        ROOT
        / "data/phase2e/manual_sources/integrated_us_multicohort_pd/manifests/download_summary.txt"
    )
    skipped_marker = (
        ROOT
        / "data/phase2e/manual_sources/integrated_us_multicohort_pd/zenodo/PD_Metagenomic_Analysis_2b2d81f.zip.partial_skipped"
    )
    final_zip = Path(
        "/Users/krishnendu/Documents/NG_work2/curli_phase2e_integrated_us_multicohort_sources.zip"
    )

    assert status_file.exists()
    status_text = status_file.read_text(encoding="utf-8")
    assert "DOWNLOAD_STATUS = DOWNLOAD_COMPLETE_WITH_WARNINGS" in status_text

    assert summary_file.exists()
    summary_text = summary_file.read_text(encoding="utf-8").lower()
    assert "raw_fastq_sra_downloaded = no" in summary_text

    assert skipped_marker.exists()
    assert final_zip.exists()


def test_alias_file_records_dataset_id_correction():
    _, rows = _read_tsv(ROOT / "metadata/phase2e_dataset_id_aliases.tsv")
    hit = [r for r in rows if r["old_dataset_id"] == "INTEGRATED_US_MULTICOOHORT_PD"]
    assert len(hit) == 1
    assert hit[0]["canonical_dataset_id"] == "INTEGRATED_US_MULTICOHORT_PD"


def test_scripts_compile_and_no_forbidden_commands():
    script = ROOT / "scripts/phase2e_02_summarize_replication_source_triage.py"
    py_compile.compile(str(script), doraise=True)
    text = script.read_text(encoding="utf-8").lower()
    forbidden = [
        "prefetch",
        "fasterq-dump",
        "fastq",
        "hmmsearch",
        "hmmbuild",
        "kraken",
        "humann",
        "metaphlan",
        "bowtie2",
        "blast",
        "diamond",
        "mafft",
    ]
    for token in forbidden:
        assert token not in text


def test_manual_download_plan_has_no_raw_read_payloads():
    _, rows = _read_tsv(ROOT / "metadata/phase2e_replication_manual_download_plan.tsv")
    bad_tokens = ("fastq", "sra", "raw", "prefetch", "fasterq")
    for r in rows:
        blob = " ".join(
            [
                r.get("resource_name", ""),
                r.get("source_url_or_repository", ""),
                r.get("verified_file_name", ""),
                r.get("local_target_path", ""),
                r.get("notes", ""),
            ]
        ).lower()
        if "raw_accession" in blob:
            continue
        for token in bad_tokens:
            assert token not in blob
