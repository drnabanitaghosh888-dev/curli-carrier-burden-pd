from __future__ import annotations

import csv
import importlib.util
import py_compile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return list(r.fieldnames or []), [dict(x) for x in r]


def test_phase2e_c_files_exist():
    required = [
        ROOT / "metadata/phase2e_cohort_source_search_targets.tsv",
        ROOT / "metadata/phase2e_cohort_specific_verification_forms.tsv",
        ROOT / "metadata/phase2e_dataset_id_aliases.tsv",
        ROOT / "docs/PHASE2E_COHORT_BY_COHORT_SOURCE_VERIFICATION_RUNBOOK.md",
        ROOT / "scripts/phase2e_03_summarize_cohort_source_forms.py",
        ROOT / "tests/test_phase2e_cohort_source_forms.py",
    ]
    missing = [str(p) for p in required if not p.exists()]
    assert not missing, f"Missing Phase 2E-C files: {missing}"


def test_all_four_cohorts_present_and_no_16s():
    _, targets = _read_tsv(ROOT / "metadata/phase2e_cohort_source_search_targets.tsv")
    ids = {r["dataset_id"] for r in targets}
    expected = {
        "INTEGRATED_US_MULTICOHORT_PD",
        "NISHIWAKI_MULTICOUNTRY_PD",
        "MAO_CENTRAL_CHINA_PD",
        "PALACIOS_PRODROMAL_PD",
    }
    assert expected.issubset(ids)
    assert "INTEGRATED_US_MULTICOOHORT_PD" not in ids
    assert "SIXTEEN_S_ONLY_PD_COHORTS" not in ids

    _, forms = _read_tsv(ROOT / "metadata/phase2e_cohort_specific_verification_forms.tsv")
    form_ids = {r["dataset_id"] for r in forms}
    assert expected.issubset(form_ids)
    assert "INTEGRATED_US_MULTICOOHORT_PD" not in form_ids
    assert "SIXTEEN_S_ONLY_PD_COHORTS" not in form_ids


def test_default_statuses_pending_and_script_safe():
    _, forms = _read_tsv(ROOT / "metadata/phase2e_cohort_specific_verification_forms.tsv")
    for r in forms:
        if r["dataset_id"] == "INTEGRATED_US_MULTICOHORT_PD":
            if r["verification_item"] in {
                "publication_page",
                "supplementary_tables",
                "data_repository",
            }:
                assert r["manual_check_status"] in {
                    "verified",
                    "pending_manual_verification",
                    "unclear",
                }
            elif r["verification_item"] in {
                "species_abundance_table",
                "sample_metadata",
                "pd_control_labels",
            }:
                assert r["manual_check_status"] in {
                    "pending_manual_verification",
                    "unclear",
                }
            elif r["verification_item"] == "license_or_access_terms":
                assert r["manual_check_status"] in {
                    "verified",
                    "pending_manual_verification",
                    "unclear",
                }
            elif r["verification_item"] == "data_dictionary_or_README":
                assert r["manual_check_status"] in {
                    "verified",
                    "pending_manual_verification",
                    "unclear",
                }
            else:
                assert r["manual_check_status"] in {
                    "pending_manual_verification",
                    "unclear",
                }
            assert r["evidence_quality"] in {"unknown", "unclear", "unavailable", "moderate", "strong"}
            assert r["ready_status_after_check"] == "pending"
            assert r["blocking_issue"] in {"source_not_checked", "unclear_source_identity", "license_unclear", "species_table_missing", "metadata_missing", "pd_control_labels_missing", "pending_manual_review"}
        else:
            assert r["manual_check_status"] == "pending_manual_verification"
            assert r["evidence_quality"] == "unknown"
            assert r["ready_status_after_check"] == "pending"
            assert r["blocking_issue"] == "source_not_checked"

    script = ROOT / "scripts/phase2e_03_summarize_cohort_source_forms.py"
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


def test_summary_writer_schema_is_explicit_and_tab_separated():
    script = ROOT / "scripts/phase2e_03_summarize_cohort_source_forms.py"
    text = script.read_text(encoding="utf-8")
    expected_cols = [
        "dataset_id",
        "total_items",
        "items_verified",
        "species_table_verified",
        "metadata_verified",
        "pd_control_labels_verified",
        "license_verified",
        "run_accession_mapping_verified",
        "ready_candidate_status",
        "major_blocking_issue",
    ]
    for col in expected_cols:
        assert col in text
    assert "SUMMARY_COLUMNS" in text
    assert "to_csv(path, sep=\"\\t\", index=False)" in text
    assert "ready_candidate_status" in text
    assert "major_blocking_issue" in text


def test_summary_columns_exact_and_separated(tmp_path: Path):
    script = ROOT / "scripts/phase2e_03_summarize_cohort_source_forms.py"
    spec = importlib.util.spec_from_file_location("phase2e_03", script)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    expected = [
        "dataset_id",
        "total_items",
        "items_verified",
        "species_table_verified",
        "metadata_verified",
        "pd_control_labels_verified",
        "license_verified",
        "run_accession_mapping_verified",
        "ready_candidate_status",
        "major_blocking_issue",
    ]
    assert mod.SUMMARY_COLUMNS == expected

    out = tmp_path / "summary.tsv"
    mod._write_summary_tsv(
        out,
        [
            {
                "dataset_id": "INTEGRATED_US_MULTICOHORT_PD",
                "total_items": "10",
                "items_verified": "0",
                "species_table_verified": "no",
                "metadata_verified": "no",
                "pd_control_labels_verified": "no",
                "license_verified": "no",
                "run_accession_mapping_verified": "no",
                "ready_candidate_status": "no",
                "major_blocking_issue": "source_not_checked",
            }
        ],
    )
    df = pd.read_csv(out, sep="\t", dtype=str)
    assert list(df.columns) == expected
    assert "ready_candidate_statu" not in df.columns
    assert df.loc[0, "ready_candidate_status"] == "no"
    assert df.loc[0, "major_blocking_issue"] == "source_not_checked"


def test_old_id_absent_from_active_metadata_and_present_in_alias_only():
    active_files = [
        ROOT / "metadata/phase2e_cohort_source_search_targets.tsv",
        ROOT / "metadata/phase2e_cohort_specific_verification_forms.tsv",
        ROOT / "metadata/phase2e_replication_source_triage.tsv",
        ROOT / "metadata/phase2e_replication_manual_download_plan.tsv",
        ROOT / "metadata/phase2e_replication_expected_files.tsv",
        ROOT / "metadata/phase2e_replication_cohort_verification_matrix.tsv",
    ]
    for p in active_files:
        text = p.read_text(encoding="utf-8")
        assert "INTEGRATED_US_MULTICOOHORT_PD" not in text
        assert "INTEGRATED_US_MULTICOHORT_PD" in text

    alias_text = (ROOT / "metadata/phase2e_dataset_id_aliases.tsv").read_text(encoding="utf-8")
    assert "INTEGRATED_US_MULTICOOHORT_PD" in alias_text
    assert "INTEGRATED_US_MULTICOHORT_PD" in alias_text


def test_integrated_us_post_download_status_and_skip_marker():
    status_file = ROOT / "data/phase2e/reports/phase2e_integrated_us_source_download_status.txt"
    summary_file = (
        ROOT
        / "data/phase2e/manual_sources/integrated_us_multicohort_pd/manifests/download_summary.txt"
    )
    skipped_marker = (
        ROOT
        / "data/phase2e/manual_sources/integrated_us_multicohort_pd/zenodo/PD_Metagenomic_Analysis_2b2d81f.zip.partial_skipped"
    )

    assert status_file.exists()
    status_text = status_file.read_text(encoding="utf-8")
    assert "DOWNLOAD_STATUS = DOWNLOAD_COMPLETE_WITH_WARNINGS" in status_text

    assert summary_file.exists()
    summary_text = summary_file.read_text(encoding="utf-8").lower()
    assert "raw_fastq_sra_downloaded = no" in summary_text

    assert skipped_marker.exists()
