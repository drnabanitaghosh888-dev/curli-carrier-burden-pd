from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_phase2h_nishiwaki_files_exist():
    required = [
        ROOT / "metadata/phase2h_nishiwaki_source_verification_targets.tsv",
        ROOT / "metadata/phase2h_nishiwaki_expected_processed_files.tsv",
        ROOT / "metadata/phase2h_nishiwaki_candidate_source_identity.tsv",
        ROOT / "metadata/phase2h_nishiwaki_processed_table_availability_checklist.tsv",
        ROOT / "metadata/phase2h_nishiwaki_manual_download_plan.tsv",
        ROOT / "docs/PHASE2H_NISHIWAKI_SOURCE_VERIFICATION_RUNBOOK.md",
        ROOT / "data/phase2h/reports/phase2h_nishiwaki_source_verification_status.txt",
        ROOT / "data/phase2h/reports/phase2h_nishiwaki_source_verification_report.tsv",
        ROOT / "data/phase2h/reports/phase2h_nishiwaki_source_verification_report.md",
        ROOT / "data/phase2h/reports/phase2h_nishiwaki_manual_download_plan.md",
        ROOT / "tests/test_phase2h_nishiwaki_source_verification_scaffold.py",
    ]
    missing = [str(p) for p in required if not p.exists()]
    assert not missing


def test_nishiwaki_dataset_and_required_processed_items_present():
    identity = (ROOT / "metadata/phase2h_nishiwaki_candidate_source_identity.tsv").read_text(encoding="utf-8")
    checklist = (ROOT / "metadata/phase2h_nishiwaki_processed_table_availability_checklist.tsv").read_text(encoding="utf-8")

    assert "NISHIWAKI_MULTICOUNTRY_PD" in identity
    assert "species_abundance_table" in checklist
    assert "sample_metadata_table" in checklist
    assert "pd_control_label_field" in checklist


def test_no_rawread_download_heavy_or_phase3_commands_in_phase2h_assets():
    combined = "\n".join(
        [
            (ROOT / "docs/PHASE2H_NISHIWAKI_SOURCE_VERIFICATION_RUNBOOK.md").read_text(encoding="utf-8", errors="ignore").lower(),
            (ROOT / "metadata/phase2h_nishiwaki_source_verification_targets.tsv").read_text(encoding="utf-8", errors="ignore").lower(),
            (ROOT / "metadata/phase2h_nishiwaki_processed_table_availability_checklist.tsv").read_text(encoding="utf-8", errors="ignore").lower(),
            (ROOT / "metadata/phase2h_nishiwaki_manual_download_plan.tsv").read_text(encoding="utf-8", errors="ignore").lower(),
            (ROOT / "data/phase2h/reports/phase2h_nishiwaki_manual_download_plan.md").read_text(encoding="utf-8", errors="ignore").lower(),
        ]
    )

    forbidden_command_patterns = [
        "wget ",
        "prefetch ",
        "fasterq-dump ",
        "hmmsearch ",
        "hmmbuild ",
        "kraken ",
        "humann ",
        "metaphlan ",
        "bowtie2 ",
        "blast ",
        "diamond ",
        "mafft ",
        "python scripts/phase3",
    ]
    for token in forbidden_command_patterns:
        assert token not in combined


def test_ready_for_manual_download_not_asserted_without_required_verification():
    status = (ROOT / "data/phase2h/reports/phase2h_nishiwaki_source_verification_status.txt").read_text(
        encoding="utf-8", errors="ignore"
    )
    report = (ROOT / "data/phase2h/reports/phase2h_nishiwaki_source_verification_report.tsv").read_text(
        encoding="utf-8", errors="ignore"
    )
    assert "READY_FOR_MANUAL_DOWNLOAD = NO" in status
    assert "\tno\tno\tno\t" in report
