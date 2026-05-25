from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return list(r.fieldnames or []), [dict(x) for x in r]


def test_expected_source_files_table_exists_and_has_required_fields():
    path = ROOT / "metadata/phase2c_wallen_expected_source_files.tsv"
    assert path.exists()
    cols, rows = _read_tsv(path)
    for c in [
        "file_id",
        "expected_file_name",
        "source_type",
        "expected_md5",
        "manual_source",
        "local_expected_path",
        "download_mode",
        "download_status",
        "inspection_status",
    ]:
        assert c in cols
    names = {r["expected_file_name"] for r in rows}
    assert "Source_Data_24Oct2022.xlsx" in names
    assert "Supplementary_Code_24Oct2022.zip" in names
    md5s = {r["expected_md5"] for r in rows}
    assert "4672da6a00ae951281441dd0a1620fdd" in md5s
    assert "732c6273349d7a2efc0bc7d13d0c7d69" in md5s


def test_runbook_and_scripts_exist():
    assert (ROOT / "docs/PHASE2C_WALLEN_SOURCE_DOWNLOAD_AND_INSPECTION_RUNBOOK.md").exists()
    for p in [
        "scripts/phase2c_02_check_wallen_source_files.py",
        "scripts/phase2c_03_inspect_wallen_source_workbook.py",
        "scripts/phase2c_04_inspect_wallen_code_zip.py",
        "scripts/phase2c_05_prepare_wallen_sheet_extraction_plan.py",
    ]:
        assert (ROOT / p).exists()


def test_no_forbidden_download_or_heavy_tool_commands_in_new_scripts():
    forbidden_exec_patterns = [
        "wget ",
        "curl ",
        "prefetch",
        "fasterq-dump",
        "sra-tools",
        "subprocess.run([\"hmmsearch",
        "subprocess.run([\"hmmbuild",
        "subprocess.run([\"kraken",
        "subprocess.run([\"humann",
        "subprocess.run([\"metaphlan",
        "subprocess.run([\"bowtie2",
    ]
    for p in [
        ROOT / "scripts/phase2c_02_check_wallen_source_files.py",
        ROOT / "scripts/phase2c_03_inspect_wallen_source_workbook.py",
        ROOT / "scripts/phase2c_04_inspect_wallen_code_zip.py",
        ROOT / "scripts/phase2c_05_prepare_wallen_sheet_extraction_plan.py",
        ]:
        text = p.read_text(encoding="utf-8").lower()
        for token in forbidden_exec_patterns:
            assert token not in text, f"Forbidden token '{token}' found in {p}"


def test_candidate_sheet_mapping_template_exists_or_can_be_used_later():
    path = ROOT / "metadata/phase2c_wallen_candidate_sheet_mapping_template.tsv"
    assert path.exists()
    cols, _ = _read_tsv(path)
    for c in [
        "standard_target_file",
        "candidate_sheet_name",
        "candidate_table_type",
        "confidence",
        "manual_review_status",
        "local_output_path_if_approved",
        "required_for_phase2d",
        "required_for_phase3",
        "notes",
    ]:
        assert c in cols
