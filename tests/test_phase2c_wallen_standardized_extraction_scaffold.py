from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return list(r.fieldnames or []), [dict(x) for x in r]


def test_scripts_exist():
    for p in [
        ROOT / "scripts/phase2c_06_extract_wallen_standardized_tables.py",
        ROOT / "scripts/phase2c_07_validate_wallen_standardized_tables.py",
        ROOT / "tests/test_phase2c_wallen_standardized_extraction_scaffold.py",
    ]:
        assert p.exists(), f"Missing scaffold file: {p}"


def test_mapping_has_corrected_approved_sheet_names():
    cols, rows = _read_tsv(ROOT / "metadata/phase2c_wallen_candidate_sheet_mapping_template.tsv")
    required_cols = [
        "standard_target_file",
        "candidate_sheet_name",
        "candidate_table_type",
        "confidence",
        "manual_review_status",
    ]
    for c in required_cols:
        assert c in cols

    by_target = {r["standard_target_file"]: r for r in rows}
    assert by_target["sample_metadata.tsv"]["candidate_sheet_name"] == "subject_metadata"
    assert by_target["clinical_metadata.tsv"]["candidate_sheet_name"] == "subject_metadata"
    assert by_target["species_abundance.tsv"]["candidate_sheet_name"] == "metaphlan_rel_ab"
    assert by_target["genus_abundance.tsv"]["candidate_sheet_name"] == "metaphlan_rel_ab"
    assert by_target["gene_family_abundance.tsv"]["candidate_sheet_name"] == "humann_KO_group_counts"
    assert by_target["pathway_abundance.tsv"]["candidate_sheet_name"] == "humann_pathway_counts"
    assert by_target["run_accession_mapping.tsv"]["candidate_sheet_name"] == "not_available_in_source_workbook"
    assert by_target["run_accession_mapping.tsv"]["candidate_sheet_name"] != "Directory"


def test_extraction_requires_execute_flag():
    text = (ROOT / "scripts/phase2c_06_extract_wallen_standardized_tables.py").read_text(encoding="utf-8")
    assert "--execute" in text
    assert "if not args.execute" in text
    assert "DRY_RUN_ONLY" in text


def test_run_accession_mapping_not_planned_as_extractable_from_workbook():
    cols, rows = _read_tsv(ROOT / "metadata/phase2c_wallen_candidate_sheet_mapping_template.tsv")
    assert "standard_target_file" in cols
    assert "candidate_sheet_name" in cols
    by_target = {r["standard_target_file"]: r for r in rows}
    run_map = by_target["run_accession_mapping.tsv"]
    assert run_map["candidate_sheet_name"] == "not_available_in_source_workbook"
    assert run_map["candidate_sheet_name"] != "Directory"


def test_script_enforces_single_not_available_run_accession_row():
    text = (ROOT / "scripts/phase2c_06_extract_wallen_standardized_tables.py").read_text(encoding="utf-8")
    assert 'if target == "run_accession_mapping.tsv":' in text
    assert 'status="PLANNED_NOT_AVAILABLE"' in text
    assert "run_accession_mapping.NOT_AVAILABLE_YET.txt" in text
    assert "blocks Phase 3 but not Phase 2D" in text
    assert "continue" in text


def test_plan_written_with_dataframe_and_expected_columns():
    text = (ROOT / "scripts/phase2c_06_extract_wallen_standardized_tables.py").read_text(encoding="utf-8")
    expected_cols = [
        "standard_target_file",
        "source_sheet",
        "output_path",
        "n_rows",
        "n_columns",
        "status",
        "notes",
    ]
    for col in expected_cols:
        assert col in text
    assert "plan_df = pd.DataFrame(rows, columns=PLAN_COLUMNS)" in text
    assert "plan_df.to_csv(plan_path, sep=\"\\t\", index=False)" in text


def test_plan_rows_and_tsv_shape_for_run_accession(tmp_path: Path):
    script_path = ROOT / "scripts/phase2c_06_extract_wallen_standardized_tables.py"
    spec = importlib.util.spec_from_file_location("phase2c_06", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    mapping_rows = [
        {"standard_target_file": "sample_metadata.tsv", "candidate_sheet_name": "subject_metadata"},
        {"standard_target_file": "run_accession_mapping.tsv", "candidate_sheet_name": "not_available_in_source_workbook"},
    ]
    outdir = ROOT / "data/phase2b/processed_tables/wallen_prjna834801"
    rows = module.build_extraction_plan_rows(mapping_rows=mapping_rows, outdir=outdir, root=ROOT)

    run_rows = [r for r in rows if r["standard_target_file"] == "run_accession_mapping.tsv"]
    assert len(run_rows) == 1
    assert run_rows[0]["status"] == "PLANNED_NOT_AVAILABLE"
    assert run_rows[0]["notes"].startswith("run/sample accession mapping")

    plan_path = tmp_path / "plan.tsv"
    module.write_extraction_plan(rows=rows, plan_path=plan_path)

    df = pd.read_csv(plan_path, sep="\t", dtype=str).fillna("")
    assert list(df.columns) == [
        "standard_target_file",
        "source_sheet",
        "output_path",
        "n_rows",
        "n_columns",
        "status",
        "notes",
    ]
    assert df.shape[1] == 7
    with plan_path.open("r", encoding="utf-8") as handle:
        lines = [ln.rstrip("\n") for ln in handle if ln.strip()]
    for line in lines:
        assert line.count("\t") == 6, f"Malformed TSV row: {line}"

    run_df = df[df["standard_target_file"] == "run_accession_mapping.tsv"]
    assert len(run_df) == 1
    row = run_df.iloc[0]
    assert row["status"] == "PLANNED_NOT_AVAILABLE"
    assert row["notes"].startswith("run/sample accession mapping")
    full_text = plan_path.read_text(encoding="utf-8")
    assert "PLANNED_NOT_AVAILABLErun/sample accession mapping" not in full_text


def test_no_fastq_sra_or_raw_read_commands_present():
    forbidden = [
        "prefetch",
        "fasterq-dump",
        "sra-tools",
        "fastq",
        "hmmsearch",
        "hmmbuild",
        "kraken",
        "metaphlan",
        "humann",
        "bowtie2",
        "blast",
        "diamond",
        "mafft",
    ]
    for p in [
        ROOT / "scripts/phase2c_06_extract_wallen_standardized_tables.py",
        ROOT / "scripts/phase2c_07_validate_wallen_standardized_tables.py",
    ]:
        text = p.read_text(encoding="utf-8").lower()
        for token in forbidden:
            assert token not in text, f"Forbidden token '{token}' found in {p}"
