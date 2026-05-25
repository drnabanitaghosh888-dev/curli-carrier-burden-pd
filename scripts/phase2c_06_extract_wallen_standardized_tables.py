#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


WORKBOOK = "data/phase2c/manual_sources/wallen_prjna834801/Source_Data_24Oct2022.xlsx"
MAPPING = "metadata/phase2c_wallen_candidate_sheet_mapping_template.tsv"
OUTDIR = "data/phase2b/processed_tables/wallen_prjna834801"
REPORT_TSV = "data/phase2c/reports/phase2c_wallen_standardized_extraction_plan.tsv"
REPORT_STATUS = "data/phase2c/reports/phase2c_wallen_standardized_extraction_status.txt"
PLAN_COLUMNS = [
    "standard_target_file",
    "source_sheet",
    "output_path",
    "n_rows",
    "n_columns",
    "status",
    "notes",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare or execute standardized extraction from Wallen workbook.")
    p.add_argument("--execute", action="store_true", help="Run extraction. Default is dry-run planning only.")
    return p.parse_args()


def read_mapping(path: Path) -> list[dict[str, str]]:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("").to_dict(orient="records")


def build_not_available_row(outdir: Path, root: Path, status: str) -> dict[str, str]:
    return {
        "standard_target_file": "run_accession_mapping.tsv",
        "source_sheet": "not_available_in_source_workbook",
        "output_path": str((outdir / "run_accession_mapping.NOT_AVAILABLE_YET.txt").relative_to(root)),
        "n_rows": "NA",
        "n_columns": "NA",
        "status": status,
        "notes": "run/sample accession mapping must be verified externally from BioProject/SRA metadata later; this blocks Phase 3 but not Phase 2D",
    }


def build_extraction_plan_rows(mapping_rows: list[dict[str, str]], outdir: Path, root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for r in mapping_rows:
        target = r["standard_target_file"]
        source_sheet = r["candidate_sheet_name"]
        if target == "run_accession_mapping.tsv":
            # Not extractable from source workbook; represented by one explicit NOT_AVAILABLE row.
            continue
        rows.append(
            {
                "standard_target_file": target,
                "source_sheet": source_sheet,
                "output_path": str((outdir / target).relative_to(root)),
                "n_rows": "NA",
                "n_columns": "NA",
                "status": "PLANNED_DRY_RUN",
                "notes": "dry-run only; no extraction performed",
            }
        )
    rows.append(build_not_available_row(outdir=outdir, root=root, status="PLANNED_NOT_AVAILABLE"))
    return rows


def write_extraction_plan(rows: list[dict[str, str]], plan_path: Path) -> None:
    plan_df = pd.DataFrame(rows, columns=PLAN_COLUMNS)
    plan_df.to_csv(plan_path, sep="\t", index=False)


def write_report(rows: list[dict[str, str]], status: str, root: Path) -> None:
    out_tsv = root / REPORT_TSV
    out_status = root / REPORT_STATUS
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    write_extraction_plan(rows=rows, plan_path=out_tsv)
    out_status.write_text(status + "\n", encoding="utf-8")
    print(f"[DONE] wrote: {out_tsv}")
    print(f"[DONE] wrote: {out_status}")


def _extract_table(df: pd.DataFrame, target: str) -> pd.DataFrame:
    # Required table shaping rules
    if target == "species_abundance.tsv":
        if "clade_name" in df.columns:
            return df[df["clade_name"].astype(str).str.contains("s__", na=False)].copy()
    if target == "genus_abundance.tsv":
        if "clade_name" in df.columns:
            m_g = df["clade_name"].astype(str).str.contains("g__", na=False)
            m_s = df["clade_name"].astype(str).str.contains("s__", na=False)
            return df[m_g & ~m_s].copy()
    return df.copy()


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    workbook = root / WORKBOOK
    mapping_file = root / MAPPING
    outdir = root / OUTDIR
    mapping_rows = read_mapping(mapping_file)

    report_rows: list[dict[str, str]] = []

    if not args.execute:
        report_rows = build_extraction_plan_rows(mapping_rows=mapping_rows, outdir=outdir, root=root)
        write_report(report_rows, "DRY_RUN_ONLY", root)
        return 0

    # Explicit execution mode
    if not workbook.exists():
        report_rows.append(
            {
                "standard_target_file": "ALL",
                "source_sheet": "NA",
                "output_path": str(outdir.relative_to(root)),
                "n_rows": "0",
                "n_columns": "0",
                "status": "MISSING_WORKBOOK",
                "notes": "workbook not found",
            }
        )
        write_report(report_rows, "FAIL_MISSING_WORKBOOK", root)
        return 1

    outdir.mkdir(parents=True, exist_ok=True)

    # Extract approved mappings only
    for r in mapping_rows:
        target = r["standard_target_file"]
        source_sheet = r["candidate_sheet_name"]
        review = r["manual_review_status"]
        out_path = outdir / target

        if target == "run_accession_mapping.tsv":
            continue

        if review != "approved_for_extraction":
            report_rows.append(
                {
                    "standard_target_file": target,
                    "source_sheet": source_sheet,
                    "output_path": str(out_path.relative_to(root)),
                    "n_rows": "0",
                    "n_columns": "0",
                    "status": "SKIPPED_NOT_APPROVED",
                    "notes": "manual_review_status not approved_for_extraction",
                }
            )
            continue

        df = pd.read_excel(workbook, sheet_name=source_sheet)
        out_df = _extract_table(df, target)
        out_df.to_csv(out_path, sep="\t", index=False)
        report_rows.append(
            {
                "standard_target_file": target,
                "source_sheet": source_sheet,
                "output_path": str(out_path.relative_to(root)),
                "n_rows": str(len(out_df)),
                "n_columns": str(len(out_df.columns)),
                "status": "EXTRACTED",
                "notes": "sheet extracted to standardized table",
            }
        )

    # Explicitly do not generate run_accession_mapping.tsv from workbook
    marker = outdir / "run_accession_mapping.NOT_AVAILABLE_YET.txt"
    marker.write_text(
        "run/sample accession mapping is not available in the source workbook and must be verified separately.\n",
        encoding="utf-8",
    )
    report_rows.append(build_not_available_row(outdir=outdir, root=root, status="NOT_AVAILABLE_YET"))

    write_report(report_rows, "EXECUTION_COMPLETE", root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
