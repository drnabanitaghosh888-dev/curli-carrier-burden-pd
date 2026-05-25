#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


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


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), [dict(r) for r in reader]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    manifest_path = root / "metadata" / "pd_metagenome_sources.tsv"
    report_dir = root / "data" / "phase2a" / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    report_tsv = report_dir / "phase2a_manifest_validation.tsv"
    status_txt = report_dir / "phase2a_status.txt"

    rows_out: list[dict[str, str]] = []
    warnings: list[str] = []
    errors: list[str] = []

    if not manifest_path.exists():
        errors.append(f"missing manifest: {manifest_path}")
        cols, rows = [], []
    else:
        cols, rows = read_tsv(manifest_path)

    missing_cols = [c for c in REQUIRED_COLUMNS if c not in cols]
    if missing_cols:
        errors.append(f"missing required columns: {missing_cols}")

    has_high_priority_shotgun_pd = False
    for r in rows:
        seq_type = (r.get("sequencing_type", "") or "").strip().lower()
        disease = (r.get("disease", "") or "").strip().lower()
        priority = (r.get("priority", "") or "").strip().lower()
        status = (r.get("download_status", "") or "").strip()
        intended_use = (r.get("intended_phase2_use", "") or "").strip().lower()

        if (
            "shotgun" in seq_type
            and "parkinson" in disease
            and priority in {"highest", "high"}
        ):
            has_high_priority_shotgun_pd = True

        if status not in VALID_DOWNLOAD_STATUS:
            errors.append(
                f"invalid download_status for {r.get('dataset_id', 'UNKNOWN')}: {status}"
            )

        if "16s" in seq_type and "raw" in intended_use and "csg" in intended_use:
            errors.append(
                f"16S-only dataset marked for raw csg screening: {r.get('dataset_id', 'UNKNOWN')}"
            )

    if not has_high_priority_shotgun_pd:
        errors.append("no high-priority shotgun PD dataset found")

    if any((r.get("download_status", "") == "candidate_pending_manual_verification") for r in rows):
        warnings.append("one or more candidate datasets require manual verification")

    if errors:
        status = "FAIL"
    elif warnings:
        status = "PASS_WITH_WARNINGS"
    else:
        status = "PASS"

    for msg in errors:
        rows_out.append({"level": "ERROR", "check": "manifest", "message": msg})
    for msg in warnings:
        rows_out.append({"level": "WARNING", "check": "manifest", "message": msg})
    if not errors and not warnings:
        rows_out.append({"level": "OK", "check": "manifest", "message": "all checks passed"})

    with report_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["level", "check", "message"], delimiter="\t")
        writer.writeheader()
        for row in rows_out:
            writer.writerow(row)

    status_txt.write_text(status + "\n", encoding="utf-8")
    print(f"[DONE] wrote: {report_tsv}")
    print(f"[DONE] wrote: {status_txt}")
    print(f"[STATUS] {status}")
    return 0 if status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
