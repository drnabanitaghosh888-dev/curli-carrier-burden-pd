#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


REQ_COLS = [
    "dataset_id",
    "cohort_role",
    "sequencing_type",
    "disease_context",
    "sample_type",
    "processed_species_table_status",
    "processed_genus_table_status",
    "processed_gene_family_table_status",
    "processed_pathway_table_status",
    "sample_metadata_status",
    "pd_control_label_status",
    "clinical_metadata_status",
    "prodromal_or_body_first_metadata_status",
    "run_accession_mapping_status",
    "source_repository_or_url",
    "accession_or_identifier",
    "access_restriction",
    "manual_verification_status",
    "ready_for_replication_processed_analysis",
    "ready_for_raw_read_validation",
    "priority",
    "notes",
]

EXPECTED_COHORTS = {
    "PALACIOS_PRODROMAL_PD",
    "NISHIWAKI_MULTICOUNTRY_PD",
    "MAO_CENTRAL_CHINA_PD",
    "INTEGRATED_US_MULTICOHORT_PD",
}


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return list(r.fieldnames or []), [dict(x) for x in r]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    matrix_path = root / "metadata/phase2e_replication_cohort_verification_matrix.tsv"
    expected_files_path = root / "metadata/phase2e_replication_expected_files.tsv"
    rules_path = root / "metadata/phase2e_replication_readiness_rules.tsv"
    vocab_path = root / "metadata/phase2c_status_vocabulary.tsv"

    out_dir = root / "data/phase2e/reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "phase2e_replication_matrix_validation.tsv"
    out_status = out_dir / "phase2e_status.txt"

    issues: list[dict[str, str]] = []
    warnings = 0
    errors = 0

    if not matrix_path.exists():
        issues.append({"level": "ERROR", "check": "exists", "message": f"missing {matrix_path}"})
        errors += 1
        cols, rows = [], []
    else:
        cols, rows = _read_tsv(matrix_path)
    for p in [expected_files_path, rules_path, vocab_path]:
        if not p.exists():
            issues.append({"level": "ERROR", "check": "exists", "message": f"missing {p}"})
            errors += 1

    missing_cols = [c for c in REQ_COLS if c not in cols]
    if missing_cols:
        issues.append({"level": "ERROR", "check": "columns", "message": f"missing columns: {missing_cols}"})
        errors += 1

    if rows:
        present = {r.get("dataset_id", "") for r in rows}
        missing = sorted(EXPECTED_COHORTS - present)
        if missing:
            issues.append({"level": "ERROR", "check": "cohorts", "message": f"missing cohorts: {missing}"})
            errors += 1

        for r in rows:
            ds = r.get("dataset_id", "UNKNOWN")
            seq = r.get("sequencing_type", "")
            if "16S" in ds or seq.lower() == "16s":
                issues.append({"level": "ERROR", "check": "16s_exclusion", "message": f"16S-only dataset in active matrix: {ds}"})
                errors += 1
            if seq != "shotgun_metagenomics":
                issues.append({"level": "ERROR", "check": "sequencing_type", "message": f"non-shotgun active dataset: {ds}"})
                errors += 1

            # Ready flags can only be yes when core evidence exists
            if r.get("ready_for_replication_processed_analysis") == "yes":
                needs = [
                    ("processed_species_table_status", "available"),
                    ("sample_metadata_status", "available"),
                    ("pd_control_label_status", "available"),
                ]
                for field, val in needs:
                    if r.get(field) != val:
                        issues.append({"level": "ERROR", "check": "processed_ready_gate", "message": f"{ds} marked ready but {field}!={val}"})
                        errors += 1

            if r.get("ready_for_raw_read_validation") == "yes":
                needs = [
                    ("run_accession_mapping_status", "available"),
                    ("sample_metadata_status", "available"),
                    ("pd_control_label_status", "available"),
                ]
                for field, val in needs:
                    if r.get(field) != val:
                        issues.append({"level": "ERROR", "check": "raw_ready_gate", "message": f"{ds} marked raw-read ready but {field}!={val}"})
                        errors += 1

            # Pending verification is expected at Phase 2E-A
            if r.get("manual_verification_status") == "candidate_pending_manual_verification":
                warnings += 1

    status = "PASS"
    if errors:
        status = "FAIL"
    elif warnings:
        status = "PASS_WITH_WARNINGS"
    else:
        status = "PASS"

    if not issues:
        issues.append({"level": "OK", "check": "summary", "message": "matrix scaffold checks passed"})

    with out_tsv.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["level", "check", "message"], delimiter="\t")
        w.writeheader()
        w.writerows(issues)
    out_status.write_text(status + "\n", encoding="utf-8")
    print(f"[DONE] {out_tsv}")
    print(f"[DONE] {out_status}")
    print(f"[STATUS] {status}")
    return 0 if status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
