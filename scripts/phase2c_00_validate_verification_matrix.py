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
    "vagal_body_first_metadata_status",
    "run_accession_mapping_status",
    "source_url_or_repository",
    "accession_or_identifier",
    "access_restriction",
    "manual_verification_status",
    "ready_for_phase2d_processed_analysis",
    "ready_for_phase3_raw_read_subset_selection",
    "priority",
    "notes",
]


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return list(r.fieldnames or []), [dict(x) for x in r]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    m = root / "metadata/phase2c_processed_table_verification_matrix.tsv"
    v = root / "metadata/phase2c_status_vocabulary.tsv"
    out_dir = root / "data/phase2c/reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "phase2c_verification_matrix_validation.tsv"
    out_status = out_dir / "phase2c_status.txt"

    errors = []
    warnings = []
    out_rows = []

    if not m.exists():
        errors.append(f"missing file: {m}")
        cols, rows = [], []
    else:
        cols, rows = read_tsv(m)
    if not v.exists():
        errors.append(f"missing file: {v}")
        voc_cols, voc_rows = [], []
    else:
        voc_cols, voc_rows = read_tsv(v)

    missing_cols = [c for c in REQ_COLS if c not in cols]
    if missing_cols:
        errors.append(f"missing required columns: {missing_cols}")

    allowed_status = {r.get("status", "") for r in voc_rows}
    status_fields = [
        "processed_species_table_status",
        "processed_genus_table_status",
        "processed_gene_family_table_status",
        "processed_pathway_table_status",
        "sample_metadata_status",
        "pd_control_label_status",
        "clinical_metadata_status",
        "vagal_body_first_metadata_status",
        "run_accession_mapping_status",
        "manual_verification_status",
        "ready_for_phase2d_processed_analysis",
        "ready_for_phase3_raw_read_subset_selection",
    ]

    if rows:
        ids = {r.get("dataset_id", "") for r in rows}
        if "WALLEN_PRJNA834801" not in ids:
            errors.append("WALLEN_PRJNA834801 missing")

        reps = [r for r in rows if "replication_candidate" in (r.get("cohort_role", ""))]
        if len(reps) < 3:
            errors.append("fewer than three replication candidates")

        for r in rows:
            ds = r.get("dataset_id", "UNKNOWN")
            seq = (r.get("sequencing_type", "") or "").lower()
            if "16s" in ds.lower() or "16s" in seq:
                errors.append(f"16S-only dataset present in active matrix: {ds}")
            if seq != "shotgun_metagenomics":
                errors.append(f"non-shotgun active dataset: {ds}")

            for f in status_fields:
                val = r.get(f, "")
                if val not in allowed_status:
                    errors.append(f"status value '{val}' for {f} not in vocabulary ({ds})")

            # Phase2D readiness constraints
            if r.get("ready_for_phase2d_processed_analysis") == "yes":
                if r.get("processed_species_table_status") != "available":
                    errors.append(f"{ds}: phase2d=yes but species table not available")
                if r.get("sample_metadata_status") != "available":
                    errors.append(f"{ds}: phase2d=yes but sample metadata not available")
                if r.get("pd_control_label_status") != "available":
                    errors.append(f"{ds}: phase2d=yes but pd/control labels not available")

            # Phase3 readiness constraints
            if r.get("ready_for_phase3_raw_read_subset_selection") == "yes":
                if r.get("run_accession_mapping_status") != "available":
                    errors.append(f"{ds}: phase3=yes but run accession mapping not available")
                if r.get("sample_metadata_status") != "available":
                    errors.append(f"{ds}: phase3=yes but sample metadata not available")
                if r.get("pd_control_label_status") != "available":
                    errors.append(f"{ds}: phase3=yes but pd/control labels not available")

            # Candidate cohorts should remain pending until evidence exists
            if "candidate" in ds.lower() and r.get("manual_verification_status") != "candidate_pending_manual_verification":
                warnings.append(f"{ds}: candidate cohort should remain candidate_pending_manual_verification")

    status = "PASS"
    if errors:
        status = "FAIL"
    elif warnings:
        status = "PASS_WITH_WARNINGS"

    for e in errors:
        out_rows.append({"level": "ERROR", "check": "phase2c_matrix", "message": e})
    for w in warnings:
        out_rows.append({"level": "WARNING", "check": "phase2c_matrix", "message": w})
    if not errors and not warnings:
        out_rows.append({"level": "OK", "check": "phase2c_matrix", "message": "all checks passed"})

    with out_tsv.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["level", "check", "message"], delimiter="\t")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    out_status.write_text(status + "\n", encoding="utf-8")
    print(f"[DONE] wrote: {out_tsv}")
    print(f"[DONE] wrote: {out_status}")
    print(f"[STATUS] {status}")
    return 0 if status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
