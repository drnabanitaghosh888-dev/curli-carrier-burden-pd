#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


SOURCES_REQUIRED = [
    "dataset_id",
    "study_name",
    "cohort_role",
    "disease",
    "sample_type",
    "sequencing_type",
    "geography",
    "accession_or_source",
    "processed_table_source",
    "processed_table_type",
    "raw_data_source",
    "sample_metadata_source",
    "clinical_metadata_source",
    "vagal_body_first_metadata_status",
    "expected_sample_count",
    "local_expected_path",
    "access_status",
    "manual_verification_status",
    "download_status",
    "priority",
    "notes",
]


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return list(r.fieldnames or []), [dict(x) for x in r]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    m_sources = root / "metadata/phase2b_processed_table_sources.tsv"
    m_schema = root / "metadata/phase2b_required_table_schema.tsv"
    m_queue = root / "metadata/phase2b_cohort_verification_queue.tsv"
    m_excluded = root / "metadata/phase2b_excluded_or_context_only_datasets.tsv"
    report_dir = root / "data/phase2b/reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = report_dir / "phase2b_processed_table_plan_validation.tsv"
    out_status = report_dir / "phase2b_status.txt"

    rows_out: list[dict[str, str]] = []
    errors: list[str] = []
    warnings: list[str] = []

    for p in [m_sources, m_schema, m_queue, m_excluded]:
        if not p.exists():
            errors.append(f"missing required file: {p}")

    if m_sources.exists():
        cols, rows = read_tsv(m_sources)
    else:
        cols, rows = [], []

    missing_cols = [c for c in SOURCES_REQUIRED if c not in cols]
    if missing_cols:
        errors.append(f"phase2b_processed_table_sources missing columns: {missing_cols}")

    excluded_ids = set()
    if m_excluded.exists():
        _, ex_rows = read_tsv(m_excluded)
        excluded_ids = {r.get("dataset_id", "") for r in ex_rows}

    if rows:
        discovery = [r for r in rows if r.get("cohort_role") == "discovery_anchor"]
        if not discovery:
            errors.append("no discovery_anchor dataset found")

        replication = [
            r for r in rows if "replication_candidate" in (r.get("cohort_role") or "")
        ]
        if len(replication) < 2:
            errors.append("fewer than two replication candidates found")

        wallen = [r for r in rows if r.get("dataset_id") == "WALLEN_PRJNA834801"]
        if not wallen:
            errors.append("WALLEN_PRJNA834801 missing")
        else:
            if (wallen[0].get("priority") or "").lower() != "highest":
                errors.append("WALLEN_PRJNA834801 is not priority highest")

        for r in rows:
            ds = r.get("dataset_id", "UNKNOWN")
            seq = (r.get("sequencing_type") or "").lower()
            notes = (r.get("notes") or "").lower()
            # Active registry must be shotgun-only.
            if "16s" in seq:
                errors.append(f"16S-only dataset present in active Phase 2B registry: {ds}")
            if "shotgun_metagenomics" not in seq:
                errors.append(f"active dataset is not shotgun_metagenomics: {ds}")

            # If a 16S dataset id appears in active rows, it must only be excluded list.
            if ds in excluded_ids:
                errors.append(f"excluded/context-only dataset appears in active registry: {ds}")

            if (r.get("download_status") or "").lower() == "downloaded":
                expected = (r.get("local_expected_path") or "").strip()
                if not expected:
                    errors.append(f"downloaded row missing local_expected_path: {ds}")

            forbidden = ["fastq", "sra", "prefetch", "fasterq", "kraken", "metaphlan", "humann", "bowtie2"]
            hay = " ".join(
                [
                    r.get("processed_table_source", ""),
                    r.get("raw_data_source", ""),
                    r.get("notes", ""),
                ]
            ).lower()
            if any(tok in hay for tok in forbidden):
                warnings.append(f"check wording for heavy-step token in sources/notes: {ds}")

            # Fail if any row indicates active/raw screening or primary association testing in this phase.
            if "raw_csg_screening" in hay or "primary_association_testing" in hay:
                errors.append(f"forbidden active marker in Phase 2B active registry row: {ds}")

            manual_status = r.get("manual_verification_status", "")
            accession = r.get("accession_or_source", "")
            if "candidate" in ds.lower() and "pending_manual_verification" not in manual_status:
                warnings.append(
                    f"candidate cohort should remain pending_manual_verification until verified: {ds}"
                )
            if "candidate" in ds.lower() and accession == "candidate_pending_manual_verification":
                pass

    # 16S-only datasets must exist only in excluded/context table when present.
    if m_excluded.exists():
        _, ex_rows = read_tsv(m_excluded)
        for r in ex_rows:
            ds = r.get("dataset_id", "")
            st = (r.get("study_type", "") or "").lower()
            not_allowed = (r.get("not_allowed_use", "") or "").lower()
            if "16s" in st:
                required_forbidden = [
                    "curli_gene_detection",
                    "csg_operon_analysis",
                    "raw_csg_screening",
                    "strain_reconstruction",
                    "primary_association_testing",
                ]
                for token in required_forbidden:
                    if token not in not_allowed:
                        errors.append(
                            f"excluded 16S dataset missing not_allowed_use token '{token}': {ds}"
                        )

    status = "PASS"
    if errors:
        status = "FAIL"
    elif warnings:
        status = "PASS_WITH_WARNINGS"

    for msg in errors:
        rows_out.append({"level": "ERROR", "check": "phase2b_plan", "message": msg})
    for msg in warnings:
        rows_out.append({"level": "WARNING", "check": "phase2b_plan", "message": msg})
    if not errors and not warnings:
        rows_out.append({"level": "OK", "check": "phase2b_plan", "message": "all checks passed"})

    with out_tsv.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["level", "check", "message"], delimiter="\t")
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    out_status.write_text(status + "\n", encoding="utf-8")

    print(f"[DONE] wrote: {out_tsv}")
    print(f"[DONE] wrote: {out_status}")
    print(f"[STATUS] {status}")
    return 0 if status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
