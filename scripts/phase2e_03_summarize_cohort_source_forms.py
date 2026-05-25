#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


SUMMARY_COLUMNS = [
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


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return [dict(x) for x in r]


def _write_summary_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    df = pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
    df.to_csv(path, sep="\t", index=False)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    forms_path = root / "metadata/phase2e_cohort_specific_verification_forms.tsv"
    out_dir = root / "data/phase2e/reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_summary = out_dir / "phase2e_cohort_source_form_summary.tsv"
    out_status = out_dir / "phase2e_cohort_source_form_status.txt"

    if not forms_path.exists():
        _write_summary_tsv(
            out_summary,
            [
                {
                    "dataset_id": "NA",
                    "total_items": "0",
                    "items_verified": "0",
                    "species_table_verified": "no",
                    "metadata_verified": "no",
                    "pd_control_labels_verified": "no",
                    "license_verified": "no",
                    "run_accession_mapping_verified": "no",
                    "ready_candidate_status": "no",
                    "major_blocking_issue": "missing_form_file",
                }
            ],
        )
        out_status.write_text("FAIL\n", encoding="utf-8")
        print(f"[DONE] {out_summary}")
        print(f"[DONE] {out_status}")
        return 1

    rows = _read_tsv(forms_path)
    by_ds: dict[str, list[dict[str, str]]] = {}
    for r in rows:
        by_ds.setdefault(r.get("dataset_id", "UNKNOWN"), []).append(r)

    summary_rows = []
    ready_candidates = 0
    for ds, ds_rows in sorted(by_ds.items()):
        total = len(ds_rows)
        verified = sum(1 for r in ds_rows if r.get("manual_check_status") == "verified")
        def _flag(item: str) -> str:
            for r in ds_rows:
                if r.get("verification_item") == item and r.get("manual_check_status") == "verified":
                    return "yes"
            return "no"

        species_ok = _flag("species_abundance_table")
        meta_ok = _flag("sample_metadata")
        labels_ok = _flag("pd_control_labels")
        license_ok = _flag("license_or_access_terms")
        runmap_ok = _flag("run_accession_mapping")

        ready = "yes" if all(x == "yes" for x in [species_ok, meta_ok, labels_ok, license_ok]) else "no"
        if ready == "yes":
            ready_candidates += 1

        blockers = [r.get("blocking_issue", "") for r in ds_rows if r.get("blocking_issue", "source_not_checked") != ""]
        major_blocker = blockers[0] if blockers else "none"

        summary_rows.append(
            {
                "dataset_id": ds,
                "total_items": str(total),
                "items_verified": str(verified),
                "species_table_verified": species_ok,
                "metadata_verified": meta_ok,
                "pd_control_labels_verified": labels_ok,
                "license_verified": license_ok,
                "run_accession_mapping_verified": runmap_ok,
                "ready_candidate_status": ready,
                "major_blocking_issue": major_blocker,
            }
        )

    _write_summary_tsv(out_summary, summary_rows)

    status = "PASS_WITH_WARNINGS" if ready_candidates == 0 else "PASS"
    out_status.write_text(status + "\n", encoding="utf-8")
    print(f"[DONE] {out_summary}")
    print(f"[DONE] {out_status}")
    print(f"[STATUS] {status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
