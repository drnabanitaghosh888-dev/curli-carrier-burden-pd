#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        return [dict(x) for x in r]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    expected_path = root / "metadata/phase2e_replication_expected_files.tsv"
    out_dir = root / "data/phase2e/reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_inv = out_dir / "phase2e_local_replication_processed_table_inventory.tsv"
    out_status = out_dir / "phase2e_local_file_check_status.txt"

    if not expected_path.exists():
        out_status.write_text("FAIL\nmissing metadata/phase2e_replication_expected_files.tsv\n", encoding="utf-8")
        print(f"[DONE] {out_status}")
        return 1

    rows = _read_tsv(expected_path)
    out_rows = []
    missing_n = 0
    present_n = 0
    for r in rows:
        rel = r.get("local_expected_path", "")
        p = root / rel
        exists = p.exists()
        if exists:
            present_n += 1
        else:
            missing_n += 1
        out_rows.append(
            {
                "dataset_id": r.get("dataset_id", ""),
                "expected_file": r.get("expected_file", ""),
                "local_expected_path": rel,
                "exists": "yes" if exists else "no",
                "required_for_replication_analysis": r.get("required_for_replication_analysis", ""),
                "required_for_raw_read_validation": r.get("required_for_raw_read_validation", ""),
                "notes": r.get("notes", ""),
            }
        )

    with out_inv.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(
            h,
            fieldnames=[
                "dataset_id",
                "expected_file",
                "local_expected_path",
                "exists",
                "required_for_replication_analysis",
                "required_for_raw_read_validation",
                "notes",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(out_rows)

    status = "PASS" if missing_n == 0 else "PASS_WITH_WARNINGS"
    out_status.write_text(
        f"{status}\nmissing_files={missing_n}\npresent_files={present_n}\n",
        encoding="utf-8",
    )
    print(f"[DONE] {out_inv}")
    print(f"[DONE] {out_status}")
    print(f"[STATUS] {status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
