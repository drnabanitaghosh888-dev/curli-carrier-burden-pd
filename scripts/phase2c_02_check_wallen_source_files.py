#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
from pathlib import Path


def md5sum(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        while True:
            data = fh.read(chunk_size)
            if not data:
                break
            h.update(data)
    return h.hexdigest()


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    source_tsv = root / "metadata/phase2c_wallen_expected_source_files.tsv"
    out_dir = root / "data/phase2c/reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "phase2c_wallen_source_file_check.tsv"
    out_status = out_dir / "phase2c_wallen_source_file_status.txt"

    with source_tsv.open("r", encoding="utf-8") as h:
        rows = list(csv.DictReader(h, delimiter="\t"))

    out_rows = []
    status = "PASS"
    for r in rows:
        expected_md5 = (r.get("expected_md5") or "").strip().lower()
        local_path = root / (r.get("local_expected_path") or "")
        exists = local_path.exists()
        size = local_path.stat().st_size if exists else 0
        computed_md5 = md5sum(local_path).lower() if exists else ""

        if not exists:
            file_status = "missing"
            status = "PASS_WITH_WARNINGS"
        elif expected_md5 and computed_md5 == expected_md5:
            file_status = "present_md5_match"
        elif expected_md5 and computed_md5 != expected_md5:
            file_status = "present_md5_mismatch"
            status = "FAIL"
        else:
            file_status = "present_no_md5_expected"

        out_rows.append(
            {
                "file_id": r.get("file_id", ""),
                "local_expected_path": str(local_path.relative_to(root)),
                "exists": "yes" if exists else "no",
                "file_size_bytes": str(size),
                "expected_md5": expected_md5,
                "computed_md5": computed_md5,
                "file_status": file_status,
            }
        )

    with out_tsv.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(
            h,
            fieldnames=[
                "file_id",
                "local_expected_path",
                "exists",
                "file_size_bytes",
                "expected_md5",
                "computed_md5",
                "file_status",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for row in out_rows:
            w.writerow(row)

    out_status.write_text(status + "\n", encoding="utf-8")
    print(f"[DONE] wrote: {out_tsv}")
    print(f"[DONE] wrote: {out_status}")
    print(f"[STATUS] {status}")
    return 0 if status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
