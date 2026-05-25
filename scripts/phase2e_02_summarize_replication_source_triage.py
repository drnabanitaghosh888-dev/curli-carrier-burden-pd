#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as h:
        reader = csv.DictReader(h, delimiter="\t")
        return [dict(r) for r in reader]


def _write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    triage_path = root / "metadata/phase2e_replication_source_triage.tsv"
    out_dir = root / "data/phase2e/reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_summary = out_dir / "phase2e_replication_source_triage_summary.tsv"
    out_status = out_dir / "phase2e_source_triage_status.txt"

    if not triage_path.exists():
        _write_tsv(
            out_summary,
            [{"metric": "error", "value": f"missing {triage_path}"}],
            ["metric", "value"],
        )
        out_status.write_text("FAIL\n", encoding="utf-8")
        print(f"[DONE] {out_summary}")
        print(f"[DONE] {out_status}")
        return 1

    rows = _read_tsv(triage_path)
    total = len(rows)
    ready_dl = sum(1 for r in rows if r.get("ready_for_manual_download", "").lower() == "yes")
    ready_rep = sum(1 for r in rows if r.get("ready_for_processed_replication_analysis", "").lower() == "yes")
    ready_raw = sum(1 for r in rows if r.get("ready_for_raw_read_validation", "").lower() == "yes")
    pending = sum(1 for r in rows if r.get("source_status", "") == "pending_manual_source_check")

    blocking = {}
    for r in rows:
        b = r.get("blocking_issue", "unspecified")
        blocking[b] = blocking.get(b, 0) + 1

    summary_rows = [
        {"metric": "n_cohorts", "value": str(total)},
        {"metric": "n_pending_manual_source_check", "value": str(pending)},
        {"metric": "n_ready_for_manual_download", "value": str(ready_dl)},
        {"metric": "n_ready_for_processed_replication_analysis", "value": str(ready_rep)},
        {"metric": "n_ready_for_raw_read_validation", "value": str(ready_raw)},
    ]
    for k, v in sorted(blocking.items()):
        summary_rows.append({"metric": f"blocking_issue::{k}", "value": str(v)})

    _write_tsv(out_summary, summary_rows, ["metric", "value"])
    status = "PASS_WITH_WARNINGS" if ready_rep == 0 else "PASS"
    out_status.write_text(status + "\n", encoding="utf-8")
    print(f"[DONE] {out_summary}")
    print(f"[DONE] {out_status}")
    print(f"[STATUS] {status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
