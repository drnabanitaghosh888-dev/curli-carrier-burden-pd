#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
from pathlib import Path


GENES = ["csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create Phase 1 curli reference report.")
    parser.add_argument("--curated-dir", default="data/phase1/curated")
    parser.add_argument("--hmm-dir", default="data/phase1/hmms")
    parser.add_argument("--validation-dir", default="data/phase1/validation")
    parser.add_argument("--report-dir", default="data/phase1/reports")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(r) for r in reader]


def status_from_checks(checks: list[tuple[str, bool]]) -> str:
    missing = [name for name, ok in checks if not ok]
    if missing:
        return "FAIL"
    return "PASS"


def write_status(path: Path, value: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(value + "\n")


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_markdown(path: Path, lines: list[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    curated_dir = (root / args.curated_dir).resolve()
    hmm_dir = (root / args.hmm_dir).resolve()
    validation_dir = (root / args.validation_dir).resolve()
    report_dir = (root / args.report_dir).resolve()
    report_dir.mkdir(parents=True, exist_ok=True)

    report_md = report_dir / "phase1_curli_reference_report.md"
    report_tsv = report_dir / "phase1_curli_reference_report.tsv"
    status_txt = report_dir / "phase1_status.txt"

    for p in [report_md, report_tsv, status_txt]:
        if p.exists() and not args.force:
            raise FileExistsError(f"Output exists: {p}. Use --force to overwrite.")

    checks: list[tuple[str, bool]] = []
    warnings: list[str] = []

    # Required curated FASTA files
    for gene in GENES:
        checks.append((f"curated_fasta_{gene}", (curated_dir / f"{gene}.curated.faa").exists()))

    # Required per-gene HMM files
    for gene in GENES:
        checks.append((f"hmm_{gene}", (hmm_dir / f"{gene}.hmm").exists()))

    # Required combined artifacts
    combined_hmm = hmm_dir / "curli_csg_combined.hmm"
    checks.append(("combined_hmm", combined_hmm.exists()))
    for ext in [".h3f", ".h3i", ".h3m", ".h3p"]:
        checks.append((f"hmmpress_index{ext}", (hmm_dir / f"curli_csg_combined.hmm{ext}").exists()))

    # Required metadata and validation summary
    curated_meta = curated_dir / "curli_reference_metadata.tsv"
    checks.append(("curated_metadata_table", curated_meta.exists()))
    validation_summary = validation_dir / "hmm_validation_summary.tsv"
    checks.append(("validation_summary_table", validation_summary.exists()))

    # pass/fail baseline
    base_status = status_from_checks(checks)

    # Optional warnings from validation classification
    if validation_summary.exists():
        vrows = read_tsv(validation_summary)
        weak = [r.get("gene", "") for r in vrows if r.get("classification", "").strip().lower() == "weak"]
        borderline = [
            r.get("gene", "")
            for r in vrows
            if r.get("classification", "").strip().lower() == "borderline"
        ]
        if weak:
            warnings.append(f"Weak HMM classification: {', '.join(weak)}")
        if borderline:
            warnings.append(f"Borderline HMM classification: {', '.join(borderline)}")

    final_status = base_status
    if base_status == "PASS" and warnings:
        final_status = "PASS_WITH_WARNINGS"

    # Build TSV report
    tsv_rows: list[dict[str, str]] = []
    for name, ok in checks:
        tsv_rows.append(
            {
                "check_name": name,
                "status": "OK" if ok else "MISSING",
                "detail": "",
            }
        )
    for w in warnings:
        tsv_rows.append(
            {
                "check_name": "warning",
                "status": "WARN",
                "detail": w,
            }
        )
    tsv_fields = ["check_name", "status", "detail"]
    write_tsv(report_tsv, tsv_rows, tsv_fields)

    # Build Markdown report
    now = dt.datetime.now().isoformat(timespec="seconds")
    md_lines = [
        "# Phase 1 Curli Reference Report",
        "",
        f"- Generated: {now}",
        f"- Status: **{final_status}**",
        "",
        "## Required checks",
        "",
        "| check_name | status |",
        "|---|---|",
    ]
    for name, ok in checks:
        md_lines.append(f"| {name} | {'OK' if ok else 'MISSING'} |")
    if warnings:
        md_lines.extend(["", "## Warnings", ""])
        for w in warnings:
            md_lines.append(f"- {w}")
    else:
        md_lines.extend(["", "## Warnings", "", "- None"])

    write_markdown(report_md, md_lines)
    write_status(status_txt, final_status)

    print(f"[DONE] wrote: {report_md}")
    print(f"[DONE] wrote: {report_tsv}")
    print(f"[DONE] wrote: {status_txt}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
