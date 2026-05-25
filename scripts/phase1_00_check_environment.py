#!/usr/bin/env python3
from __future__ import annotations

import datetime as dt
import importlib
import os
import shutil
import socket
import subprocess
import sys
import urllib.request
from pathlib import Path


REQUIRED_COMMANDS = [
    "python",
    "mafft",
    "hmmbuild",
    "hmmpress",
    "hmmsearch",
    "seqkit",
    "cd-hit",
    "diamond",
    "blastp",
]

REQUIRED_PYTHON_PACKAGES = [
    "Bio",
    "pandas",
    "numpy",
    "requests",
    "yaml",
    "tqdm",
]


def check_command(cmd: str) -> tuple[str, str, str]:
    path = shutil.which(cmd)
    if path:
        return ("command", cmd, f"OK: {path}")
    return ("command", cmd, "MISSING")


def check_python_package(pkg: str) -> tuple[str, str, str]:
    try:
        importlib.import_module(pkg)
        return ("python_package", pkg, "OK")
    except Exception as exc:  # pragma: no cover
        return ("python_package", pkg, f"MISSING: {exc}")


def check_internet_head(url: str = "https://rest.uniprot.org", timeout: int = 5) -> tuple[str, str, str]:
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:  # nosec B310
            code = getattr(resp, "status", "unknown")
        return ("internet", url, f"OK: HEAD status {code}")
    except Exception as exc:
        return ("internet", url, f"FAILED: {exc}")


def check_dns(hostname: str = "rest.uniprot.org") -> tuple[str, str, str]:
    try:
        resolved = socket.gethostbyname(hostname)
        return ("dns", hostname, f"OK: {resolved}")
    except Exception as exc:
        return ("dns", hostname, f"FAILED: {exc}")


def check_disk_space(path: Path) -> tuple[str, str, str]:
    usage = shutil.disk_usage(path)
    free_gb = usage.free / (1024**3)
    total_gb = usage.total / (1024**3)
    return ("disk", str(path), f"OK: free_gb={free_gb:.2f}; total_gb={total_gb:.2f}")


def check_python_runtime() -> tuple[str, str, str]:
    version = sys.version.replace("\n", " ")
    return ("python_runtime", "python", f"OK: {version}")


def write_tsv(rows: list[tuple[str, str, str]], out_tsv: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8") as handle:
        handle.write("check_type\tname\tstatus\n")
        for check_type, name, status in rows:
            handle.write(f"{check_type}\t{name}\t{status}\n")


def summarize_status(rows: list[tuple[str, str, str]]) -> tuple[int, int]:
    ok = 0
    fail = 0
    for _, _, status in rows:
        if status.startswith("OK"):
            ok += 1
        else:
            fail += 1
    return ok, fail


def write_markdown(
    rows: list[tuple[str, str, str]],
    out_md: Path,
    tsv_path: Path,
    started_at: str,
    finished_at: str,
) -> None:
    ok, fail = summarize_status(rows)
    overall = "PASS" if fail == 0 else "PASS_WITH_WARNINGS"
    out_md.parent.mkdir(parents=True, exist_ok=True)
    with out_md.open("w", encoding="utf-8") as handle:
        handle.write("# Phase 1 Environment Check\n\n")
        handle.write(f"- Started: {started_at}\n")
        handle.write(f"- Finished: {finished_at}\n")
        handle.write(f"- Overall status: {overall}\n")
        handle.write(f"- OK checks: {ok}\n")
        handle.write(f"- Non-OK checks: {fail}\n")
        handle.write(f"- TSV: `{tsv_path}`\n\n")
        handle.write("| check_type | name | status |\n")
        handle.write("|---|---|---|\n")
        for check_type, name, status in rows:
            safe_status = status.replace("|", "\\|")
            handle.write(f"| {check_type} | {name} | {safe_status} |\n")


def main() -> int:
    script_path = Path(__file__).resolve()
    project_root = script_path.parents[1]
    reports_dir = project_root / "data" / "phase1" / "reports"
    out_tsv = reports_dir / "environment_check.tsv"
    out_md = reports_dir / "environment_check.md"

    started_at = dt.datetime.now().isoformat(timespec="seconds")
    rows: list[tuple[str, str, str]] = []

    rows.append(check_python_runtime())
    rows.append(check_disk_space(project_root))
    rows.append(check_dns("rest.uniprot.org"))
    rows.append(check_internet_head("https://rest.uniprot.org"))

    for cmd in REQUIRED_COMMANDS:
        rows.append(check_command(cmd))

    for pkg in REQUIRED_PYTHON_PACKAGES:
        rows.append(check_python_package(pkg))

    finished_at = dt.datetime.now().isoformat(timespec="seconds")
    write_tsv(rows, out_tsv)
    write_markdown(rows, out_md, out_tsv, started_at, finished_at)

    ok, fail = summarize_status(rows)
    print(f"[phase1_00_check_environment] done: ok={ok}, non_ok={fail}")
    print(f"[phase1_00_check_environment] wrote: {out_tsv}")
    print(f"[phase1_00_check_environment] wrote: {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
