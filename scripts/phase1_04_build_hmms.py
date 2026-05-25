#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


GENES = ["csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build and press curli HMMs.")
    parser.add_argument("--input-dir", default="data/phase1/alignments")
    parser.add_argument("--outdir", default="data/phase1/hmms")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def write_summary(rows: list[dict[str, str]], path: Path) -> None:
    fields = [
        "gene",
        "alignment_path",
        "hmm_path",
        "status",
        "message",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    in_dir = (root / args.input_dir).resolve()
    outdir = (root / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    hmmbuild = shutil.which("hmmbuild")
    hmmpress = shutil.which("hmmpress")
    if hmmbuild is None:
        raise RuntimeError("hmmbuild not found in PATH.")
    if hmmpress is None:
        raise RuntimeError("hmmpress not found in PATH.")

    summary_rows: list[dict[str, str]] = []
    built_hmms: list[Path] = []

    for gene in GENES:
        aln = in_dir / f"{gene}.mafft.aln.faa"
        hmm = outdir / f"{gene}.hmm"

        if not aln.exists():
            summary_rows.append(
                {
                    "gene": gene,
                    "alignment_path": str(aln),
                    "hmm_path": str(hmm),
                    "status": "SKIPPED_MISSING_ALIGNMENT",
                    "message": "alignment file missing",
                }
            )
            continue

        if hmm.exists() and not args.force:
            raise FileExistsError(f"HMM already exists for {gene}: {hmm}. Use --force.")

        cmd = [hmmbuild, str(hmm), str(aln)]
        try:
            proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
            msg = "hmmbuild OK"
            if proc.stderr.strip():
                msg = f"{msg}; stderr={proc.stderr.strip()[:300]}"
            summary_rows.append(
                {
                    "gene": gene,
                    "alignment_path": str(aln),
                    "hmm_path": str(hmm),
                    "status": "OK",
                    "message": msg,
                }
            )
            built_hmms.append(hmm)
        except subprocess.CalledProcessError as exc:
            summary_rows.append(
                {
                    "gene": gene,
                    "alignment_path": str(aln),
                    "hmm_path": str(hmm),
                    "status": "FAILED",
                    "message": f"hmmbuild failed: {exc.stderr.strip()}",
                }
            )

    combined = outdir / "curli_csg_combined.hmm"
    if combined.exists() and not args.force:
        raise FileExistsError(f"Combined HMM exists: {combined}. Use --force.")

    with combined.open("w", encoding="utf-8") as handle:
        for hmm in built_hmms:
            handle.write(hmm.read_text(encoding="utf-8"))

    if built_hmms:
        press_cmd = [hmmpress, str(combined)]
        try:
            subprocess.run(press_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            summary_rows.append(
                {
                    "gene": "COMBINED",
                    "alignment_path": "NA",
                    "hmm_path": str(combined),
                    "status": "FAILED",
                    "message": f"hmmpress failed: {exc.stderr.strip()}",
                }
            )
        else:
            summary_rows.append(
                {
                    "gene": "COMBINED",
                    "alignment_path": "NA",
                    "hmm_path": str(combined),
                    "status": "OK",
                    "message": "combined HMM created and hmmpress completed",
                }
            )
    else:
        summary_rows.append(
            {
                "gene": "COMBINED",
                "alignment_path": "NA",
                "hmm_path": str(combined),
                "status": "SKIPPED",
                "message": "no gene HMMs built; combined HMM not meaningful",
            }
        )

    summary_path = outdir / "hmm_build_summary.tsv"
    write_summary(summary_rows, summary_path)
    print(f"[DONE] wrote: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
