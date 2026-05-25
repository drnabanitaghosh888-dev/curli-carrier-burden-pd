#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path
from typing import Any

import yaml
from Bio import SeqIO


GENES = ["csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build per-gene MAFFT alignments for curli genes.")
    parser.add_argument("--config", default="config/phase1_config.yml")
    parser.add_argument("--input-dir", default="data/phase1/curated")
    parser.add_argument("--outdir", default="data/phase1/alignments")
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--mafft-mode",
        default="--auto",
        choices=["--auto", "--localpair", "--globalpair"],
        help="MAFFT strategy flag.",
    )
    return parser.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Invalid YAML: {path}")
    return data


def count_sequences(path: Path) -> int:
    return sum(1 for _ in SeqIO.parse(str(path), "fasta"))


def parse_fasta_lengths(path: Path) -> list[int]:
    return [len(rec.seq) for rec in SeqIO.parse(str(path), "fasta")]


def write_summary(rows: list[dict[str, str]], path: Path) -> None:
    fields = [
        "gene",
        "input_fasta",
        "output_alignment",
        "status",
        "input_sequences",
        "aligned_sequences",
        "alignment_length",
        "mean_gap_fraction",
        "warning",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def alignment_gap_fraction(path: Path) -> float:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        return 1.0
    total_chars = 0
    total_gaps = 0
    for rec in records:
        s = str(rec.seq)
        total_chars += len(s)
        total_gaps += s.count("-")
    if total_chars == 0:
        return 1.0
    return total_gaps / total_chars


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg = load_yaml(root / args.config)
    input_dir = (root / args.input_dir).resolve()
    outdir = (root / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    threads = int(cfg.get("threads", 4))
    threads = max(1, min(threads, 4))
    mafft = shutil.which("mafft")
    if mafft is None:
        raise RuntimeError("mafft not found in PATH.")

    summary_rows: list[dict[str, str]] = []
    for gene in GENES:
        in_faa = input_dir / f"{gene}.hmm_training.nr95.faa"
        if not in_faa.exists():
            summary_rows.append(
                {
                    "gene": gene,
                    "input_fasta": str(in_faa),
                    "output_alignment": str(outdir / f"{gene}.mafft.aln.faa"),
                    "status": "SKIPPED_MISSING_INPUT",
                    "input_sequences": "0",
                    "aligned_sequences": "0",
                    "alignment_length": "0",
                    "mean_gap_fraction": "NA",
                    "warning": "input FASTA missing",
                }
            )
            continue

        out_aln = outdir / f"{gene}.mafft.aln.faa"
        if out_aln.exists() and not args.force:
            raise FileExistsError(f"Alignment exists for {gene}: {out_aln}. Use --force.")

        input_n = count_sequences(in_faa)
        if input_n < 2:
            summary_rows.append(
                {
                    "gene": gene,
                    "input_fasta": str(in_faa),
                    "output_alignment": str(out_aln),
                    "status": "SKIPPED_TOO_FEW_SEQUENCES",
                    "input_sequences": str(input_n),
                    "aligned_sequences": "0",
                    "alignment_length": "0",
                    "mean_gap_fraction": "NA",
                    "warning": "need at least 2 sequences",
                }
            )
            continue

        cmd = [
            mafft,
            args.mafft_mode,
            "--thread",
            str(threads),
            str(in_faa),
        ]
        try:
            with out_aln.open("w", encoding="utf-8") as fout:
                subprocess.run(cmd, check=True, stdout=fout, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as exc:
            summary_rows.append(
                {
                    "gene": gene,
                    "input_fasta": str(in_faa),
                    "output_alignment": str(out_aln),
                    "status": "FAILED",
                    "input_sequences": str(input_n),
                    "aligned_sequences": "0",
                    "alignment_length": "0",
                    "mean_gap_fraction": "NA",
                    "warning": f"mafft failed: {exc.stderr.strip()}",
                }
            )
            continue

        aligned_n = count_sequences(out_aln)
        lengths = parse_fasta_lengths(out_aln)
        aln_len = max(lengths) if lengths else 0
        gap_frac = alignment_gap_fraction(out_aln)
        warning = ""
        if aln_len < 50:
            warning = "very short alignment"
        if gap_frac > 0.7:
            warning = (warning + "; " if warning else "") + "extreme gap fraction"
        if aligned_n < 10:
            warning = (warning + "; " if warning else "") + "too few sequences"

        summary_rows.append(
            {
                "gene": gene,
                "input_fasta": str(in_faa),
                "output_alignment": str(out_aln),
                "status": "OK",
                "input_sequences": str(input_n),
                "aligned_sequences": str(aligned_n),
                "alignment_length": str(aln_len),
                "mean_gap_fraction": f"{gap_frac:.4f}",
                "warning": warning,
            }
        )

    summary_path = outdir / "alignment_summary.tsv"
    write_summary(summary_rows, summary_path)
    print(f"[DONE] wrote: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
