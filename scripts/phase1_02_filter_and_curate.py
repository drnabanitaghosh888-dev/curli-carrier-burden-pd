#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

import yaml
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


GENES = ["csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter and curate curli reference sequences.")
    parser.add_argument("--config", default="config/phase1_config.yml")
    parser.add_argument("--roles-config", default="config/curli_gene_roles.yml")
    parser.add_argument("--raw-dir", default="data/phase1/raw")
    parser.add_argument("--outdir", default="data/phase1/curated")
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--skip-nr95",
        action="store_true",
        help="Skip nr95 clustering even if cd-hit is available.",
    )
    return parser.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Invalid YAML at {path}")
    return data


def safe_lower(x: str) -> str:
    return (x or "").strip().lower()


def is_fragment_like(text: str) -> bool:
    t = safe_lower(text)
    keys = ["fragment", "partial", "truncated"]
    return any(k in t for k in keys)


def fraction_x(seq: str) -> float:
    if not seq:
        return 1.0
    sx = seq.upper()
    return sx.count("X") / len(sx)


def read_metadata(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(r) for r in reader]


def read_fasta(path: Path) -> list[SeqRecord]:
    return list(SeqIO.parse(str(path), "fasta"))


def write_fasta(records: list[SeqRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, str(path), "fasta")


def write_tsv(rows: list[dict[str, str]], fields: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_cdhit_nr95(input_faa: Path, output_faa: Path, identity: float) -> tuple[bool, str]:
    cdhit = shutil.which("cd-hit")
    if cdhit is None:
        return False, "cd-hit not found; nr95 clustering skipped."
    cmd = [
        cdhit,
        "-i",
        str(input_faa),
        "-o",
        str(output_faa),
        "-c",
        f"{identity:.2f}",
        "-n",
        "5",
        "-d",
        "0",
        "-M",
        "0",
        "-T",
        "1",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        msg = f"cd-hit failed; nr95 clustering skipped. stderr={exc.stderr.strip()}"
        return False, msg
    return True, "cd-hit nr95 clustering completed."


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg = load_yaml(root / args.config)
    roles = load_yaml(root / args.roles_config)
    raw_dir = (root / args.raw_dir).resolve()
    outdir = (root / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    max_fraction_x = float(cfg.get("sequence_quality", {}).get("max_fraction_X", 0.05))
    remove_fragments = bool(cfg.get("sequence_quality", {}).get("remove_fragments", True))
    remove_partial = bool(cfg.get("sequence_quality", {}).get("remove_partial", True))
    nr95_identity = float(cfg.get("cluster_identity_for_hmm_training", 0.95))
    min_warn = int(cfg.get("min_sequences_per_gene_warning", 30))

    combined_records: list[SeqRecord] = []
    combined_metadata: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []

    for gene in GENES:
        gene_meta_path = raw_dir / f"{gene}_uniprot_raw.tsv"
        gene_faa_path = raw_dir / f"{gene}_uniprot_raw.faa"
        if not gene_meta_path.exists() or not gene_faa_path.exists():
            summary_rows.append(
                {
                    "gene": gene,
                    "input_sequences": "0",
                    "removed_fragment_partial": "0",
                    "removed_length": "0",
                    "removed_fraction_X": "0",
                    "removed_exact_duplicates": "0",
                    "curated_sequences": "0",
                    "nr95_sequences": "0",
                    "status": "MISSING_INPUT",
                    "warning": "raw FASTA or metadata missing",
                }
            )
            continue

        metadata = read_metadata(gene_meta_path)
        records = read_fasta(gene_faa_path)
        seq_by_acc = {}
        for rec in records:
            acc = rec.id.split("|")[0]
            seq_by_acc[acc] = rec

        length_min, length_max = roles.get(gene, {}).get("expected_length_range_aa", [1, 10_000])
        seen_seq: set[str] = set()
        removed = Counter()

        curated_gene_records: list[SeqRecord] = []
        curated_gene_meta: list[dict[str, str]] = []

        for row in metadata:
            acc = row.get("accession", "")
            rec = seq_by_acc.get(acc)
            if rec is None:
                continue
            seq = str(rec.seq).strip().upper()
            prot = row.get("protein_name", "")
            assigned = row.get("assigned_gene", "")

            if remove_fragments and is_fragment_like(prot):
                removed["fragment_partial"] += 1
                continue
            if remove_partial and is_fragment_like(assigned):
                removed["fragment_partial"] += 1
                continue

            if len(seq) < int(length_min) or len(seq) > int(length_max):
                removed["length"] += 1
                continue

            if fraction_x(seq) > max_fraction_x:
                removed["fraction_x"] += 1
                continue

            if seq in seen_seq:
                removed["exact_dup"] += 1
                continue
            seen_seq.add(seq)

            new_id = f"{acc}|{gene}|{row.get('organism', '').replace(' ', '_')}"
            out_rec = SeqRecord(seq=rec.seq, id=new_id, description="")
            curated_gene_records.append(out_rec)

            out_row = dict(row)
            out_row["assigned_gene"] = gene
            curated_gene_meta.append(out_row)

        curated_faa = outdir / f"{gene}.curated.faa"
        nr95_faa = outdir / f"{gene}.hmm_training.nr95.faa"
        if (curated_faa.exists() or nr95_faa.exists()) and not args.force:
            raise FileExistsError(
                f"{gene} curated outputs exist; use --force to overwrite: {curated_faa}, {nr95_faa}"
            )

        write_fasta(curated_gene_records, curated_faa)

        nr95_count = len(curated_gene_records)
        warning = ""
        if args.skip_nr95:
            write_fasta(curated_gene_records, nr95_faa)
            warning = "nr95 clustering skipped by --skip-nr95"
        else:
            with tempfile.TemporaryDirectory(prefix=f"{gene}_nr95_") as tmpd:
                tmp_out = Path(tmpd) / f"{gene}.nr95.tmp.faa"
                ok, msg = run_cdhit_nr95(curated_faa, tmp_out, nr95_identity)
                if ok and tmp_out.exists():
                    nr95_records = read_fasta(tmp_out)
                    write_fasta(nr95_records, nr95_faa)
                    nr95_count = len(nr95_records)
                    warning = msg
                else:
                    write_fasta(curated_gene_records, nr95_faa)
                    warning = msg

        status = "OK"
        if len(curated_gene_records) < min_warn:
            status = "LOW_COUNT_WARNING"
            warning = (warning + "; " if warning else "") + f"curated count < {min_warn}"

        summary_rows.append(
            {
                "gene": gene,
                "input_sequences": str(len(metadata)),
                "removed_fragment_partial": str(removed["fragment_partial"]),
                "removed_length": str(removed["length"]),
                "removed_fraction_X": str(removed["fraction_x"]),
                "removed_exact_duplicates": str(removed["exact_dup"]),
                "curated_sequences": str(len(curated_gene_records)),
                "nr95_sequences": str(nr95_count),
                "status": status,
                "warning": warning,
            }
        )

        combined_records.extend(curated_gene_records)
        combined_metadata.extend(curated_gene_meta)

    combined_faa = outdir / "curli_reference_proteins.faa"
    combined_tsv = outdir / "curli_reference_metadata.tsv"
    summary_tsv = outdir / "curli_filtering_summary.tsv"

    write_fasta(combined_records, combined_faa)
    meta_fields = [
        "sequence_id",
        "accession",
        "gene_symbol_query",
        "assigned_gene",
        "protein_name",
        "organism",
        "taxon_id",
        "lineage",
        "reviewed_status",
        "sequence_length",
        "source_database",
        "source_url_or_query",
        "retrieval_date",
        "md5_sequence_hash",
    ]
    if combined_metadata:
        # Normalize to declared fields and preserve missing keys as blank.
        norm_rows = [{k: row.get(k, "") for k in meta_fields} for row in combined_metadata]
    else:
        norm_rows = []
    write_tsv(norm_rows, meta_fields, combined_tsv)

    summary_fields = [
        "gene",
        "input_sequences",
        "removed_fragment_partial",
        "removed_length",
        "removed_fraction_X",
        "removed_exact_duplicates",
        "curated_sequences",
        "nr95_sequences",
        "status",
        "warning",
    ]
    write_tsv(summary_rows, summary_fields, summary_tsv)

    print(f"[DONE] wrote: {combined_faa}")
    print(f"[DONE] wrote: {combined_tsv}")
    print(f"[DONE] wrote: {summary_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
