#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


EXPECTED_GENES = {"csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"}
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Validate tiny artificial curli mock reference.")
    p.add_argument("--fasta", default="data/mock_reference/curli_mock_proteins.faa")
    p.add_argument("--metadata", default="data/mock_reference/curli_mock_metadata.tsv")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def inside_project(root: Path, p: Path) -> bool:
    r = root.resolve()
    t = p.resolve()
    return t == r or r in t.parents


def parse_header(header: str) -> dict[str, str]:
    parts = header.strip().split("|")
    out: dict[str, str] = {}
    for part in parts:
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def parse_fasta(path: Path) -> list[dict[str, str]]:
    records = []
    header = None
    seq = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                meta = parse_header(header[1:])
                records.append({"header": header, "sequence": "".join(seq), **meta})
            header = line
            seq = []
        else:
            seq.append(line.strip())
    if header is not None:
        meta = parse_header(header[1:])
        records.append({"header": header, "sequence": "".join(seq), **meta})
    return records


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    fasta = (root / args.fasta).resolve()
    meta = (root / args.metadata).resolve()
    if not inside_project(root, fasta):
        raise ValueError(f"FASTA path escapes project: {fasta}")
    if not inside_project(root, meta):
        raise ValueError(f"Metadata path escapes project: {meta}")

    if args.dry_run:
        print("[DRY-RUN] validate_mock_curli_reference")
        print(f"[DRY-RUN] fasta={fasta}")
        print(f"[DRY-RUN] metadata={meta}")
        return 0

    if not fasta.exists() or not meta.exists():
        raise FileNotFoundError("Mock FASTA and metadata must both exist.")

    fasta_records = parse_fasta(fasta)
    if len(fasta_records) != 7:
        raise ValueError(f"Expected 7 FASTA records, found {len(fasta_records)}")

    fasta_genes = set()
    fasta_acc = set()
    for rec in fasta_records:
        gene = rec.get("gene_symbol", "")
        acc = rec.get("mock_accession", "")
        seq = rec.get("sequence", "")
        if not gene or not acc:
            raise ValueError("Missing gene_symbol or mock_accession in FASTA header.")
        fasta_genes.add(gene)
        if acc in fasta_acc:
            raise ValueError(f"Duplicate mock_accession in FASTA: {acc}")
        fasta_acc.add(acc)
        if not seq:
            raise ValueError(f"Empty sequence for {gene}")
        bad = set(seq.upper()) - VALID_AA
        if bad:
            raise ValueError(f"Invalid amino acid symbols for {gene}: {sorted(bad)}")

    if fasta_genes != EXPECTED_GENES:
        raise ValueError(f"FASTA gene set mismatch: {sorted(fasta_genes)}")

    with meta.open("r", encoding="utf-8") as h:
        reader = csv.DictReader(h, delimiter="\t")
        rows = [dict(r) for r in reader]

    if len(rows) != 7:
        raise ValueError(f"Expected 7 metadata rows, found {len(rows)}")

    meta_genes = {r["gene_symbol"] for r in rows}
    if meta_genes != EXPECTED_GENES:
        raise ValueError(f"Metadata gene set mismatch: {sorted(meta_genes)}")

    seen_acc = set()
    by_gene = {r["gene_symbol"]: r for r in rows}
    for rec in fasta_records:
        gene = rec["gene_symbol"]
        acc = rec["mock_accession"]
        if gene not in by_gene:
            raise ValueError(f"Gene missing in metadata: {gene}")
        r = by_gene[gene]
        if r["mock_accession"] != acc:
            raise ValueError(f"Accession mismatch for {gene}: fasta={acc}, meta={r['mock_accession']}")
        if acc in seen_acc:
            raise ValueError(f"Duplicate mock_accession in metadata: {acc}")
        seen_acc.add(acc)
        if int(r["sequence_length"]) != len(rec["sequence"]):
            raise ValueError(f"Length mismatch for {gene}")
        if not r["expected_role"].strip():
            raise ValueError(f"Empty expected_role for {gene}")

    print("[OK] mock curli reference validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
