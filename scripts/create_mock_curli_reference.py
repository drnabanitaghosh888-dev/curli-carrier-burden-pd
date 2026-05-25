#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


MOCK_ROWS = [
    ("csgA", "Mockus enterica", "MOCK_CSGA_001", "major curli fiber subunit"),
    ("csgB", "Mockus enterica", "MOCK_CSGB_001", "nucleator protein for CsgA polymerization"),
    ("csgC", "Mockus enterica", "MOCK_CSGC_001", "periplasmic chaperone/anti-aggregation factor"),
    ("csgD", "Mockus enterica", "MOCK_CSGD_001", "transcriptional regulator of curli and biofilm genes"),
    ("csgE", "Mockus enterica", "MOCK_CSGE_001", "secretion/assembly accessory protein"),
    ("csgF", "Mockus enterica", "MOCK_CSGF_001", "curli assembly accessory protein and CsgB localization factor"),
    ("csgG", "Mockus enterica", "MOCK_CSGG_001", "outer membrane secretion pore for curli subunits"),
]


MOCK_SEQUENCES = {
    "csgA": "MNNATNQATGNNQATNNQATNNQATNNQATNNQATNNQATNNQATNNQATNNQAT",
    "csgB": "MNKQATNNQATNNQATNNQATNNQATNNQATNNQATNNQATNNQATNNQATNNQA",
    "csgC": "MKKLAVAVALLAAGTANAQGQGQGQGQGQGQGQGQGQGQGQGQGQGQGQGQGQGQ",
    "csgD": "MSEKRFLLKQNNNQQQQQQAAATTTGGGSSSNNNQQQAAATTTGGGSSSNNNQQQAA",
    "csgE": "MQQQQTTTAAAGGGNNNQQQTTTAAAGGGNNNQQQTTTAAAGGGNNNQQQTTTAAA",
    "csgF": "MNNQQQAAATTTGGGNNNQQQAAATTTGGGNNNQQQAAATTTGGGNNNQQQAAATT",
    "csgG": "MKKLLPLLLAAALAACSQASNNNQQQAAATTTGGGNNNQQQAAATTTGGGNNNQQQA",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Create tiny artificial curli mock reference files.")
    p.add_argument("--fasta", default="data/mock_reference/curli_mock_proteins.faa")
    p.add_argument("--metadata", default="data/mock_reference/curli_mock_metadata.tsv")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def inside_project(root: Path, p: Path) -> bool:
    r = root.resolve()
    t = p.resolve()
    return t == r or r in t.parents


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
        print("[DRY-RUN] create_mock_curli_reference")
        print(f"[DRY-RUN] fasta={fasta}")
        print(f"[DRY-RUN] metadata={meta}")
        return 0

    fasta.parent.mkdir(parents=True, exist_ok=True)
    meta.parent.mkdir(parents=True, exist_ok=True)

    with fasta.open("w", encoding="utf-8") as fh:
        for gene, species, accession, _role in MOCK_ROWS:
            seq = MOCK_SEQUENCES[gene]
            header = f">gene_symbol={gene}|mock_species={species.replace(' ', '_')}|mock_accession={accession}"
            fh.write(f"{header}\n{seq}\n")

    with meta.open("w", encoding="utf-8", newline="") as mh:
        writer = csv.DictWriter(
            mh,
            fieldnames=[
                "gene_symbol",
                "mock_species",
                "mock_accession",
                "sequence_length",
                "expected_role",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for gene, species, accession, role in MOCK_ROWS:
            writer.writerow(
                {
                    "gene_symbol": gene,
                    "mock_species": species,
                    "mock_accession": accession,
                    "sequence_length": str(len(MOCK_SEQUENCES[gene])),
                    "expected_role": role,
                }
            )

    print(f"[DONE] wrote: {fasta}")
    print(f"[DONE] wrote: {meta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
