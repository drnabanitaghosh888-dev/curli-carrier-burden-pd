#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


ALL_GENES = {"csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"}
CORE_GENES = {"csgA", "csgB", "csgD", "csgG"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute toy curli completeness/core scores from mock FASTA.")
    p.add_argument("--fasta", default="data/mock_reference/curli_mock_proteins.faa")
    p.add_argument("--out-json", default=None)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def parse_headers_genes(path: Path) -> set[str]:
    genes = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            body = line[1:]
            for part in body.split("|"):
                if part.startswith("gene_symbol="):
                    genes.add(part.split("=", 1)[1])
    return genes


def compute_scores(detected_genes: set[str]) -> dict[str, float]:
    completeness = len(detected_genes & ALL_GENES) / 7.0
    core_score = len(detected_genes & CORE_GENES) / 4.0
    return {
        "Curli_Operon_Completeness": round(completeness, 6),
        "Core_Curli_Score": round(core_score, 6),
    }


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    fasta = (root / args.fasta).resolve()
    if args.out_json:
        out_json = (root / args.out_json).resolve()
    else:
        out_json = root / "data/mock_reference/mock_curli_scores.json"

    if args.dry_run:
        print("[DRY-RUN] score_mock_curli_operon")
        print(f"[DRY-RUN] fasta={fasta}")
        print(f"[DRY-RUN] out_json={out_json}")
        return 0

    if not fasta.exists():
        raise FileNotFoundError(f"FASTA missing: {fasta}")
    genes = parse_headers_genes(fasta)
    scores = compute_scores(genes)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(scores, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(scores, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
