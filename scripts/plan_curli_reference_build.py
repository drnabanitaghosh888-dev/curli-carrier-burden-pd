#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import yaml


EXPECTED_GENES = ["csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"]
PLAN_COLUMNS = [
    "gene_symbol",
    "reference_group",
    "target_taxa",
    "intended_reference_source",
    "expected_output_fasta",
    "expected_output_alignment",
    "expected_output_hmm",
    "requires_download",
    "requires_alignment",
    "requires_hmmbuild",
    "heavy_step",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate deterministic Phase 2B curli reference plan.")
    p.add_argument("--phase2-config", default="config/phase2_reference_config.yml")
    p.add_argument("--output", default=None)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as h:
        d = yaml.safe_load(h)
    if not isinstance(d, dict):
        raise ValueError(f"Invalid YAML mapping: {path}")
    return d


def inside_project(root: Path, target: Path) -> bool:
    r = root.resolve()
    t = target.resolve()
    return t == r or r in t.parents


def build_rows(cfg: dict[str, Any]) -> list[dict[str, str]]:
    p2 = cfg.get("phase2_reference", {})
    rp = p2.get("reference_plan", {})
    outputs = rp.get("outputs", {})
    group = str(rp.get("reference_group", "curli_operon_reference"))
    taxa = str(rp.get("target_taxa", "Escherichia coli; Salmonella enterica; Enterobacteriaceae"))
    source = str(rp.get("intended_reference_source", "planned_reference_source"))
    fasta_dir = str(outputs.get("fasta_dir", "data/reference_build/fasta"))
    aln_dir = str(outputs.get("alignment_dir", "data/reference_build/alignments"))
    hmm_dir = str(outputs.get("hmm_dir", "data/reference_build/hmms"))

    rows: list[dict[str, str]] = []
    for gene in sorted(EXPECTED_GENES):
        rows.append(
            {
                "gene_symbol": gene,
                "reference_group": group,
                "target_taxa": taxa,
                "intended_reference_source": source,
                "expected_output_fasta": f"{fasta_dir}/{gene}.reference.faa",
                "expected_output_alignment": f"{aln_dir}/{gene}.reference.aln.faa",
                "expected_output_hmm": f"{hmm_dir}/{gene}.reference.hmm",
                "requires_download": "yes",
                "requires_alignment": "yes",
                "requires_hmmbuild": "yes",
                "heavy_step": "yes",
            }
        )
    return rows


def write_plan(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=PLAN_COLUMNS, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg = load_yaml(root / args.phase2_config)
    default_out = cfg.get("phase2_reference", {}).get("reference_plan", {}).get(
        "output_plan_tsv", "data/reference_manifest/curli_reference_plan.tsv"
    )
    out = (root / (args.output or default_out)).resolve()
    if not inside_project(root, out):
        raise ValueError(f"Refusing output path outside project: {out}")

    rows = build_rows(cfg)
    if args.dry_run:
        print("[DRY-RUN] Phase 2B reference plan")
        print(f"[DRY-RUN] output_plan={out}")
        print(f"[DRY-RUN] genes={','.join(sorted(EXPECTED_GENES))}")
        print("[DRY-RUN] heavy steps are declared only, not executed")
        return 0

    write_plan(out, rows)
    print(f"[DONE] wrote: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
