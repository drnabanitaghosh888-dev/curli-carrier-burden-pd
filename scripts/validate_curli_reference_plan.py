#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import yaml


EXPECTED_GENES = {"csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"}
REQUIRED_COLUMNS = {
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
}
BOOL_LIKE = {"yes", "no", "true", "false", "1", "0"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Validate Phase 2B curli reference plan TSV.")
    p.add_argument("--phase2-config", default="config/phase2_reference_config.yml")
    p.add_argument("--plan", default=None)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as h:
        d = yaml.safe_load(h)
    if not isinstance(d, dict):
        raise ValueError(f"Invalid YAML mapping: {path}")
    return d


def inside_project(root: Path, p: Path) -> bool:
    r = root.resolve()
    t = p.resolve()
    return t == r or r in t.parents


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg = load_yaml(root / args.phase2_config)
    default_plan = cfg.get("phase2_reference", {}).get("reference_plan", {}).get(
        "output_plan_tsv", "data/reference_manifest/curli_reference_plan.tsv"
    )
    plan = (root / (args.plan or default_plan)).resolve()
    if not inside_project(root, plan):
        raise ValueError(f"Plan path escapes project root: {plan}")

    if args.dry_run:
        print("[DRY-RUN] validate_curli_reference_plan")
        print(f"[DRY-RUN] plan={plan}")
        print("[DRY-RUN] validating declared heavy-step fields only; no heavy execution")
        return 0

    if not plan.exists():
        raise FileNotFoundError(f"Plan file missing: {plan}")

    with plan.open("r", encoding="utf-8") as h:
        reader = csv.DictReader(h, delimiter="\t")
        cols = set(reader.fieldnames or [])
        rows = [dict(r) for r in reader]

    missing_cols = sorted(REQUIRED_COLUMNS - cols)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    genes = {r.get("gene_symbol", "").strip() for r in rows}
    if genes != EXPECTED_GENES:
        missing = sorted(EXPECTED_GENES - genes)
        extra = sorted(genes - EXPECTED_GENES)
        raise ValueError(f"Gene set mismatch. missing={missing}, extra={extra}")

    for row in rows:
        for key in ["requires_download", "requires_alignment", "requires_hmmbuild", "heavy_step"]:
            val = str(row.get(key, "")).strip().lower()
            if val not in BOOL_LIKE:
                raise ValueError(f"Invalid boolean-like value for {key}: {val}")
        for pkey in ["expected_output_fasta", "expected_output_alignment", "expected_output_hmm"]:
            rel = str(row.get(pkey, "")).strip()
            p = (root / rel).resolve()
            if not inside_project(root, p):
                raise ValueError(f"Output path escapes project root for {row.get('gene_symbol')}: {p}")

    print("[OK] Phase 2B reference plan validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
