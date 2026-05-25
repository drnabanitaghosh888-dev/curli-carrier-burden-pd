#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import yaml


BOOLEAN_LIKE = {"yes", "no", "true", "false", "1", "0"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate curli gene manifest for Phase 2A.")
    parser.add_argument("--phase2-config", default="config/phase2_reference_config.yml")
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--allow-heavy",
        action="store_true",
        help="Explicit heavy execution flag. Heavy commands remain blocked in this script.",
    )
    return parser.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        obj = yaml.safe_load(handle)
    if not isinstance(obj, dict):
        raise ValueError(f"YAML must be mapping: {path}")
    return obj


def ensure_inside_project(project_root: Path, p: Path) -> None:
    root = project_root.resolve()
    path = p.resolve()
    if root not in path.parents and path != root:
        raise ValueError(f"Path outside project root is not allowed: {p}")


def read_manifest(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        cols = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]
    return cols, rows


def validate_rows(
    rows: list[dict[str, str]],
    expected_genes: list[str],
    required_cols: list[str],
    valid_weights: list[str],
) -> list[str]:
    errors: list[str] = []
    expected = set(expected_genes)
    seen = set()
    for idx, row in enumerate(rows, start=1):
        gene = str(row.get("gene_symbol", "")).strip()
        role = str(row.get("functional_role", "")).strip()
        weight = str(row.get("biological_weight", "")).strip().lower()
        core = str(row.get("core_score", "")).strip().lower()
        operon = str(row.get("operon_completeness", "")).strip().lower()

        if not gene:
            errors.append(f"row {idx}: empty gene_symbol")
            continue
        seen.add(gene)
        if gene not in expected:
            errors.append(f"row {idx}: unexpected gene_symbol={gene}")
        if not role:
            errors.append(f"row {idx}: empty functional_role for {gene}")
        if weight not in set(valid_weights):
            errors.append(f"row {idx}: invalid biological_weight={weight} for {gene}")
        if core not in BOOLEAN_LIKE:
            errors.append(f"row {idx}: invalid core_score={core} for {gene}")
        if operon not in BOOLEAN_LIKE:
            errors.append(f"row {idx}: invalid operon_completeness={operon} for {gene}")

    missing = sorted(expected - seen)
    if missing:
        errors.append(f"missing expected genes: {missing}")
    return errors


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parents[1]
    cfg = load_yaml(project_root / args.phase2_config)
    p2 = cfg.get("phase2_reference", {})
    expected_genes = list(p2.get("expected_genes", []))
    required_cols = list(p2.get("required_columns", []))
    valid_weights = [str(x).lower() for x in p2.get("valid_weights", ["high", "moderate", "low"])]

    manifest_rel = args.manifest or p2.get("paths", {}).get("manifest_tsv")
    if not manifest_rel:
        raise ValueError("Manifest path not set.")
    manifest = (project_root / manifest_rel).resolve()
    ensure_inside_project(project_root, manifest)

    if args.dry_run:
        print("[DRY-RUN] validate_curli_reference_manifest")
        print(f"[DRY-RUN] manifest: {manifest}")
        print(f"[DRY-RUN] expected_genes: {expected_genes}")
        print(f"[DRY-RUN] required_columns: {required_cols}")
        print("[DRY-RUN] heavy tools blocked (no downloads, no HMMER/BLAST/DIAMOND/CD-HIT/MAFFT).")
        return 0

    if not manifest.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest}")
    cols, rows = read_manifest(manifest)

    missing_cols = [c for c in required_cols if c not in cols]
    if missing_cols:
        raise ValueError(f"Manifest missing required columns: {missing_cols}")

    errors = validate_rows(rows, expected_genes, required_cols, valid_weights)
    if errors:
        for err in errors:
            print(f"[ERROR] {err}")
        return 1

    print("[OK] Manifest validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
