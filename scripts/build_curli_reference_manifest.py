#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build deterministic curli gene manifest.")
    parser.add_argument("--genes-config", default="config/curli_genes.yml")
    parser.add_argument("--phase2-config", default="config/phase2_reference_config.yml")
    parser.add_argument("--output", default=None, help="Override output TSV path.")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--allow-heavy",
        action="store_true",
        help="Explicit heavy execution flag. Not required for this lightweight script.",
    )
    return parser.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        obj = yaml.safe_load(handle)
    if not isinstance(obj, dict):
        raise ValueError(f"YAML must be a mapping: {path}")
    return obj


def ensure_inside_project(project_root: Path, output_path: Path) -> None:
    root = project_root.resolve()
    out = output_path.resolve()
    if root not in out.parents and out != root:
        raise ValueError(f"Refusing to write outside project root: {output_path}")


def as_bool_like(value: Any) -> str:
    return "yes" if bool(value) else "no"


def build_rows(genes_cfg: dict[str, Any], expected_genes: list[str]) -> list[dict[str, str]]:
    genes = genes_cfg.get("curli_genes", {})
    if not isinstance(genes, dict):
        raise ValueError("config/curli_genes.yml missing mapping: curli_genes")
    rows: list[dict[str, str]] = []
    for gene in sorted(expected_genes):
        meta = genes.get(gene, {})
        if not isinstance(meta, dict):
            raise ValueError(f"Gene metadata must be mapping: {gene}")
        rows.append(
            {
                "gene_symbol": gene,
                "functional_role": str(meta.get("role", "")).strip(),
                "biological_weight": str(meta.get("weight", "")).strip().lower(),
                "core_score": as_bool_like(meta.get("core_score", False)),
                "operon_completeness": as_bool_like(meta.get("operon_completeness", False)),
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parents[1]
    genes_cfg = load_yaml(project_root / args.genes_config)
    phase2_cfg = load_yaml(project_root / args.phase2_config)

    p2 = phase2_cfg.get("phase2_reference", {})
    expected_genes = list(p2.get("expected_genes", []))
    required_columns = list(p2.get("required_columns", []))
    output_rel = args.output or p2.get("paths", {}).get("manifest_tsv")
    if not output_rel:
        raise ValueError("No output path configured.")
    output_path = (project_root / output_rel).resolve()
    ensure_inside_project(project_root, output_path)

    rows = build_rows(genes_cfg, expected_genes)
    if args.dry_run:
        print("[DRY-RUN] build_curli_reference_manifest")
        print(f"[DRY-RUN] output: {output_path}")
        print(f"[DRY-RUN] rows: {len(rows)}")
        print(f"[DRY-RUN] columns: {required_columns}")
        return 0

    write_tsv(output_path, rows, required_columns)
    print(f"[DONE] wrote: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
