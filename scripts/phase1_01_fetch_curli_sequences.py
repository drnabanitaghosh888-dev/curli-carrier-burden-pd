#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import requests
import yaml


GENES = ["csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"]
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


@dataclass
class FetchPaths:
    fasta: Path
    tsv: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch csgA-csgG reference sequences from UniProt (Phase 1)."
    )
    parser.add_argument(
        "--config",
        default="config/phase1_config.yml",
        help="Path to phase1 config YAML.",
    )
    parser.add_argument(
        "--query-config",
        default="config/curli_query_terms.yml",
        help="Path to query terms YAML.",
    )
    parser.add_argument(
        "--outdir",
        default="data/phase1/raw",
        help="Output directory for raw files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned UniProt queries and output paths without downloading.",
    )
    parser.add_argument(
        "--max-sequences",
        type=int,
        default=None,
        help="Override maximum sequences per gene (default from config).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files.",
    )
    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=0.3,
        help="Sleep between requests to avoid hammering API.",
    )
    return parser.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing YAML file: {path}")
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Invalid YAML dict: {path}")
    return data


def build_query(gene: str, include_terms: list[str], scope: list[str], reviewed: bool) -> str:
    reviewed_term = "reviewed:true" if reviewed else "reviewed:false"
    include_chunks = [f'({term})' for term in include_terms if term]
    if not include_chunks:
        include_chunks = [f"(gene:{gene})"]
    include_block = " OR ".join(include_chunks)
    scope_block = " OR ".join([f'(taxonomy_name:"{x}")' for x in scope]) if scope else ""
    query = f"({include_block}) AND ({reviewed_term})"
    if scope_block:
        query += f" AND ({scope_block})"
    return query


def make_paths(outdir: Path, gene: str) -> FetchPaths:
    return FetchPaths(
        fasta=outdir / f"{gene}_uniprot_raw.faa",
        tsv=outdir / f"{gene}_uniprot_raw.tsv",
    )


def ensure_writable(paths: FetchPaths, force: bool) -> None:
    existing = [p for p in [paths.fasta, paths.tsv] if p.exists()]
    if existing and not force:
        joined = ", ".join(str(x) for x in existing)
        raise FileExistsError(
            f"Output file(s) already exist ({joined}). Use --force to overwrite."
        )


def request_uniprot(query: str, size: int, cursor: str | None = None) -> dict[str, Any]:
    params = {
        "query": query,
        "format": "json",
        "size": min(size, 500),
        "fields": ",".join(
            [
                "accession",
                "id",
                "gene_names",
                "protein_name",
                "organism_name",
                "organism_id",
                "lineage",
                "reviewed",
                "length",
                "sequence",
            ]
        ),
    }
    if cursor:
        params["cursor"] = cursor
    resp = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=30)
    resp.raise_for_status()
    return resp.json()


def extract_gene_name(entry: dict[str, Any], fallback: str) -> str:
    genes = entry.get("genes", [])
    for gene in genes:
        if "geneName" in gene and "value" in gene["geneName"]:
            return str(gene["geneName"]["value"])
    return fallback


def extract_protein_name(entry: dict[str, Any]) -> str:
    pd = entry.get("proteinDescription", {})
    rec = pd.get("recommendedName", {})
    full = rec.get("fullName", {})
    if "value" in full:
        return str(full["value"])
    return ""


def extract_lineage(entry: dict[str, Any]) -> str:
    lineage = entry.get("organism", {}).get("lineage", [])
    if isinstance(lineage, list):
        return "; ".join(str(x) for x in lineage)
    return ""


def iterate_results(query: str, max_sequences: int, sleep_seconds: float) -> list[dict[str, Any]]:
    collected: list[dict[str, Any]] = []
    cursor: str | None = None
    while len(collected) < max_sequences:
        batch_size = min(500, max_sequences - len(collected))
        payload = request_uniprot(query=query, size=batch_size, cursor=cursor)
        results = payload.get("results", [])
        if not results:
            break
        if not isinstance(results, list):
            break
        collected.extend(results)
        cursor = payload.get("nextPageCursor")
        if not cursor:
            break
        time.sleep(max(0.0, sleep_seconds))
    return collected[:max_sequences]


def map_entries(
    entries: list[dict[str, Any]],
    gene_symbol_query: str,
    source_query: str,
    reviewed_status: str,
) -> list[dict[str, str]]:
    retrieval_date = dt.date.today().isoformat()
    mapped: list[dict[str, str]] = []
    for entry in entries:
        seq_block = entry.get("sequence", {}) or {}
        sequence = str(seq_block.get("value", "")).strip()
        if not sequence:
            continue
        accession = str(entry.get("primaryAccession", "")).strip()
        seq_id = str(entry.get("uniProtkbId", "")).strip()
        protein_name = extract_protein_name(entry)
        organism = str((entry.get("organism", {}) or {}).get("scientificName", "")).strip()
        taxon_id = str((entry.get("organism", {}) or {}).get("taxonId", "")).strip()
        lineage = extract_lineage(entry)
        length = str(seq_block.get("length", len(sequence)))
        assigned_gene = extract_gene_name(entry, fallback=gene_symbol_query)
        md5_hash = hashlib.md5(sequence.encode("utf-8")).hexdigest()
        mapped.append(
            {
                "sequence_id": seq_id or accession,
                "accession": accession,
                "gene_symbol_query": gene_symbol_query,
                "assigned_gene": assigned_gene,
                "protein_name": protein_name,
                "organism": organism,
                "taxon_id": taxon_id,
                "lineage": lineage,
                "reviewed_status": reviewed_status,
                "sequence_length": str(length),
                "source_database": "UniProtKB",
                "source_url_or_query": source_query,
                "retrieval_date": retrieval_date,
                "md5_sequence_hash": md5_hash,
                "sequence": sequence,
            }
        )
    return mapped


def write_fasta(records: list[dict[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for rec in records:
            header = (
                f">{rec['accession']}|{rec['assigned_gene']}|{rec['organism'].replace(' ', '_')}"
            )
            handle.write(f"{header}\n{rec['sequence']}\n")


def write_gene_tsv(records: list[dict[str, str]], path: Path) -> None:
    fields = [
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
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for rec in records:
            row = {k: rec[k] for k in fields}
            writer.writerow(row)


def write_manifest(manifest_rows: list[dict[str, str]], path: Path) -> None:
    fields = [
        "gene",
        "planned_query_reviewed",
        "planned_query_unreviewed",
        "reviewed_count",
        "unreviewed_count",
        "total_count",
        "fasta_path",
        "tsv_path",
        "status",
        "message",
        "timestamp",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in manifest_rows:
            writer.writerow(row)


def write_errors(errors: list[dict[str, str]], path: Path) -> None:
    fields = ["gene", "query_type", "query", "error", "timestamp"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in errors:
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parents[1]
    cfg = load_yaml(project_root / args.config)
    query_cfg = load_yaml(project_root / args.query_config)
    outdir = (project_root / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    max_sequences = int(args.max_sequences or cfg.get("max_sequences_per_gene", 1000))
    scope = cfg.get("allowed_taxonomic_scope", ["Enterobacterales"])
    if not isinstance(scope, list):
        scope = ["Enterobacterales"]

    manifest_rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []

    for gene in GENES:
        paths = make_paths(outdir, gene)
        include_terms = (query_cfg.get(gene, {}) or {}).get("include_terms", [f"gene:{gene}"])
        if not isinstance(include_terms, list):
            include_terms = [f"gene:{gene}"]

        reviewed_query = build_query(gene, include_terms, scope, reviewed=True)
        unreviewed_query = build_query(gene, include_terms, scope, reviewed=False)

        if args.dry_run:
            print(f"[DRY-RUN] gene={gene}")
            print(f"[DRY-RUN] reviewed_query: {reviewed_query}")
            print(f"[DRY-RUN] unreviewed_query: {unreviewed_query}")
            print(f"[DRY-RUN] outputs: {paths.fasta} | {paths.tsv}")
            manifest_rows.append(
                {
                    "gene": gene,
                    "planned_query_reviewed": reviewed_query,
                    "planned_query_unreviewed": unreviewed_query,
                    "reviewed_count": "0",
                    "unreviewed_count": "0",
                    "total_count": "0",
                    "fasta_path": str(paths.fasta),
                    "tsv_path": str(paths.tsv),
                    "status": "DRY_RUN",
                    "message": "No download executed.",
                    "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
                }
            )
            continue

        try:
            ensure_writable(paths, force=args.force)
            reviewed_entries = iterate_results(
                query=reviewed_query,
                max_sequences=max_sequences,
                sleep_seconds=args.sleep_seconds,
            )
            reviewed_records = map_entries(
                reviewed_entries, gene, reviewed_query, reviewed_status="reviewed"
            )

            remaining = max(0, max_sequences - len(reviewed_records))
            unreviewed_records: list[dict[str, str]] = []
            if remaining > 0:
                unreviewed_entries = iterate_results(
                    query=unreviewed_query,
                    max_sequences=remaining,
                    sleep_seconds=args.sleep_seconds,
                )
                unreviewed_records = map_entries(
                    unreviewed_entries, gene, unreviewed_query, reviewed_status="unreviewed"
                )

            all_records = reviewed_records + unreviewed_records
            write_fasta(all_records, paths.fasta)
            write_gene_tsv(all_records, paths.tsv)

            manifest_rows.append(
                {
                    "gene": gene,
                    "planned_query_reviewed": reviewed_query,
                    "planned_query_unreviewed": unreviewed_query,
                    "reviewed_count": str(len(reviewed_records)),
                    "unreviewed_count": str(len(unreviewed_records)),
                    "total_count": str(len(all_records)),
                    "fasta_path": str(paths.fasta),
                    "tsv_path": str(paths.tsv),
                    "status": "OK",
                    "message": "",
                    "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
                }
            )
            print(
                f"[OK] {gene}: reviewed={len(reviewed_records)} "
                f"unreviewed={len(unreviewed_records)} total={len(all_records)}"
            )
        except Exception as exc:
            message = str(exc)
            manifest_rows.append(
                {
                    "gene": gene,
                    "planned_query_reviewed": reviewed_query,
                    "planned_query_unreviewed": unreviewed_query,
                    "reviewed_count": "0",
                    "unreviewed_count": "0",
                    "total_count": "0",
                    "fasta_path": str(paths.fasta),
                    "tsv_path": str(paths.tsv),
                    "status": "FAILED",
                    "message": message,
                    "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
                }
            )
            errors.append(
                {
                    "gene": gene,
                    "query_type": "reviewed_or_unreviewed",
                    "query": json.dumps(
                        {"reviewed": reviewed_query, "unreviewed": unreviewed_query}
                    ),
                    "error": message,
                    "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
                }
            )
            print(f"[FAILED] {gene}: {message}", file=sys.stderr)

    manifest_path = outdir / "uniprot_download_manifest.tsv"
    write_manifest(manifest_rows, manifest_path)

    if errors:
        error_path = outdir / "uniprot_failed_queries.tsv"
        write_errors(errors, error_path)
        print(f"[WARN] failed query log written: {error_path}")

    print(f"[DONE] manifest written: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
