#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml


RELIABLE = {"csgD", "csgE", "csgF", "csgG"}
CAUTIOUS = {"csgA", "csgB", "csgC"}
ALL_GENES = ["csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"]


def infer_column(columns: list[str], candidates: list[str]) -> str | None:
    lower = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    for c in columns:
        cl = c.lower()
        if any(cand in cl for cand in candidates):
            return c
    return None


def infer_gene_from_rows(df: pd.DataFrame) -> pd.Series:
    gene_col = infer_column(list(df.columns), ["assigned_gene", "gene_symbol_query", "gene", "gene_symbol"])
    if gene_col:
        vals = df[gene_col].astype(str)
    else:
        vals = pd.Series([""] * len(df), index=df.index)
    # normalize to csg[ABCDEFG]
    norm = vals.str.extract(r"(csg[ABCDEFG])", expand=False).fillna("").str.lower()
    return norm


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    cfg_path = root / "config/phase2b_config.yml"
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    meta_path = root / cfg["phase1_metadata"]
    out_path = root / "metadata/curli_candidate_taxa_from_phase1.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not meta_path.exists():
        raise FileNotFoundError(f"Phase 1 metadata not found: {meta_path}")

    df = pd.read_csv(meta_path, sep="\t", dtype=str).fillna("")
    org_col = infer_column(list(df.columns), ["organism", "taxon_name", "taxon", "species"])
    if not org_col:
        raise ValueError("Could not infer organism/taxon column from Phase 1 metadata.")
    gene_series = infer_gene_from_rows(df)
    if (gene_series == "").all():
        raise ValueError("Could not infer curli gene symbols from Phase 1 metadata.")

    work = pd.DataFrame(
        {
            "taxon_name": df[org_col].astype(str).str.strip(),
            "gene": gene_series,
        }
    )
    work = work[(work["taxon_name"] != "") & (work["gene"].isin([g.lower() for g in ALL_GENES]))]
    if work.empty:
        raise ValueError("No usable taxa+gene rows found in Phase 1 metadata.")

    rows = []
    for taxon, grp in work.groupby("taxon_name"):
        counts = {g: int((grp["gene"] == g.lower()).sum()) for g in ALL_GENES}
        n_rel = sum(1 for g in RELIABLE if counts[g] > 0)
        n_cau = sum(1 for g in CAUTIOUS if counts[g] > 0)
        if n_rel >= 2 and n_cau >= 1:
            conf = "high"
        elif n_rel >= 1 or n_cau >= 2:
            conf = "medium"
        else:
            conf = "low"
        rows.append(
            {
                "taxon_name": taxon,
                "n_curli_reference_sequences": int(len(grp)),
                "n_csgA": counts["csgA"],
                "n_csgB": counts["csgB"],
                "n_csgC": counts["csgC"],
                "n_csgD": counts["csgD"],
                "n_csgE": counts["csgE"],
                "n_csgF": counts["csgF"],
                "n_csgG": counts["csgG"],
                "has_reliable_core_markers": "yes" if n_rel >= 1 else "no",
                "has_cautious_markers": "yes" if n_cau >= 1 else "no",
                "curli_candidate_confidence": conf,
                "notes": "processed-table triage proxy only; not proof of expression or operon integrity",
            }
        )

    out_df = pd.DataFrame(rows).sort_values(
        by=["curli_candidate_confidence", "n_curli_reference_sequences", "taxon_name"],
        ascending=[True, False, True],
    )
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[DONE] wrote: {out_path}")
    print(f"[COUNT] candidate_taxa={len(out_df)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
