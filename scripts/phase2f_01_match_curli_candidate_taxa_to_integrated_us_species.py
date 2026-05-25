from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CAND = ROOT / "metadata/curli_candidate_taxa_from_phase1.tsv"
RULES = ROOT / "metadata/phase2d_curli_taxon_matching_rules.tsv"
SP = ROOT / "data/phase2e/processed_tables/integrated_us_multicohort_pd/species_abundance.tsv"
OUT_DIR = ROOT / "results/phase2f_integrated_us_replication"
REP = ROOT / "data/phase2f/reports"
OUT_DIR.mkdir(parents=True, exist_ok=True)
REP.mkdir(parents=True, exist_ok=True)

MATCHES = OUT_DIR / "integrated_us_curli_candidate_species_matches.tsv"
PRIMARY_SET = OUT_DIR / "integrated_us_curli_burden_species_set.tsv"
SUM_TSV = REP / "phase2f_integrated_us_curli_taxon_matching_summary.tsv"
STATUS_TXT = REP / "phase2f_integrated_us_curli_taxon_matching_status.txt"
DIAG_TSV = REP / "phase2f_integrated_us_species_name_diagnostics.tsv"

TAX_PREFIX_RE = re.compile(r"^[kpcofgs]__", flags=re.IGNORECASE)
WS_RE = re.compile(r"\s+")
PARENS_RE = re.compile(r"\([^)]*\)")
STRAIN_TAIL_RE = re.compile(
    r"\b(?:subsp\.?|str\.?|strain|serovar|complex|sp\.?|group)\b.*$",
    flags=re.IGNORECASE,
)
NON_ALPHA_RE = re.compile(r"[^a-z\s]")


def _clean_token(t: str) -> str:
    t = TAX_PREFIX_RE.sub("", str(t or "").strip())
    t = t.replace("_", " ")
    t = PARENS_RE.sub(" ", t)
    t = WS_RE.sub(" ", t).strip()
    return t


def _normalized_text(t: str) -> str:
    t = _clean_token(t).lower()
    t = t.replace("[", " ").replace("]", " ")
    t = WS_RE.sub(" ", t).strip()
    return t


def _to_binomial(norm: str) -> str:
    # Remove trailing decorations before extracting first two words.
    s = STRAIN_TAIL_RE.sub("", norm)
    s = NON_ALPHA_RE.sub(" ", s)
    s = WS_RE.sub(" ", s).strip()
    parts = [p for p in s.split(" ") if p]
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    if len(parts) == 1:
        return parts[0]
    return ""


def parse_metaphlan_clade_name(clade_name: str) -> dict[str, str]:
    clade = str(clade_name or "")
    final_field = clade.split("|")[-1] if "|" in clade else clade
    species_raw = final_field
    if "s__" in clade:
        species_raw = clade.rsplit("s__", 1)[-1]
    species_name = _clean_token(species_raw)
    species_norm = _normalized_text(species_name)
    binom_norm = _to_binomial(species_norm)
    genus = binom_norm.split(" ")[0] if binom_norm else (species_norm.split(" ")[0] if species_norm else "")
    return {
        "metaphlan_species_name": species_name,
        "metaphlan_species_normalized": species_norm,
        "metaphlan_binomial_normalized": binom_norm,
        "metaphlan_genus_normalized": genus,
    }


def parse_candidate_taxon(candidate_name: str) -> dict[str, str]:
    full_norm = _normalized_text(candidate_name)
    binom = _to_binomial(full_norm)
    genus = binom.split(" ")[0] if binom else (full_norm.split(" ")[0] if full_norm else "")
    return {
        "candidate_full_normalized": full_norm,
        "candidate_binomial_normalized": binom,
        "candidate_genus": genus,
    }


def _match_candidate_to_species(cand: dict[str, str], sp_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, str]] = []

    exact_hits = sp_df[sp_df["metaphlan_species_normalized"] == cand["candidate_full_normalized"]]
    if not exact_hits.empty:
        for _, hit in exact_hits.iterrows():
            rows.append(
                {
                    "candidate_taxon_name": cand["candidate_taxon_name"],
                    "candidate_full_normalized": cand["candidate_full_normalized"],
                    "candidate_binomial_normalized": cand["candidate_binomial_normalized"],
                    "candidate_genus": cand["candidate_genus"],
                    "metaphlan_clade_name": hit["clade_name"],
                    "metaphlan_species_name": hit["metaphlan_species_name"],
                    "metaphlan_species_normalized": hit["metaphlan_species_normalized"],
                    "metaphlan_binomial_normalized": hit["metaphlan_binomial_normalized"],
                    "metaphlan_genus_normalized": hit["metaphlan_genus_normalized"],
                    "matching_method": "exact_species_match",
                    "include_for_primary_burden": "yes",
                    "confidence_weight": "high",
                    "notes": "exact normalized species match",
                }
            )
        return pd.DataFrame(rows)

    binom_hits = sp_df[
        (sp_df["metaphlan_binomial_normalized"] != "")
        & (sp_df["metaphlan_binomial_normalized"] == cand["candidate_binomial_normalized"])
    ]
    if not binom_hits.empty:
        for _, hit in binom_hits.iterrows():
            rows.append(
                {
                    "candidate_taxon_name": cand["candidate_taxon_name"],
                    "candidate_full_normalized": cand["candidate_full_normalized"],
                    "candidate_binomial_normalized": cand["candidate_binomial_normalized"],
                    "candidate_genus": cand["candidate_genus"],
                    "metaphlan_clade_name": hit["clade_name"],
                    "metaphlan_species_name": hit["metaphlan_species_name"],
                    "metaphlan_species_normalized": hit["metaphlan_species_normalized"],
                    "metaphlan_binomial_normalized": hit["metaphlan_binomial_normalized"],
                    "metaphlan_genus_normalized": hit["metaphlan_genus_normalized"],
                    "matching_method": "binomial_species_match",
                    "include_for_primary_burden": "yes",
                    "confidence_weight": "medium_high",
                    "notes": "binomial species match after normalization",
                }
            )
        return pd.DataFrame(rows)

    genus_hits = sp_df[(sp_df["metaphlan_genus_normalized"] != "") & (sp_df["metaphlan_genus_normalized"] == cand["candidate_genus"])]
    if not genus_hits.empty:
        for _, hit in genus_hits.iterrows():
            rows.append(
                {
                    "candidate_taxon_name": cand["candidate_taxon_name"],
                    "candidate_full_normalized": cand["candidate_full_normalized"],
                    "candidate_binomial_normalized": cand["candidate_binomial_normalized"],
                    "candidate_genus": cand["candidate_genus"],
                    "metaphlan_clade_name": hit["clade_name"],
                    "metaphlan_species_name": hit["metaphlan_species_name"],
                    "metaphlan_species_normalized": hit["metaphlan_species_normalized"],
                    "metaphlan_binomial_normalized": hit["metaphlan_binomial_normalized"],
                    "metaphlan_genus_normalized": hit["metaphlan_genus_normalized"],
                    "matching_method": "genus_fallback_exploratory",
                    "include_for_primary_burden": "no",
                    "confidence_weight": "low",
                    "notes": "genus-only context match; excluded from primary burden",
                }
            )
        return pd.DataFrame(rows)

    return pd.DataFrame(
        [
            {
                "candidate_taxon_name": cand["candidate_taxon_name"],
                "candidate_full_normalized": cand["candidate_full_normalized"],
                "candidate_binomial_normalized": cand["candidate_binomial_normalized"],
                "candidate_genus": cand["candidate_genus"],
                "metaphlan_clade_name": "",
                "metaphlan_species_name": "",
                "metaphlan_species_normalized": "",
                "metaphlan_binomial_normalized": "",
                "metaphlan_genus_normalized": "",
                "matching_method": "no_match",
                "include_for_primary_burden": "no",
                "confidence_weight": "none",
                "notes": "no species/binomial/genus match",
            }
        ]
    )


def _build_species_table(sp: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c in sp["clade_name"].astype(str):
        parsed = parse_metaphlan_clade_name(c)
        rows.append({"clade_name": c, **parsed})
    return pd.DataFrame(rows)


def _write_diagnostics(sp_df: pd.DataFrame) -> None:
    out = sp_df.rename(
        columns={
            "metaphlan_species_name": "parsed_species_name",
            "metaphlan_species_normalized": "parsed_species_normalized",
            "metaphlan_binomial_normalized": "parsed_binomial_normalized",
            "metaphlan_genus_normalized": "parsed_genus",
        }
    ).copy()
    out["notes"] = "parsed from clade_name"
    out[["clade_name", "parsed_species_name", "parsed_species_normalized", "parsed_binomial_normalized", "parsed_genus", "notes"]].to_csv(
        DIAG_TSV, sep="\t", index=False
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase2F taxon matching scaffold")
    parser.add_argument("--execute", action="store_true", help="Write output tables")
    args = parser.parse_args()

    if not args.execute:
        STATUS_TXT.write_text("STATUS = DRY_RUN_ONLY\n", encoding="utf-8")
        return

    cand = pd.read_csv(CAND, sep="\t")
    _ = pd.read_csv(RULES, sep="\t")
    sp = pd.read_csv(SP, sep="\t")

    sp_df = _build_species_table(sp)
    _write_diagnostics(sp_df)

    ccol = "taxon_name" if "taxon_name" in cand.columns else cand.columns[0]
    match_frames = []
    for name in cand[ccol].astype(str):
        c = parse_candidate_taxon(name)
        c["candidate_taxon_name"] = name
        match_frames.append(_match_candidate_to_species(c, sp_df))

    mdf = pd.concat(match_frames, ignore_index=True) if match_frames else pd.DataFrame()

    primary_rows = mdf[mdf["include_for_primary_burden"] == "yes"].copy()
    n_dup_collapsed = int(len(primary_rows) - primary_rows["metaphlan_clade_name"].nunique()) if len(primary_rows) else 0

    primary = (
        primary_rows.groupby("metaphlan_clade_name", as_index=False)
        .agg(
            metaphlan_species_name=("metaphlan_species_name", "first"),
            metaphlan_species_normalized=("metaphlan_species_normalized", "first"),
            matching_method=("matching_method", lambda x: "|".join(sorted(set(x)))),
            n_candidate_taxa_collapsed=("candidate_taxon_name", "nunique"),
            confidence_weight=("confidence_weight", "first"),
            notes=("notes", lambda x: "collapsed primary matches"),
        )
        .rename(columns={"metaphlan_clade_name": "clade_name"})
    )

    mdf.to_csv(MATCHES, sep="\t", index=False)
    primary.to_csv(PRIMARY_SET, sep="\t", index=False)

    n_exact = int((mdf["matching_method"] == "exact_species_match").sum())
    n_binom = int((mdf["matching_method"] == "binomial_species_match").sum())
    n_genus = int((mdf["matching_method"] == "genus_fallback_exploratory").sum())
    n_no = int((mdf["matching_method"] == "no_match").sum())
    n_primary = int(primary.shape[0])
    status = "PASS" if n_primary > 0 else "NO_PRIMARY_SPECIES_MATCHES"

    with SUM_TSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "n_candidates",
                "n_species_rows",
                "n_exact_species_matches",
                "n_binomial_species_matches",
                "n_genus_fallback_exploratory",
                "n_no_match",
                "n_primary_species_unique",
                "n_duplicate_primary_rows_collapsed",
                "status",
            ]
        )
        w.writerow(
            [
                str(len(cand)),
                str(len(sp_df)),
                str(n_exact),
                str(n_binom),
                str(n_genus),
                str(n_no),
                str(n_primary),
                str(n_dup_collapsed),
                status,
            ]
        )

    status_note = "WARNING_BURDEN_COMPUTATION_MUST_NOT_PROCEED" if status == "NO_PRIMARY_SPECIES_MATCHES" else ""
    STATUS_TXT.write_text(
        f"STATUS = {status}\n"
        f"N_PRIMARY_SPECIES_UNIQUE = {n_primary}\n"
        f"{status_note}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
