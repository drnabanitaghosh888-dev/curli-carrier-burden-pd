from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ROOT / "metadata/curli_candidate_taxa_from_phase1.tsv"
SPECIES = ROOT / "data/phase2i/processed_tables/mao_central_china_pd/species_abundance.tsv"
RULES = ROOT / "metadata/phase2d_curli_taxon_matching_rules.tsv"
OUT_DIR = ROOT / "results/phase2i_mao_replication"
REP = ROOT / "data/phase2i/reports"
OUT_DIR.mkdir(parents=True, exist_ok=True)
REP.mkdir(parents=True, exist_ok=True)

MATCHES = OUT_DIR / "mao_curli_candidate_species_matches.tsv"
BURDEN_SET = OUT_DIR / "mao_curli_burden_species_set.tsv"
SUMMARY = REP / "phase2i_mao_curli_taxon_matching_summary.tsv"
STATUS = REP / "phase2i_mao_curli_taxon_matching_status.txt"
DIAG = REP / "phase2i_mao_species_name_diagnostics.tsv"

SPACE_RE = re.compile(r"\s+")
PAREN_RE = re.compile(r"\([^)]*\)")
RANK_PREFIX_RE = re.compile(r"^(?:[kpcofgs]__|[kgps]|species)\s+", flags=re.IGNORECASE)
RANK_FIELD_RE = re.compile(r"^[kpcofgs]__", flags=re.IGNORECASE)
NON_WORD_RE = re.compile(r"[^a-z\s]")
TAIL_RE = re.compile(r"\b(?:subsp\.?|strain|str\.?|serovar|complex|sp\.?|group)\b.*$", flags=re.IGNORECASE)


def _strip_rank_prefix(value: str) -> str:
    s = str(value or "").strip()
    if "|" in s:
        s = s.split("|")[-1]
    s = RANK_FIELD_RE.sub("", s)
    s = s.replace("_", " ")
    s = SPACE_RE.sub(" ", s).strip()
    previous = None
    while previous != s:
        previous = s
        s = RANK_PREFIX_RE.sub("", s).strip()
        s = SPACE_RE.sub(" ", s).strip()
    return s


def normalize_taxon_name(value: str) -> str:
    s = _strip_rank_prefix(value)
    s = PAREN_RE.sub(" ", s)
    s = s.replace("[", " ").replace("]", " ")
    s = SPACE_RE.sub(" ", s).strip().lower()
    return s


def binomial_from_normalized(value: str) -> str:
    s = TAIL_RE.sub("", value)
    s = NON_WORD_RE.sub(" ", s)
    s = SPACE_RE.sub(" ", s).strip()
    parts = [p for p in s.split(" ") if p]
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return parts[0] if parts else ""


def genus_from_binomial(value: str) -> str:
    return value.split(" ")[0] if value else ""


def _candidate_rows() -> pd.DataFrame:
    candidates = pd.read_csv(CANDIDATES, sep="\t")
    col = "taxon_name" if "taxon_name" in candidates.columns else candidates.columns[0]
    out = candidates[[col]].rename(columns={col: "candidate_taxon"}).copy()
    out["normalized_candidate_taxon"] = out["candidate_taxon"].map(normalize_taxon_name)
    out["candidate_binomial"] = out["normalized_candidate_taxon"].map(binomial_from_normalized)
    out["candidate_genus"] = out["candidate_binomial"].map(genus_from_binomial)
    return out


def _species_rows() -> pd.DataFrame:
    species = pd.read_csv(SPECIES, sep="\t", usecols=["clade_name"])
    species["normalized_mao_clade_name"] = species["clade_name"].map(normalize_taxon_name)
    species["mao_binomial"] = species["normalized_mao_clade_name"].map(binomial_from_normalized)
    species["mao_genus"] = species["mao_binomial"].map(genus_from_binomial)
    return species


def _match_one(candidate: pd.Series, species: pd.DataFrame) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    exact = species[species["normalized_mao_clade_name"] == candidate["normalized_candidate_taxon"]]
    method = "exact_species_match"
    hits = exact

    if hits.empty:
        method = "binomial_species_match"
        hits = species[
            (species["mao_binomial"] != "")
            & (species["mao_binomial"] == candidate["candidate_binomial"])
        ]

    if hits.empty:
        method = "genus_fallback_exploratory"
        hits = species[
            (species["mao_genus"] != "")
            & (species["mao_genus"] == candidate["candidate_genus"])
        ]

    if hits.empty:
        return [
            {
                "candidate_taxon": candidate["candidate_taxon"],
                "matched_clade_name": "",
                "normalized_candidate_taxon": candidate["normalized_candidate_taxon"],
                "normalized_mao_clade_name": "",
                "matching_method": "no_match",
                "include_for_primary_burden": "no",
                "notes": "no species/binomial/genus match",
            }
        ]

    include = "yes" if method in {"exact_species_match", "binomial_species_match"} else "no"
    notes = "primary burden eligible" if include == "yes" else "genus fallback excluded from primary burden"
    for _, hit in hits.iterrows():
        rows.append(
            {
                "candidate_taxon": candidate["candidate_taxon"],
                "matched_clade_name": hit["clade_name"],
                "normalized_candidate_taxon": candidate["normalized_candidate_taxon"],
                "normalized_mao_clade_name": hit["normalized_mao_clade_name"],
                "matching_method": method,
                "include_for_primary_burden": include,
                "notes": notes,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase2I Mao curli taxon matching")
    parser.add_argument("--execute", action="store_true", help="Write matching outputs")
    args = parser.parse_args()

    if not args.execute:
        STATUS.write_text("STATUS = DRY_RUN_ONLY\n", encoding="utf-8")
        return

    if RULES.exists():
        _ = pd.read_csv(RULES, sep="\t")

    candidates = _candidate_rows()
    species = _species_rows()

    diagnostics = species.rename(
        columns={
            "clade_name": "original_clade_name",
            "mao_binomial": "normalized_mao_binomial",
            "mao_genus": "normalized_mao_genus",
        }
    )
    diagnostics.to_csv(DIAG, sep="\t", index=False)

    rows: list[dict[str, str]] = []
    for _, candidate in candidates.iterrows():
        rows.extend(_match_one(candidate, species))
    matches = pd.DataFrame(rows)

    primary_rows = matches[matches["include_for_primary_burden"] == "yes"].copy()
    duplicate_collapsed = int(len(primary_rows) - primary_rows["matched_clade_name"].nunique()) if len(primary_rows) else 0
    burden_set = (
        primary_rows.groupby("matched_clade_name", as_index=False)
        .agg(
            normalized_mao_clade_name=("normalized_mao_clade_name", "first"),
            matching_method=("matching_method", lambda x: "|".join(sorted(set(x)))),
            n_candidate_taxa_collapsed=("candidate_taxon", "nunique"),
            notes=("notes", lambda _: "deduplicated primary burden species"),
        )
        if len(primary_rows)
        else pd.DataFrame(
            columns=[
                "matched_clade_name",
                "normalized_mao_clade_name",
                "matching_method",
                "n_candidate_taxa_collapsed",
                "notes",
            ]
        )
    )

    matches.to_csv(MATCHES, sep="\t", index=False)
    burden_set.to_csv(BURDEN_SET, sep="\t", index=False)

    n_primary = int(burden_set.shape[0])
    status = "PASS" if n_primary > 0 else "NO_PRIMARY_SPECIES_MATCHES"
    with SUMMARY.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
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
        writer.writerow(
            [
                str(len(candidates)),
                str(len(species)),
                str(int((matches["matching_method"] == "exact_species_match").sum())),
                str(int((matches["matching_method"] == "binomial_species_match").sum())),
                str(int((matches["matching_method"] == "genus_fallback_exploratory").sum())),
                str(int((matches["matching_method"] == "no_match").sum())),
                str(n_primary),
                str(duplicate_collapsed),
                status,
            ]
        )

    STATUS.write_text(
        f"STATUS = {status}\nN_PRIMARY_SPECIES_UNIQUE = {n_primary}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
