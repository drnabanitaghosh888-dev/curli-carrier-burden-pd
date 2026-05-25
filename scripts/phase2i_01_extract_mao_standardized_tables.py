from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/phase2i/manual_sources/mao_central_china_pd/supplementary_manual_verified"
T1 = SRC / "Table_1.XLSX"
T2 = SRC / "Table_2.XLSX"
T3 = SRC / "Table_3.XLSX"

OUT = ROOT / "data/phase2i/processed_tables/mao_central_china_pd"
REP = ROOT / "data/phase2i/reports"
OUT.mkdir(parents=True, exist_ok=True)
REP.mkdir(parents=True, exist_ok=True)

PLAN = REP / "phase2i_mao_standardized_extraction_plan.tsv"
STATUS = REP / "phase2i_mao_standardized_extraction_status.txt"
SUMMARY = REP / "phase2i_mao_standardized_extraction_summary.tsv"
TAX_REP = REP / "phase2i_mao_species_taxonomy_extraction_report.tsv"


def _safe_col(name: str) -> str:
    return (
        str(name)
        .replace("\xa0", " ")
        .replace("(", "")
        .replace(")", "")
        .replace("&", "and")
        .replace("-", "_")
        .replace(" ", "_")
        .replace("%", "pct")
        .replace("/", "_")
    )


def _sample_cols(cols: list[str]) -> list[str]:
    out = []
    for c in cols:
        s = str(c)
        if s.startswith("PD_") or s.startswith("SP_"):
            out.append(s)
    return out


def _norm_tax(v: object) -> str:
    if pd.isna(v):
        return ""
    return str(v).strip().replace("_", " ")


def _species_clade(genus: str, species: str) -> tuple[str, str]:
    g = _norm_tax(genus)
    s = _norm_tax(species)
    if s and len(s.split()) >= 2:
        return s, "binomial_as_is"
    if g and s:
        return f"{g} {s}", "genus_plus_epithet"
    deepest = s or g
    if deepest:
        return deepest, "unresolved_species_rank"
    return "unresolved_species_rank", "unresolved_species_rank"


def _write_plan() -> None:
    rows = [
        ["sample_metadata.tsv", str(T1), str(OUT / "sample_metadata.tsv"), "execute_required"],
        ["species_abundance.tsv", str(T3), str(OUT / "species_abundance.tsv"), "execute_required"],
        ["genus_abundance.tsv", str(T3), str(OUT / "genus_abundance.tsv"), "execute_required"],
        ["README_or_data_dictionary.txt", str(SRC), str(OUT / "README_or_data_dictionary.txt"), "execute_required"],
        ["license_or_access_terms.txt", str(SRC), str(OUT / "license_or_access_terms.txt"), "execute_required"],
        ["run_accession_mapping.NOT_AVAILABLE_YET.txt", "PRJNA588035_note", str(OUT / "run_accession_mapping.NOT_AVAILABLE_YET.txt"), "execute_required"],
    ]
    with PLAN.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["output_file", "source", "target", "mode"])
        w.writerows(rows)


def _metadata_from_t1_t2(sample_cols: list[str]) -> pd.DataFrame:
    df = pd.DataFrame({"sample_name": sample_cols})
    df["Case_status"] = np.where(df["sample_name"].str.startswith("PD_"), "PD", "Control")
    df["PD_binary"] = np.where(df["Case_status"] == "PD", 1, 0)
    df["donor_group"] = np.where(df["sample_name"].str.startswith("PD_"), "PD", "SP")

    t1 = pd.read_excel(T1)
    t2 = pd.read_excel(T2)

    t1 = t1.rename(columns={c: f"source_Table1_{_safe_col(c)}" for c in t1.columns if str(c) != "Sample"})
    t2 = t2.rename(columns={c: f"source_Table2_{_safe_col(c)}" for c in t2.columns if str(c) != "Sample"})

    if "Sample" in t1.columns:
        t1["Sample"] = t1["Sample"].astype(str)
        df = df.merge(t1, how="left", left_on="sample_name", right_on="Sample").drop(columns=["Sample"], errors="ignore")
    if "Sample" in t2.columns:
        t2["Sample"] = t2["Sample"].astype(str)
        df = df.merge(t2, how="left", left_on="sample_name", right_on="Sample").drop(columns=["Sample"], errors="ignore")

    preferred = [
        "sample_name",
        "Case_status",
        "PD_binary",
        "donor_group",
        "source_Table1_Gender",
        "source_Table1_Age",
        "source_Table1_Age_of_onset",
        "source_Table1_Duration_of_PD",
        "source_Table1_Hoehn_and_Yahr_HandY_stage",
        "source_Table1_UPDRS_III_score",
        "source_Table2_Platform",
        "source_Table2_Data_Volume",
        "source_Table2_Raw_Reads",
        "source_Table2_Clean_ReadsEliminating_host_reads",
    ]
    cols = [c for c in preferred if c in df.columns] + [c for c in df.columns if c not in preferred]
    return df[cols]


def _abundance_from_sheet(df: pd.DataFrame, rank: str) -> tuple[pd.DataFrame, list[str], pd.DataFrame | None]:
    samples = _sample_cols(list(df.columns))
    tax_cols = [c for c in ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"] if c in df.columns]

    if rank == "Species":
        tax_report_rows = []
        clades = []
        statuses = []
        for i, r in df.iterrows():
            cl, st = _species_clade(r.get("Genus", ""), r.get("Species", ""))
            clades.append(cl)
            statuses.append(st)
            tax_report_rows.append(
                {
                    "row_index": int(i),
                    "Kingdom": r.get("Kingdom", np.nan),
                    "Phylum": r.get("Phylum", np.nan),
                    "Class": r.get("Class", np.nan),
                    "Order": r.get("Order", np.nan),
                    "Family": r.get("Family", np.nan),
                    "Genus": r.get("Genus", np.nan),
                    "Species": r.get("Species", np.nan),
                    "clade_name": cl,
                    "species_rank_status": st,
                }
            )
        abd = df[samples].copy()
        abd.insert(0, "clade_name", clades)
        return abd, samples, pd.DataFrame(tax_report_rows)

    # Genus table
    genus = df["Genus"].astype(str).str.strip().replace({"": np.nan}) if "Genus" in df.columns else pd.Series([np.nan] * len(df))
    valid = genus.notna()
    out = df.loc[valid, samples].copy()
    out.insert(0, "clade_name", genus[valid].values)
    return out, samples, None


def _write_docs() -> None:
    (OUT / "README_or_data_dictionary.txt").write_text(
        "Mao et al. 2021 supplementary processed-table extraction\n"
        "- Table_1.XLSX: clinical/sample metadata (primarily PD cohort fields)\n"
        "- Table_2.XLSX: sequencing/sample information\n"
        "- Table_3.XLSX: taxonomic relative abundance by sample\n"
        "- Species sheet is primary source for species_abundance.tsv\n"
        "- Phenotype mapping from sample IDs: PD_* -> PD/1, SP_* -> Control/0\n"
        "- Processed-table taxonomic proxy only; not raw-read validation\n",
        encoding="utf-8",
    )
    (OUT / "license_or_access_terms.txt").write_text(
        "Source context: Frontiers article supplementary processed tables (open-access article context).\n"
        "This Phase 2I extraction uses processed-table supplementary materials only.\n"
        "No raw-read processing or raw-read validation claim is made here.\n",
        encoding="utf-8",
    )
    (OUT / "run_accession_mapping.NOT_AVAILABLE_YET.txt").write_text(
        "PRJNA588035 is a raw sequence repository accession.\n"
        "Run-level accession mapping is not used in Phase 2I processed-table analysis.\n"
        "Raw-read validation remains not performed.\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase2I Mao standardized extraction")
    parser.add_argument("--execute", action="store_true", help="Write standardized tables")
    args = parser.parse_args()

    _write_plan()
    n_valid = sum(1 for p in [T1, T2, T3] if p.exists())

    if not args.execute:
        STATUS.write_text(
            "MODE = DRY_RUN_ONLY\n"
            "STATUS = DRY_RUN_READY\n"
            f"N_VALID_XLSX_SOURCES = {n_valid}\n"
            f"SOURCE_USED = {SRC}\n",
            encoding="utf-8",
        )
        return

    sp = pd.read_excel(T3, sheet_name="Species")
    gen = pd.read_excel(T3, sheet_name="Genus")

    species_abd, sample_cols, tax_rep = _abundance_from_sheet(sp, "Species")
    genus_abd, sample_cols_gen, _ = _abundance_from_sheet(gen, "Genus")
    md = _metadata_from_t1_t2(sample_cols)

    OUT.mkdir(parents=True, exist_ok=True)
    md.to_csv(OUT / "sample_metadata.tsv", sep="\t", index=False)
    species_abd.to_csv(OUT / "species_abundance.tsv", sep="\t", index=False)
    genus_abd.to_csv(OUT / "genus_abundance.tsv", sep="\t", index=False)
    if tax_rep is not None:
        tax_rep.to_csv(TAX_REP, sep="\t", index=False)

    _write_docs()

    n_pd = int((md["Case_status"] == "PD").sum())
    n_ct = int((md["Case_status"] == "Control").sum())
    match = set(md["sample_name"]) == set(sample_cols)

    with SUMMARY.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "n_samples",
            "n_PD",
            "n_Control",
            "n_species_rows",
            "n_genus_rows",
            "sample_id_match_status",
            "n_valid_xlsx_sources",
            "source_used",
        ])
        w.writerow([
            str(len(sample_cols)),
            str(n_pd),
            str(n_ct),
            str(species_abd.shape[0]),
            str(genus_abd.shape[0]),
            "PASS" if match else "FAIL",
            str(n_valid),
            str(SRC),
        ])

    STATUS.write_text(
        "MODE = EXECUTE\n"
        "STATUS = PASS\n"
        f"N_VALID_XLSX_SOURCES = {n_valid}\n"
        f"SOURCE_USED = {SRC}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
