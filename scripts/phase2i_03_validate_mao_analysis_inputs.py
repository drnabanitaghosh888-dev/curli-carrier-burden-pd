from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data/phase2i/processed_tables/mao_central_china_pd"
REP = ROOT / "data/phase2i/reports"
REP.mkdir(parents=True, exist_ok=True)

META = PROC / "sample_metadata.tsv"
SPECIES = PROC / "species_abundance.tsv"
VALIDATION = REP / "phase2i_mao_analysis_input_validation.tsv"
STATUS = REP / "phase2i_mao_analysis_input_status.txt"


def _add(checks: list[tuple[str, str, str]], name: str, ok: bool, notes: str = "") -> bool:
    checks.append((name, "PASS" if ok else "FAIL", notes))
    return ok


def main() -> None:
    checks: list[tuple[str, str, str]] = []
    ok_all = True

    ok_all &= _add(checks, "sample_metadata_exists", META.exists(), str(META))
    ok_all &= _add(checks, "species_abundance_exists", SPECIES.exists(), str(SPECIES))

    if META.exists() and SPECIES.exists():
        meta = pd.read_csv(META, sep="\t")
        species = pd.read_csv(SPECIES, sep="\t")

        for col in ["sample_name", "Case_status", "PD_binary", "donor_group"]:
            ok_all &= _add(checks, f"sample_metadata_has_{col}", col in meta.columns)

        ok_all &= _add(
            checks,
            "species_first_col_clade_name",
            len(species.columns) > 0 and species.columns[0] == "clade_name",
            species.columns[0] if len(species.columns) else "missing_columns",
        )

        if "sample_name" in meta.columns:
            sample_cols = list(species.columns[1:])
            meta_samples = list(meta["sample_name"].astype(str))
            ok_all &= _add(checks, "sample_ids_match_exactly", sample_cols == meta_samples)
            ok_all &= _add(checks, "no_duplicate_sample_ids", not pd.Series(meta_samples).duplicated().any())

        if "Case_status" in meta.columns:
            n_pd = int((meta["Case_status"].astype(str) == "PD").sum())
            n_control = int((meta["Case_status"].astype(str) == "Control").sum())
            ok_all &= _add(checks, "n_pd_is_39", n_pd == 39, str(n_pd))
            ok_all &= _add(checks, "n_control_is_39", n_control == 39, str(n_control))

        if len(species.columns) > 1:
            values = species.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
            ok_all &= _add(checks, "abundance_numeric", values.notna().all().all())
            ok_all &= _add(checks, "species_table_not_empty", species.shape[0] > 0 and values.shape[1] > 0)

    with VALIDATION.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["check", "result", "notes"])
        writer.writerows(checks)

    status = "PASS" if ok_all else "FAIL"
    STATUS.write_text(f"STATUS = {status}\n", encoding="utf-8")
    if not ok_all:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
