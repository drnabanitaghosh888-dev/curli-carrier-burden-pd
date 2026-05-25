from __future__ import annotations

import csv
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data/phase2e/processed_tables/integrated_us_multicohort_pd"
REP = ROOT / "data/phase2e/reports"
REP.mkdir(parents=True, exist_ok=True)

VAL_TSV = REP / "phase2e_integrated_us_standardized_table_validation.tsv"
VAL_STATUS = REP / "phase2e_integrated_us_standardized_table_status.txt"


def main() -> None:
    checks: list[tuple[str, str, str]] = []

    sm = PROC / "sample_metadata.tsv"
    sp = PROC / "species_abundance.tsv"
    rd = PROC / "README_or_data_dictionary.txt"
    lc = PROC / "license_or_access_terms.txt"

    checks.append(("sample_metadata_exists", "PASS" if sm.exists() else "FAIL", str(sm)))
    checks.append(("species_abundance_exists", "PASS" if sp.exists() else "FAIL", str(sp)))
    checks.append(("readme_exists", "PASS" if rd.exists() else "FAIL", str(rd)))
    checks.append(("license_exists", "PASS" if lc.exists() else "FAIL", str(lc)))

    status = "PASS"
    if not (sm.exists() and sp.exists() and rd.exists() and lc.exists()):
        status = "FAIL"
    else:
        smdf = pd.read_csv(sm, sep="\t")
        spdf = pd.read_csv(sp, sep="\t")

        req_cols = ["sample_name", "original_sample_id", "Case_status", "PD_binary"]
        for c in req_cols:
            ok = c in smdf.columns
            checks.append((f"sample_metadata_has_{c}", "PASS" if ok else "FAIL", ""))
            if not ok:
                status = "FAIL"

        if "PD_binary" in smdf.columns:
            pd_ok = set(smdf["PD_binary"].dropna().astype(str).unique()).issubset({"0", "1"})
            checks.append(("pd_binary_is_0_1", "PASS" if pd_ok else "FAIL", ""))
            if not pd_ok:
                status = "FAIL"
        if "Case_status" in smdf.columns:
            cs_vals = set(smdf["Case_status"].dropna().astype(str).unique())
            cs_ok = cs_vals == {"PD", "Control"}
            checks.append(("case_status_is_pd_control", "PASS" if cs_ok else "FAIL", ",".join(sorted(cs_vals))))
            if not cs_ok:
                status = "FAIL"

        if {"Case_status", "PD_binary"}.issubset(set(smdf.columns)):
            agree = (
                ((smdf["Case_status"].astype(str) == "PD") & (smdf["PD_binary"].astype(str) == "1"))
                | ((smdf["Case_status"].astype(str) == "Control") & (smdf["PD_binary"].astype(str) == "0"))
            ).all()
            checks.append(("case_status_pd_binary_agree", "PASS" if agree else "FAIL", ""))
            if not agree:
                status = "FAIL"

        if {"donor_group", "PD"}.issubset(set(smdf.columns)):
            dg = smdf["donor_group"].astype(str).str.strip()
            pdv = smdf["PD"].astype(str).str.strip()
            dg_ok = (
                ((dg == "PD") & (pdv == "Yes"))
                | (dg.isin(["PC", "HC"]) & (pdv == "No"))
            ).all()
            checks.append(("donor_group_pd_consistency", "PASS" if dg_ok else "FAIL", ""))
            if not dg_ok:
                status = "FAIL"

        if spdf.columns[0] != "clade_name":
            checks.append(("species_first_col_clade_name", "FAIL", str(spdf.columns[0])))
            status = "FAIL"
        else:
            checks.append(("species_first_col_clade_name", "PASS", ""))

        sp_samples = list(spdf.columns[1:])
        md_samples = list(smdf["sample_name"]) if "sample_name" in smdf.columns else []
        match = set(sp_samples) == set(md_samples)
        checks.append(("sample_names_match", "PASS" if match else "FAIL", ""))
        if not match:
            status = "FAIL"

    with VAL_TSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["check", "result", "notes"])
        w.writerows(checks)

    VAL_STATUS.write_text(f"STATUS = {status}\n", encoding="utf-8")
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
