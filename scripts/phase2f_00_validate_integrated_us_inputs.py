from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data/phase2e/processed_tables/integrated_us_multicohort_pd"
REP = ROOT / "data/phase2f/reports"
REP.mkdir(parents=True, exist_ok=True)

OUT_TSV = REP / "phase2f_integrated_us_input_validation.tsv"
OUT_STATUS = REP / "phase2f_integrated_us_input_status.txt"


def main() -> None:
    checks: list[tuple[str, str, str]] = []
    warnings = 0
    failures = 0

    sm = PROC / "sample_metadata.tsv"
    sp = PROC / "species_abundance.tsv"
    runmap = PROC / "run_accession_mapping.NOT_AVAILABLE_YET.txt"

    for name, p in [("sample_metadata_exists", sm), ("species_abundance_exists", sp)]:
        ok = p.exists()
        checks.append((name, "PASS" if ok else "FAIL", str(p)))
        failures += 0 if ok else 1

    if failures == 0:
        smdf = pd.read_csv(sm, sep="\t")
        spdf = pd.read_csv(sp, sep="\t")

        req_cols = ["sample_name", "Case_status", "PD_binary", "donor_group"]
        for c in req_cols:
            ok = c in smdf.columns
            checks.append((f"sample_metadata_has_{c}", "PASS" if ok else "FAIL", ""))
            failures += 0 if ok else 1

        for cov in ["host_age", "sex", "host_body_mass_index"]:
            ok = cov in smdf.columns
            checks.append((f"sample_metadata_has_{cov}", "PASS" if ok else "WARN", ""))
            warnings += 0 if ok else 1

        first_col_ok = len(spdf.columns) > 0 and spdf.columns[0] == "clade_name"
        checks.append(("species_first_column_is_clade_name", "PASS" if first_col_ok else "FAIL", ""))
        failures += 0 if first_col_ok else 1

        sp_samples = set(spdf.columns[1:])
        md_samples = set(smdf["sample_name"].astype(str)) if "sample_name" in smdf.columns else set()
        match = sp_samples == md_samples
        checks.append(("sample_id_set_match", "PASS" if match else "FAIL", ""))
        failures += 0 if match else 1

        case_vals = set(smdf["Case_status"].dropna().astype(str)) if "Case_status" in smdf.columns else set()
        case_ok = case_vals == {"PD", "Control"}
        checks.append(("case_status_values_pd_control", "PASS" if case_ok else "FAIL", ",".join(sorted(case_vals))))
        failures += 0 if case_ok else 1

        pd_vals = set(smdf["PD_binary"].dropna().astype(str)) if "PD_binary" in smdf.columns else set()
        pd_ok = pd_vals == {"0", "1"}
        checks.append(("pd_binary_values_0_1", "PASS" if pd_ok else "FAIL", ",".join(sorted(pd_vals))))
        failures += 0 if pd_ok else 1

        dg_vals = set(smdf["donor_group"].dropna().astype(str)) if "donor_group" in smdf.columns else set()
        dg_ok = dg_vals == {"PD", "PC", "HC"}
        checks.append(("donor_group_values_pd_pc_hc", "PASS" if dg_ok else "FAIL", ",".join(sorted(dg_vals))))
        failures += 0 if dg_ok else 1

    runmap_ok = runmap.exists()
    checks.append(("run_accession_mapping_not_available", "PASS" if runmap_ok else "WARN", str(runmap)))
    warnings += 0 if runmap_ok else 1

    status = "FAIL" if failures > 0 else ("PASS_WITH_WARNINGS" if warnings > 0 else "PASS")

    with OUT_TSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["check", "result", "notes"])
        w.writerows(checks)

    OUT_STATUS.write_text(
        f"STATUS = {status}\n"
        f"FAILURE_COUNT = {failures}\n"
        f"WARNING_COUNT = {warnings}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
