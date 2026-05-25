from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data/phase2i/processed_tables/mao_central_china_pd"
REP = ROOT / "data/phase2i/reports"
REP.mkdir(parents=True, exist_ok=True)

VAL = REP / "phase2i_mao_standardized_table_validation.tsv"
STS = REP / "phase2i_mao_standardized_table_status.txt"


def main() -> None:
    checks: list[tuple[str, str, str]] = []
    status = "PASS"

    sm = PROC / "sample_metadata.tsv"
    sp = PROC / "species_abundance.tsv"

    for name, p in [("sample_metadata_exists", sm), ("species_abundance_exists", sp)]:
        ok = p.exists()
        checks.append((name, "PASS" if ok else "FAIL", str(p)))
        if not ok:
            status = "FAIL"

    if status == "PASS":
        smdf = pd.read_csv(sm, sep="\t")
        spdf = pd.read_csv(sp, sep="\t")

        for c in ["sample_name", "Case_status", "PD_binary", "donor_group"]:
            ok = c in smdf.columns
            checks.append((f"sample_metadata_has_{c}", "PASS" if ok else "FAIL", ""))
            if not ok:
                status = "FAIL"

        cset = set(smdf["Case_status"].dropna().astype(str)) if "Case_status" in smdf.columns else set()
        checks.append(("case_status_contains_pd_control", "PASS" if cset == {"PD", "Control"} else "FAIL", ",".join(sorted(cset))))
        if cset != {"PD", "Control"}:
            status = "FAIL"

        bset = set(smdf["PD_binary"].dropna().astype(str)) if "PD_binary" in smdf.columns else set()
        checks.append(("pd_binary_contains_0_1", "PASS" if bset == {"0", "1"} else "FAIL", ",".join(sorted(bset))))
        if bset != {"0", "1"}:
            status = "FAIL"

        first = spdf.columns[0] if len(spdf.columns) else ""
        ok_first = first == "clade_name"
        checks.append(("species_first_col_clade_name", "PASS" if ok_first else "FAIL", first))
        if not ok_first:
            status = "FAIL"

        sp_samples = list(spdf.columns[1:])
        md_samples = list(smdf["sample_name"]) if "sample_name" in smdf.columns else []
        match = set(sp_samples) == set(md_samples)
        checks.append(("sample_id_match", "PASS" if match else "FAIL", ""))
        if not match:
            status = "FAIL"

        dup = smdf["sample_name"].duplicated().any() if "sample_name" in smdf.columns else True
        checks.append(("no_duplicate_sample_ids", "PASS" if not dup else "FAIL", ""))
        if dup:
            status = "FAIL"

        if sp_samples:
            num_ok = spdf[sp_samples].apply(pd.to_numeric, errors="coerce").notna().all().all()
        else:
            num_ok = False
        checks.append(("abundance_values_numeric", "PASS" if num_ok else "FAIL", ""))
        if not num_ok:
            status = "FAIL"

        nonempty = (len(spdf) > 0) and (len(sp_samples) > 0)
        checks.append(("species_table_not_empty", "PASS" if nonempty else "FAIL", ""))
        if not nonempty:
            status = "FAIL"

        n_pd = int((smdf["Case_status"] == "PD").sum()) if "Case_status" in smdf.columns else -1
        n_ct = int((smdf["Case_status"] == "Control").sum()) if "Case_status" in smdf.columns else -1
        checks.append(("n_pd_samples", "PASS", str(n_pd)))
        checks.append(("n_control_samples", "PASS", str(n_ct)))

    with VAL.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["check", "result", "notes"])
        w.writerows(checks)

    STS.write_text(f"STATUS = {status}\n", encoding="utf-8")
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
