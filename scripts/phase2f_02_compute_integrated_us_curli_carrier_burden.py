from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data/phase2e/processed_tables/integrated_us_multicohort_pd"
IN_META = PROC / "sample_metadata.tsv"
IN_AB = PROC / "species_abundance.tsv"
IN_SET = ROOT / "results/phase2f_integrated_us_replication/integrated_us_curli_burden_species_set.tsv"
OUT_DIR = ROOT / "results/phase2f_integrated_us_replication"
REP = ROOT / "data/phase2f/reports"
OUT_DIR.mkdir(parents=True, exist_ok=True)
REP.mkdir(parents=True, exist_ok=True)

OUT = OUT_DIR / "integrated_us_curli_carrier_burden.tsv"
SUM_TSV = REP / "phase2f_integrated_us_curli_carrier_burden_summary.tsv"
STATUS_TXT = REP / "phase2f_integrated_us_curli_carrier_burden_status.txt"


def _median(s: pd.Series) -> float:
    if len(s) == 0:
        return float("nan")
    return float(np.median(s.astype(float)))


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase2F burden computation scaffold")
    parser.add_argument("--execute", action="store_true", help="Write final burden table")
    args = parser.parse_args()

    if not args.execute:
        STATUS_TXT.write_text("STATUS = DRY_RUN_ONLY\n", encoding="utf-8")
        return

    md = pd.read_csv(IN_META, sep="\t")
    ab = pd.read_csv(IN_AB, sep="\t")
    burden_set = pd.read_csv(IN_SET, sep="\t")

    use = set(burden_set["clade_name"].astype(str).tolist())
    ab_sub = ab[ab["clade_name"].astype(str).isin(use)]

    sample_cols = [c for c in ab.columns if c != "clade_name"]
    burden = ab_sub[sample_cols].sum(axis=0)
    bdf = pd.DataFrame({"sample_name": burden.index, "CurliCarrierBurden": burden.values})
    out = md.merge(bdf, on="sample_name", how="left")
    out["CurliCarrierBurden"] = out["CurliCarrierBurden"].fillna(0.0)
    out["CurliCarrierBurden_log1p"] = np.log1p(out["CurliCarrierBurden"].astype(float))
    out["CurliCarrierPresence"] = (out["CurliCarrierBurden"].astype(float) > 0).astype(int)

    cols = [
        "sample_name",
        "Case_status",
        "PD_binary",
        "donor_group",
        "host_age",
        "sex",
        "host_body_mass_index",
        "CurliCarrierBurden",
        "CurliCarrierBurden_log1p",
        "CurliCarrierPresence",
    ]
    for c in cols:
        if c not in out.columns:
            out[c] = np.nan
    out = out[cols]
    out.to_csv(OUT, sep="\t", index=False)

    pd_mask = out["Case_status"].astype(str) == "PD"
    ct_mask = out["Case_status"].astype(str) == "Control"

    dg = out["donor_group"].astype(str)
    burden_s = out["CurliCarrierBurden"].astype(float)

    with SUM_TSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["n_samples", str(len(out))])
        w.writerow(["n_pd", str(int(pd_mask.sum()))])
        w.writerow(["n_control", str(int(ct_mask.sum()))])
        w.writerow(["donor_group_PD", str(int((dg == "PD").sum()))])
        w.writerow(["donor_group_PC", str(int((dg == "PC").sum()))])
        w.writerow(["donor_group_HC", str(int((dg == "HC").sum()))])
        w.writerow(["burden_min", str(float(burden_s.min()))])
        w.writerow(["burden_median", str(_median(burden_s))])
        w.writerow(["burden_mean", str(float(burden_s.mean()))])
        w.writerow(["burden_max", str(float(burden_s.max()))])
        w.writerow(["pd_median_burden", str(_median(burden_s[pd_mask]))])
        w.writerow(["control_median_burden", str(_median(burden_s[ct_mask]))])
        w.writerow(["pc_median_burden", str(_median(burden_s[dg == "PC"]))])
        w.writerow(["hc_median_burden", str(_median(burden_s[dg == "HC"]))])
        w.writerow(["zero_fraction_case_pd", str(float((burden_s[pd_mask] == 0).mean() if pd_mask.any() else np.nan))])
        w.writerow(["zero_fraction_case_control", str(float((burden_s[ct_mask] == 0).mean() if ct_mask.any() else np.nan))])
        w.writerow(["zero_fraction_donor_pd", str(float((burden_s[dg == "PD"] == 0).mean() if (dg == "PD").any() else np.nan))])
        w.writerow(["zero_fraction_donor_pc", str(float((burden_s[dg == "PC"] == 0).mean() if (dg == "PC").any() else np.nan))])
        w.writerow(["zero_fraction_donor_hc", str(float((burden_s[dg == "HC"] == 0).mean() if (dg == "HC").any() else np.nan))])
        w.writerow(["n_burden_species", str(len(use))])

    STATUS_TXT.write_text("STATUS = PASS\n", encoding="utf-8")


if __name__ == "__main__":
    main()
