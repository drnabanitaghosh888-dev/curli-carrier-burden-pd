from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data/phase2i/processed_tables/mao_central_china_pd"
META = PROC / "sample_metadata.tsv"
SPECIES = PROC / "species_abundance.tsv"
BURDEN_SPECIES = ROOT / "results/phase2i_mao_replication/mao_curli_burden_species_set.tsv"
OUT_DIR = ROOT / "results/phase2i_mao_replication"
REP = ROOT / "data/phase2i/reports"
OUT_DIR.mkdir(parents=True, exist_ok=True)
REP.mkdir(parents=True, exist_ok=True)

OUT = OUT_DIR / "mao_curli_carrier_burden.tsv"
SUMMARY = REP / "phase2i_mao_curli_carrier_burden_summary.tsv"
STATUS = REP / "phase2i_mao_curli_carrier_burden_status.txt"


def _metric(writer: csv.writer, name: str, value: object) -> None:
    writer.writerow([name, value])


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase2I Mao curli carrier burden")
    parser.add_argument("--execute", action="store_true", help="Write burden output")
    args = parser.parse_args()

    if not args.execute:
        STATUS.write_text("STATUS = DRY_RUN_ONLY\n", encoding="utf-8")
        return

    meta = pd.read_csv(META, sep="\t")
    species = pd.read_csv(SPECIES, sep="\t")
    burden_species = pd.read_csv(BURDEN_SPECIES, sep="\t")

    clade_col = "matched_clade_name" if "matched_clade_name" in burden_species.columns else burden_species.columns[0]
    selected = set(burden_species[clade_col].dropna().astype(str))
    sample_cols = [c for c in species.columns if c != "clade_name"]
    present = species[species["clade_name"].astype(str).isin(selected)].copy()

    abundance = present[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    burden = abundance.sum(axis=0)
    burden_df = pd.DataFrame({"sample_name": burden.index.astype(str), "CurliCarrierBurden": burden.values})

    out = meta.merge(burden_df, on="sample_name", how="left")
    out["CurliCarrierBurden"] = out["CurliCarrierBurden"].fillna(0.0)
    out["CurliCarrierBurden_log1p"] = np.log1p(out["CurliCarrierBurden"].astype(float))
    out["CurliCarrierPresence"] = (out["CurliCarrierBurden"].astype(float) > 0).astype(int)

    first_cols = [
        "sample_name",
        "Case_status",
        "PD_binary",
        "donor_group",
        "CurliCarrierBurden",
        "CurliCarrierBurden_log1p",
        "CurliCarrierPresence",
    ]
    ordered = first_cols + [c for c in out.columns if c not in first_cols]
    out[ordered].to_csv(OUT, sep="\t", index=False)

    burden_series = out["CurliCarrierBurden"].astype(float)
    pd_mask = out["Case_status"].astype(str) == "PD"
    control_mask = out["Case_status"].astype(str) == "Control"

    with SUMMARY.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["metric", "value"])
        _metric(writer, "n_samples", int(out.shape[0]))
        _metric(writer, "n_pd", int(pd_mask.sum()))
        _metric(writer, "n_control", int(control_mask.sum()))
        _metric(writer, "burden_min", float(burden_series.min()))
        _metric(writer, "burden_median", float(burden_series.median()))
        _metric(writer, "burden_mean", float(burden_series.mean()))
        _metric(writer, "burden_max", float(burden_series.max()))
        _metric(writer, "pd_median_burden", float(burden_series[pd_mask].median()))
        _metric(writer, "control_median_burden", float(burden_series[control_mask].median()))
        _metric(writer, "zero_fraction_case_pd", float((burden_series[pd_mask] == 0).mean()))
        _metric(writer, "zero_fraction_case_control", float((burden_series[control_mask] == 0).mean()))
        _metric(writer, "n_burden_species", int(len(selected)))
        _metric(writer, "status", "PASS")

    STATUS.write_text("STATUS = PASS\n", encoding="utf-8")


if __name__ == "__main__":
    main()
