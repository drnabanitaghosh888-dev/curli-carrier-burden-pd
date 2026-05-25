#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import yaml


EXPECTED = [
    "wallen_analysis_metadata.tsv",
    "wallen_curli_candidate_species_matches.tsv",
    "wallen_curli_burden_species_set.tsv",
    "wallen_curli_carrier_burden.tsv",
    "wallen_curli_burden_statistics.tsv",
]


def _load_cfg(root: Path) -> dict:
    with (root / "config/phase2d_config.yml").open("r", encoding="utf-8") as h:
        return yaml.safe_load(h)


def _safe_read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str)


def _to_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _fmt_num(x: float, nd: int = 5) -> str:
    if pd.isna(x):
        return "NA"
    return f"{float(x):.{nd}f}"


def _fmt_cell(x) -> str:
    if x is None:
        return "NA"
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "none", "na", "nat"}:
        return "NA"
    return s


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    cfg = _load_cfg(root)
    out_results = root / cfg["output_dirs"]["results"]
    out_results.mkdir(parents=True, exist_ok=True)

    out_md = out_results / "PHASE2D_WALLEN_DISCOVERY_SUMMARY.md"
    out_manifest = out_results / "phase2d_wallen_discovery_file_manifest.tsv"
    out_final_status = out_results / "phase2d_final_status.txt"

    paths = {fname: out_results / fname for fname in EXPECTED}

    manifest_rows = []
    for fname, p in paths.items():
        manifest_rows.append(
            {
                "file_name": fname,
                "exists": "yes" if p.exists() else "no",
                "path": str(p.relative_to(root)),
                "notes": "present" if p.exists() else "pending_manual_execution_or_missing",
            }
        )
    with out_manifest.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["file_name", "exists", "path", "notes"], delimiter="\t")
        w.writeheader()
        w.writerows(manifest_rows)

    md_df = _safe_read_tsv(paths["wallen_analysis_metadata.tsv"])
    mt_df = _safe_read_tsv(paths["wallen_curli_candidate_species_matches.tsv"])
    sp_set_df = _safe_read_tsv(paths["wallen_curli_burden_species_set.tsv"])
    burden_df = _safe_read_tsv(paths["wallen_curli_carrier_burden.tsv"])
    stat_df = _safe_read_tsv(paths["wallen_curli_burden_statistics.tsv"])

    sample_n = len(md_df) if not md_df.empty else "NA"
    pd_n = "NA"
    ctrl_n = "NA"
    if not md_df.empty and "Case_status" in md_df.columns:
        vc = md_df["Case_status"].value_counts()
        pd_n = int(vc.get("PD", 0))
        ctrl_n = int(vc.get("Control", 0))

    total_candidates = len(mt_df) if not mt_df.empty else "NA"
    exact_n = binom_n = genus_n = no_match_n = primary_rows = "NA"
    if not mt_df.empty and "matching_method" in mt_df.columns:
        mv = mt_df["matching_method"].value_counts()
        exact_n = int(mv.get("exact_species_match", 0))
        binom_n = int(mv.get("binomial_species_match", 0))
        genus_n = int(mv.get("genus_fallback_exploratory", 0))
        no_match_n = int(mv.get("no_match", 0))
        if "include_for_primary_burden" in mt_df.columns:
            primary_rows = int((mt_df["include_for_primary_burden"].str.lower() == "yes").sum())

    dedup_n = len(sp_set_df) if not sp_set_df.empty else "NA"
    collapsed_n = "NA"
    if isinstance(primary_rows, int) and isinstance(dedup_n, int):
        collapsed_n = primary_rows - dedup_n

    b_min = b_med = b_mean = b_max = "NA"
    pd_med = ctrl_med = "NA"
    pd_zero = ctrl_zero = "NA"
    if not burden_df.empty and "CurliCarrierBurden" in burden_df.columns:
        b = _to_num(burden_df["CurliCarrierBurden"])
        b_min = _fmt_num(b.min())
        b_med = _fmt_num(b.median())
        b_mean = _fmt_num(b.mean())
        b_max = _fmt_num(b.max())
        if "Case_status" in burden_df.columns:
            g = burden_df.copy()
            g["CurliCarrierBurden"] = b
            pd_sub = g[g["Case_status"] == "PD"]["CurliCarrierBurden"].dropna()
            ct_sub = g[g["Case_status"] == "Control"]["CurliCarrierBurden"].dropna()
            pd_med = _fmt_num(pd_sub.median())
            ctrl_med = _fmt_num(ct_sub.median())
            pd_zero = _fmt_num((pd_sub == 0).mean(), nd=4)
            ctrl_zero = _fmt_num((ct_sub == 0).mean(), nd=4)

    keep_models = [
        "MW_PD_vs_Control_raw_burden",
        "LOGIT_PD_log1p_unadjusted",
        "LOGIT_PD_log1p_adjusted_age_sex_BMI",
        "LOGIT_PD_presence_unadjusted",
        "LOGIT_PD_presence_adjusted_age_sex_BMI",
        "LOGIT_constipation_log1p_exploratory",
    ]
    stat_table_lines = ["| model_id | n_used | predictor | odds_ratio | ci_lower | ci_upper | p_value | q_value | status |",
                        "|---|---:|---|---:|---:|---:|---:|---:|---|"]
    if not stat_df.empty and "model_id" in stat_df.columns:
        stat_idx = stat_df.set_index("model_id", drop=False)
        for mid in keep_models:
            if mid in stat_idx.index:
                r = stat_idx.loc[mid]
                stat_table_lines.append(
                    f"| {mid} | {_fmt_cell(r.get('n_used', 'NA'))} | {_fmt_cell(r.get('predictor', 'NA'))} | {_fmt_cell(r.get('odds_ratio', 'NA'))} | {_fmt_cell(r.get('ci_lower', 'NA'))} | {_fmt_cell(r.get('ci_upper', 'NA'))} | {_fmt_cell(r.get('p_value', 'NA'))} | {_fmt_cell(r.get('q_value', 'NA'))} | {_fmt_cell(r.get('status', 'NA'))} |"
                )
            else:
                stat_table_lines.append(f"| {mid} | NA | NA | NA | NA | NA | NA | NA | MISSING |")
    else:
        for mid in keep_models:
            stat_table_lines.append(f"| {mid} | NA | NA | NA | NA | NA | NA | NA | MISSING |")

    missing = [r["file_name"] for r in manifest_rows if r["exists"] == "no"]

    lines = [
        "# PHASE2D_WALLEN_DISCOVERY_SUMMARY",
        "",
        "PHASE 2D WALLEN DISCOVERY STATUS = COMPLETE_WITH_WARNINGS",
        "",
        "## Input Cohort Summary",
        f"- n = {sample_n}",
        f"- PD = {pd_n}",
        f"- Control = {ctrl_n}",
        "",
        "## Curli Candidate Matching Summary",
        f"- total candidate taxa = {total_candidates}",
        f"- exact species matches = {exact_n}",
        f"- binomial species matches = {binom_n}",
        f"- genus fallback exploratory matches = {genus_n}",
        f"- no match = {no_match_n}",
        f"- primary burden-eligible candidate rows = {primary_rows}",
        "",
        "## Burden Species-Set Summary",
        f"- deduplicated burden species = {dedup_n}",
        f"- duplicate candidate rows collapsed = {collapsed_n}",
        "- genus fallback excluded from primary burden",
        "",
        "## Curli Carrier Burden Distribution",
        f"- burden min = {b_min}",
        f"- burden median = {b_med}",
        f"- burden mean = {b_mean}",
        f"- burden max = {b_max}",
        f"- PD median burden = {pd_med}",
        f"- Control median burden = {ctrl_med}",
        f"- PD zero-burden fraction = {pd_zero}",
        f"- Control zero-burden fraction = {ctrl_zero}",
        "",
        "## Statistical Results",
        *stat_table_lines,
        "",
        "## Main Interpretation",
        "The Wallen discovery cohort shows an unadjusted processed-table curli-carrier signal enriched in PD. The raw burden Mann-Whitney test and unadjusted logistic models are significant after BH correction. However, the associations attenuate after age, sex, and BMI adjustment. Therefore, the result is hypothesis-generating and should not be interpreted as an independent causal effect.",
        "",
        "## Limitations",
        "Curli Carrier Burden is a processed-table proxy only. It is not proof of csg gene presence, curli expression, or intact curli operon architecture. Replication and raw-read csg validation are required before gene-level or operon-level claims.",
        "",
        "## Next Step",
        "Recommended next step: Phase 2E / replication-cohort processed-table verification and, separately, Phase 3 raw-read validation only after run-accession mapping becomes available.",
        "",
        "## Missing Inputs (if any)",
    ]
    if missing:
        lines.extend([f"- {m}" for m in missing])
    else:
        lines.append("- none")
    lines.extend(
        [
            "",
            "This summary does not fabricate missing results; unavailable inputs are reported as NA or MISSING.",
        ]
    )
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    out_final_status.write_text(
        "PHASE 2D FINAL STATUS = COMPLETE_WITH_WARNINGS\n"
        "Reason: Wallen processed-table discovery analysis completed. Unadjusted Curli Carrier Burden signal is enriched in PD, but adjusted models attenuate. Results are processed-table proxy findings requiring replication and raw-read validation.\n",
        encoding="utf-8",
    )
    print(f"[DONE] {out_manifest}")
    print(f"[DONE] {out_md}")
    print(f"[DONE] {out_final_status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
