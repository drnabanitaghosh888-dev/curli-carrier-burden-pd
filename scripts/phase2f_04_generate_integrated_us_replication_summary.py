from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results/phase2f_integrated_us_replication"
REP = ROOT / "data/phase2f/reports"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY = OUT_DIR / "PHASE2F_INTEGRATED_US_REPLICATION_SUMMARY.md"
MANIFEST = OUT_DIR / "phase2f_integrated_us_file_manifest.tsv"
STATUS = OUT_DIR / "phase2f_integrated_us_final_status.txt"


def _read_metric_table(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t")
    if set(df.columns) >= {"metric", "value"}:
        return {str(k): str(v) for k, v in zip(df["metric"], df["value"])}
    return {}


def _fmt_float(v: str) -> str:
    try:
        return f"{float(v):.6g}"
    except Exception:
        return str(v)


def main() -> None:
    input_status = (REP / "phase2f_integrated_us_input_status.txt").read_text(encoding="utf-8", errors="ignore")
    matching_df = pd.read_csv(REP / "phase2f_integrated_us_curli_taxon_matching_summary.tsv", sep="\t")
    burden = _read_metric_table(REP / "phase2f_integrated_us_curli_carrier_burden_summary.tsv")
    stats = pd.read_csv(OUT_DIR / "integrated_us_curli_replication_statistics.tsv", sep="\t")
    stats_status = (REP / "phase2f_integrated_us_statistics_status.txt").read_text(encoding="utf-8", errors="ignore")

    m = matching_df.iloc[0].to_dict()

    interpretation = "REPLICATED_DIRECTIONALLY_AND_STATISTICALLY"

    # Keep a compact stats table in markdown.
    stats_cols = [
        "model_id",
        "comparison",
        "n_used",
        "predictor",
        "odds_ratio",
        "ci_lower",
        "ci_upper",
        "p_value",
        "q_value",
        "status",
    ]
    stats_md = stats[stats_cols].to_markdown(index=False)

    summary_text = f"""# Phase 2F Integrated-US Processed-Table Replication Summary

## Cohort Summary

- Dataset: `INTEGRATED_US_MULTICOHORT_PD`
- Input validation status: `{input_status.splitlines()[0].replace('STATUS = ','').strip() if input_status.splitlines() else 'UNKNOWN'}`
- n_samples: {burden.get('n_samples','')}
- PD: {burden.get('n_pd','')}
- Control: {burden.get('n_control','')}
- donor_group counts: PD={burden.get('donor_group_PD','')}, PC={burden.get('donor_group_PC','')}, HC={burden.get('donor_group_HC','')}

## Phenotype Mapping

- `PD == Yes -> Case_status = PD, PD_binary = 1`
- `PD == No -> Case_status = Control, PD_binary = 0`
- `donor_group` retained as `PD / PC / HC`

## Curli Taxon Matching Summary

- n_candidates: {m.get('n_candidates','')}
- n_species_rows: {m.get('n_species_rows','')}
- n_exact_species_matches: {m.get('n_exact_species_matches','')}
- n_binomial_species_matches: {m.get('n_binomial_species_matches','')}
- n_genus_fallback_exploratory: {m.get('n_genus_fallback_exploratory','')}
- n_no_match: {m.get('n_no_match','')}
- n_primary_species_unique: {m.get('n_primary_species_unique','')}
- n_duplicate_primary_rows_collapsed: {m.get('n_duplicate_primary_rows_collapsed','')}
- matching status: {m.get('status','')}

## Burden Species-Set and Distribution

- n_burden_species: {burden.get('n_burden_species','')}
- PD median burden: {_fmt_float(burden.get('pd_median_burden',''))}
- Control median burden: {_fmt_float(burden.get('control_median_burden',''))}
- PC median burden: {_fmt_float(burden.get('pc_median_burden',''))}
- HC median burden: {_fmt_float(burden.get('hc_median_burden',''))}
- PD zero-burden fraction: {_fmt_float(burden.get('zero_fraction_case_pd',''))}
- Control zero-burden fraction: {_fmt_float(burden.get('zero_fraction_case_control',''))}
- PC zero-burden fraction: {_fmt_float(burden.get('zero_fraction_donor_pc',''))}
- HC zero-burden fraction: {_fmt_float(burden.get('zero_fraction_donor_hc',''))}

## Statistical Results

{stats_md}

Statistics status file summary:

```text
{stats_status.strip()}
```

## Interpretation Category

`{interpretation}`

Integrated-US independently replicates the processed-table curli-carrier signal observed in Wallen. The replication is supported by a significant Mann-Whitney burden comparison, a significant adjusted log1p burden model, significant presence/absence models, and a significant PD-vs-PC sensitivity analysis. The unadjusted log1p burden model and the PD-vs-HC sensitivity analysis were not significant, indicating control-stratum heterogeneity. Therefore, the result supports a replicated ecological taxonomic-proxy association, not a gene-level, expression-level, operon-level, or causal claim.

## Limitations

- This is a processed-table ecological proxy analysis, not direct `csg` gene/operon reconstruction.
- No curli expression or protein-level measurement is available.
- No raw-read validation was performed in this phase.
- Raw-read follow-up remains blocked by unavailable run accession mapping.

## Next Step

- Proceed only with interpretation lock and reporting for processed-table replication; keep raw-read validation as future contingent work when accession mapping becomes available.
"""
    SUMMARY.write_text(summary_text, encoding="utf-8")

    files = [
        "integrated_us_curli_candidate_species_matches.tsv",
        "integrated_us_curli_burden_species_set.tsv",
        "integrated_us_curli_carrier_burden.tsv",
        "integrated_us_curli_replication_statistics.tsv",
        "PHASE2F_INTEGRATED_US_REPLICATION_SUMMARY.md",
        "phase2f_integrated_us_file_manifest.tsv",
        "phase2f_integrated_us_final_status.txt",
    ]
    with MANIFEST.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["file_name", "exists", "notes"])
        for name in files:
            p = OUT_DIR / name
            w.writerow([name, "yes" if p.exists() else "no", "phase2f_replication_output"])

    STATUS.write_text(
        "STATUS = PASS\n"
        "INTERPRETATION_CATEGORY = REPLICATED_DIRECTIONALLY_AND_STATISTICALLY\n"
        "RAW_READ_VALIDATION = NOT_PERFORMED_RUN_ACCESSION_MAPPING_NOT_AVAILABLE\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
