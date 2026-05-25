#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import yaml


def _load_cfg(root: Path) -> dict:
    with (root / "config/phase2d_config.yml").open("r", encoding="utf-8") as h:
        return yaml.safe_load(h)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    cfg = _load_cfg(root)
    inputs = cfg["input_files"]
    sample_id_col = cfg["sample_id_column"]
    case_col = cfg["case_control_column"]
    case_label = cfg["case_label"]
    control_label = cfg["control_label"]

    out_dir = root / cfg["output_dirs"]["reports"]
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "phase2d_wallen_input_validation.tsv"
    out_status = out_dir / "phase2d_input_status.txt"

    rows: list[dict[str, str]] = []
    errors: list[str] = []
    warnings: list[str] = []

    required_files = ["sample_metadata", "species_abundance", "gene_family_abundance", "pathway_abundance", "curli_candidate_taxa"]
    for key in required_files:
        p = root / inputs[key]
        if not p.exists():
            errors.append(f"missing required input: {key} -> {p}")
        rows.append({"check": f"exists:{key}", "status": "PASS" if p.exists() else "FAIL", "detail": str(p)})

    sm_path = root / inputs["sample_metadata"]
    species_path = root / inputs["species_abundance"]
    if sm_path.exists():
        sm = pd.read_csv(sm_path, sep="\t", dtype=str)
        for c in [sample_id_col, case_col]:
            ok = c in sm.columns
            rows.append({"check": f"sample_metadata_column:{c}", "status": "PASS" if ok else "FAIL", "detail": c})
            if not ok:
                errors.append(f"sample metadata missing column: {c}")
        if case_col in sm.columns:
            vals = set(sm[case_col].dropna().astype(str))
            ok = case_label in vals and control_label in vals
            rows.append({"check": "case_control_labels", "status": "PASS" if ok else "FAIL", "detail": f"observed={sorted(vals)}"})
            if not ok:
                errors.append("Case_status does not contain both PD and Control labels")
    else:
        sm = pd.DataFrame()

    if species_path.exists() and not sm.empty and sample_id_col in sm.columns:
        sp = pd.read_csv(species_path, sep="\t", nrows=5)
        sample_cols = [c for c in sp.columns if c != "clade_name"]
        md_samples = set(sm[sample_id_col].astype(str))
        ok = md_samples.issubset(set(sample_cols))
        rows.append({"check": "sample_id_match_species_columns", "status": "PASS" if ok else "FAIL", "detail": f"metadata_n={len(md_samples)} species_cols_n={len(sample_cols)}"})
        if not ok:
            errors.append("sample IDs do not match species abundance columns")

    rows.append(
        {
            "check": "phase3_mapping_requirement",
            "status": "WARNING",
            "detail": "run_accession_mapping is not required for Phase 2D; still required for Phase 3",
        }
    )
    warnings.append("run_accession_mapping unavailable does not block Phase 2D")

    status = "PASS"
    if errors:
        status = "FAIL"
    elif warnings:
        status = "PASS_WITH_WARNINGS"

    with out_tsv.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["check", "status", "detail"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    out_status.write_text(status + "\n", encoding="utf-8")
    print(f"[DONE] {out_tsv}")
    print(f"[DONE] {out_status}")
    return 0 if status != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
