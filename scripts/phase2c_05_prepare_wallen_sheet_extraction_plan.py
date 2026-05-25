#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


TARGETS = [
    ("sample_metadata.tsv", "sample_metadata_candidate", "yes", "yes"),
    ("species_abundance.tsv", "species_abundance_candidate", "yes", "no"),
    ("genus_abundance.tsv", "genus_abundance_candidate", "no", "no"),
    ("gene_family_abundance.tsv", "gene_family_candidate", "yes", "no"),
    ("pathway_abundance.tsv", "pathway_candidate", "yes", "no"),
    ("run_accession_mapping.tsv", "unknown", "yes", "yes"),
    ("clinical_metadata.tsv", "clinical_metadata_candidate", "no", "yes"),
]


def load_inventory(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def propose_sheet(candidates: list[dict[str, str]], target_type: str) -> tuple[str, str]:
    matches = [r for r in candidates if r.get("candidate_table_type") == target_type]
    if not matches:
        return "pending_manual_review", "low"
    return matches[0].get("sheet_name", "pending_manual_review"), "medium"


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    inv = root / "data/phase2c/reports/phase2c_wallen_workbook_sheet_inventory.tsv"
    out_meta = root / "metadata/phase2c_wallen_candidate_sheet_mapping_template.tsv"
    out_status = root / "data/phase2c/reports/phase2c_wallen_sheet_mapping_plan_status.txt"
    out_status.parent.mkdir(parents=True, exist_ok=True)

    if not inv.exists():
        out_status.write_text("PENDING_WORKBOOK_INVENTORY\n", encoding="utf-8")
        print(f"[DONE] wrote: {out_status}")
        return 0

    inventory = load_inventory(inv)
    rows = []
    for fname, ctype, req2d, req3 in TARGETS:
        sheet_name, conf = propose_sheet(inventory, ctype)
        rows.append(
            {
                "standard_target_file": fname,
                "candidate_sheet_name": sheet_name,
                "candidate_table_type": ctype,
                "confidence": conf,
                "manual_review_status": "pending_manual_review",
                "local_output_path_if_approved": f"data/phase2b/processed_tables/wallen_prjna834801/{fname}",
                "required_for_phase2d": req2d,
                "required_for_phase3": req3,
                "notes": "proposed from workbook sheet inventory; requires manual confirmation",
            }
        )

    with out_meta.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(
            h,
            fieldnames=[
                "standard_target_file",
                "candidate_sheet_name",
                "candidate_table_type",
                "confidence",
                "manual_review_status",
                "local_output_path_if_approved",
                "required_for_phase2d",
                "required_for_phase3",
                "notes",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)

    out_status.write_text("SHEET_MAPPING_PLAN_READY\n", encoding="utf-8")
    print(f"[DONE] wrote: {out_meta}")
    print(f"[DONE] wrote: {out_status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
