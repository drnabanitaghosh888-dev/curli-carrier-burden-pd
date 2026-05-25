#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    outdir = root / "analysis/phase2b"
    outdir.mkdir(parents=True, exist_ok=True)

    readme = outdir / "README_CURlI_CARRIER_BURDEN.md"
    formula = outdir / "curli_carrier_burden_formula.md"
    schema = outdir / "example_curli_carrier_burden_input_schema.tsv"

    readme.write_text(
        "# Phase 2B Curli Carrier Burden Scaffold\n\n"
        "This directory contains template artifacts only.\n"
        "No real cohort computation is executed in Phase 2B.\n"
        "The burden estimate is a processed-table proxy and must later be validated\n"
        "with targeted raw-read csg screening.\n",
        encoding="utf-8",
    )

    formula.write_text(
        "# Curli Carrier Burden Formula (Template)\n\n"
        "CurliCarrierBurden_s = sum_i abundance_s,i * CurliCandidateConfidence_i\n\n"
        "Where:\n"
        "- s is sample\n"
        "- i is candidate taxon\n"
        "- abundance_s,i is processed-table abundance\n"
        "- CurliCandidateConfidence_i is confidence weight from Phase 1-derived taxa catalog\n\n"
        "Important: this is a proxy only and not proof of curli expression or operon integrity.\n",
        encoding="utf-8",
    )

    schema.write_text(
        "sample_id\tdataset_id\ttaxon_name\tabundance\tcurli_candidate_confidence\tconfidence_weight\n"
        "SAMPLE_001\tWALLEN_PRJNA834801\tEscherichia coli\t0.012\thigh\t1.0\n",
        encoding="utf-8",
    )

    print(f"[DONE] wrote: {readme}")
    print(f"[DONE] wrote: {formula}")
    print(f"[DONE] wrote: {schema}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
