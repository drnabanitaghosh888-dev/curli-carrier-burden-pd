from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/phase2d_03_compute_curli_carrier_burden.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("phase2d_03", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_dedup_collapses_duplicate_candidates_and_retains_max_weight():
    mod = _load_module()
    df = pd.DataFrame(
        [
            {
                "candidate_taxon_name": "Escherichia coli",
                "curli_candidate_confidence": "high",
                "metaphlan_clade_name": "k__Bacteria|...|s__Escherichia_coli",
                "metaphlan_species_name": "Escherichia coli",
                "matching_method": "exact_species_match",
                "include_for_primary_burden": "yes",
            },
            {
                "candidate_taxon_name": "Escherichia coli O157:H7",
                "curli_candidate_confidence": "low",
                "metaphlan_clade_name": "k__Bacteria|...|s__Escherichia_coli",
                "metaphlan_species_name": "Escherichia coli",
                "matching_method": "binomial_species_match",
                "include_for_primary_burden": "yes",
            },
            {
                "candidate_taxon_name": "Salmonella bongori",
                "curli_candidate_confidence": "medium",
                "metaphlan_clade_name": "k__Bacteria|...|s__Salmonella_enterica",
                "metaphlan_species_name": "Salmonella enterica",
                "matching_method": "genus_fallback_exploratory",
                "include_for_primary_burden": "no",
            },
        ]
    )
    species_set, stats = mod.build_primary_species_set(df)
    assert len(species_set) == 1
    row = species_set.iloc[0]
    assert row["metaphlan_species_name"] == "Escherichia coli"
    assert float(row["confidence_weight"]) == 1.0
    assert int(row["n_candidate_taxa_collapsed"]) == 2
    assert "Escherichia coli" in row["collapsed_candidate_taxa"]
    assert stats["genus_fallback_rows_excluded"] == 1


def test_default_mode_is_dry_run_and_execute_flag_present():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    assert "--execute" in text
    assert "if not args.execute" in text
    assert "dry-run only" in text


def test_no_fastq_sra_or_raw_tool_commands_present():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    forbidden = [
        "prefetch",
        "fasterq-dump",
        "fastq",
        "hmmsearch",
        "hmmbuild",
        "kraken",
        "humann",
        "bowtie2",
        "blast",
        "diamond",
        "mafft",
    ]
    for token in forbidden:
        assert token not in text
