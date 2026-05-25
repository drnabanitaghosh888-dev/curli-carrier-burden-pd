from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/phase2d_02_match_curli_candidate_taxa_to_wallen_species.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("phase2d_02", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_metaphlan_species_parsing_and_escherichia_match():
    mod = _load_module()
    s = mod._parse_metaphlan_clade("k__Bacteria|p__P|g__Escherichia|s__Escherichia_coli")
    assert s["metaphlan_species_raw"] == "s__Escherichia_coli"
    assert s["metaphlan_species_name"] == "Escherichia coli"
    assert s["metaphlan_genus_normalized"] == "escherichia"
    c = mod._normalize_candidate_taxon("Escherichia coli")
    method, include = mod._match_candidate_to_species(c, s)
    assert method == "exact_species_match"
    assert include == "yes"


def test_salmonella_binomial_fallback_and_genus_exploratory():
    mod = _load_module()
    s = mod._parse_metaphlan_clade("k__Bacteria|g__Salmonella|s__Salmonella_enterica")
    c = mod._normalize_candidate_taxon("Salmonella enterica I")
    method, include = mod._match_candidate_to_species(c, s)
    assert method == "binomial_species_match"
    assert include == "yes"

    c2 = mod._normalize_candidate_taxon("Salmonella bongori")
    method2, include2 = mod._match_candidate_to_species(c2, s)
    assert method2 == "genus_fallback_exploratory"
    assert include2 == "no"


def test_confidence_weights_and_safety_guards_present():
    mod = _load_module()
    assert mod._confidence_weight("high") == 1.0
    assert mod._confidence_weight("medium") == 0.5
    assert mod._confidence_weight("low") == 0.25
    assert mod._confidence_weight("unknown") == 0.25
    assert mod._confidence_weight("other") == 0.25

    text = SCRIPT.read_text(encoding="utf-8").lower()
    assert "--execute" in text
    assert "if not args.execute" in text
    forbidden = [
        "prefetch",
        "fasterq-dump",
        "sra",
        "fastq",
        "hmmsearch",
        "hmmbuild",
        "kraken",
        "humann",
        "metaphlan ",
        "bowtie2",
        "blast",
        "diamond",
        "mafft",
    ]
    for token in forbidden:
        assert token not in text
