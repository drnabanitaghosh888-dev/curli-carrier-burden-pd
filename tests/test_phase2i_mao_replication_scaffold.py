from __future__ import annotations

import importlib.util
import py_compile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VALIDATE = ROOT / "scripts/phase2i_03_validate_mao_analysis_inputs.py"
MATCH = ROOT / "scripts/phase2i_04_match_curli_candidate_taxa_to_mao_species.py"
BURDEN = ROOT / "scripts/phase2i_05_compute_mao_curli_carrier_burden.py"


def _load_matcher():
    spec = importlib.util.spec_from_file_location("phase2i_mao_match", MATCH)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_phase2i_mao_replication_scripts_exist_and_compile():
    for script in [VALIDATE, MATCH, BURDEN]:
        assert script.exists()
        py_compile.compile(str(script), doraise=True)


def test_matching_normalizes_rank_prefixes():
    mod = _load_matcher()
    assert mod.normalize_taxon_name("s  Methanobrevibacter smithii") == "methanobrevibacter smithii"
    assert mod.normalize_taxon_name("s__Escherichia_coli") == "escherichia coli"
    assert mod.normalize_taxon_name("species  Actinomyces   graevenitzii") == "actinomyces graevenitzii"
    assert mod.binomial_from_normalized("escherichia coli o25b h4 st131") == "escherichia coli"


def test_matching_excludes_genus_fallback_from_primary_burden():
    text = MATCH.read_text(encoding="utf-8", errors="ignore").lower()
    assert "genus_fallback_exploratory" in text
    assert "genus fallback excluded from primary burden" in text
    assert "include_for_primary_burden" in text
    assert "exact_species_match" in text
    assert "binomial_species_match" in text


def test_burden_script_creates_log1p_and_presence_and_execute_gate():
    text = BURDEN.read_text(encoding="utf-8", errors="ignore").lower()
    assert "curlicarrierburden_log1p" in text
    assert "np.log1p" in text
    assert "curlicarrierpresence" in text
    assert "if not args.execute" in text
    assert "dry_run_only" in text
    assert "mao_curli_carrier_burden.tsv" in text


def test_matching_execute_gate_and_output_paths_declared():
    text = MATCH.read_text(encoding="utf-8", errors="ignore").lower()
    assert "if not args.execute" in text
    assert "dry_run_only" in text
    assert "mao_curli_candidate_species_matches.tsv" in text
    assert "mao_curli_burden_species_set.tsv" in text
    assert "phase2i_mao_species_name_diagnostics.tsv" in text


def test_no_download_rawread_heavy_or_phase3_command_patterns():
    combined = "\n".join(
        script.read_text(encoding="utf-8", errors="ignore").lower()
        for script in [VALIDATE, MATCH, BURDEN]
    )
    forbidden = [
        "wget ",
        "curl ",
        "prefetch ",
        "fasterq-dump ",
        "hmmsearch ",
        "hmmbuild ",
        "kraken ",
        "humann ",
        "metaphlan ",
        "bowtie2 ",
        "blast ",
        "diamond ",
        "mafft ",
        "python scripts/phase3",
    ]
    for token in forbidden:
        assert token not in combined
