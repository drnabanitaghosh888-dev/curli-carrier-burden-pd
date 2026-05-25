from __future__ import annotations

import importlib.util
import py_compile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SCRIPTS = [
    ROOT / "scripts/phase2f_00_validate_integrated_us_inputs.py",
    ROOT / "scripts/phase2f_01_match_curli_candidate_taxa_to_integrated_us_species.py",
    ROOT / "scripts/phase2f_02_compute_integrated_us_curli_carrier_burden.py",
    ROOT / "scripts/phase2f_03_run_integrated_us_replication_statistics.py",
    ROOT / "scripts/phase2f_04_generate_integrated_us_replication_summary.py",
]


def _load_phase2f_matcher_module():
    script = ROOT / "scripts/phase2f_01_match_curli_candidate_taxa_to_integrated_us_species.py"
    spec = importlib.util.spec_from_file_location("phase2f_match", script)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_phase2f_scaffold_files_exist():
    required = SCRIPTS + [
        ROOT / "config/phase2f_integrated_us_config.yml",
        ROOT / "metadata/phase2f_integrated_us_model_plan.tsv",
        ROOT / "docs/PHASE2F_INTEGRATED_US_REPLICATION_ANALYSIS_PLAN.md",
        ROOT / "tests/test_phase2f_integrated_us_replication_scaffold.py",
    ]
    missing = [str(p) for p in required if not p.exists()]
    assert not missing


def test_scripts_compile():
    for p in SCRIPTS:
        py_compile.compile(str(p), doraise=True)


def test_paths_and_primary_secondary_predictors_and_sensitivity_are_declared():
    txt = "\n".join(p.read_text(encoding="utf-8", errors="ignore").lower() for p in SCRIPTS)
    assert "integrated_us_multicohort_pd" in txt
    assert "sample_metadata.tsv" in txt
    assert "species_abundance.tsv" in txt
    assert "curlicarrierburden_log1p" in txt
    assert "curlicarrierpresence" in txt
    assert "pd_vs_pc_only" in txt
    assert "pd_vs_hc_only" in txt


def test_default_modes_are_dry_run():
    stats = (ROOT / "scripts/phase2f_03_run_integrated_us_replication_statistics.py").read_text(encoding="utf-8", errors="ignore").lower()
    burden = (ROOT / "scripts/phase2f_02_compute_integrated_us_curli_carrier_burden.py").read_text(encoding="utf-8", errors="ignore").lower()
    matcher = (ROOT / "scripts/phase2f_01_match_curli_candidate_taxa_to_integrated_us_species.py").read_text(encoding="utf-8", errors="ignore").lower()
    assert "if not args.execute" in stats
    assert "if not args.execute" in burden
    assert "if not args.execute" in matcher
    assert "dry_run_only" in stats


def test_species_normalization_and_matching_rules():
    mod = _load_phase2f_matcher_module()

    parsed = mod.parse_metaphlan_clade_name(
        "k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli"
    )
    assert parsed["metaphlan_species_name"] == "Escherichia coli"
    assert parsed["metaphlan_species_normalized"] == "escherichia coli"
    assert parsed["metaphlan_binomial_normalized"] == "escherichia coli"
    assert parsed["metaphlan_genus_normalized"] == "escherichia"

    cand_exact = mod.parse_candidate_taxon("Escherichia coli")
    assert cand_exact["candidate_full_normalized"] == "escherichia coli"
    assert cand_exact["candidate_binomial_normalized"] == "escherichia coli"

    cand_complex = mod.parse_candidate_taxon("Enterobacter cloacae complex sp. Mu1197")
    assert cand_complex["candidate_binomial_normalized"] == "enterobacter cloacae"


def test_matcher_script_declares_required_methods_and_outputs_and_status_logic():
    t = (ROOT / "scripts/phase2f_01_match_curli_candidate_taxa_to_integrated_us_species.py").read_text(
        encoding="utf-8", errors="ignore"
    ).lower()

    must_have = [
        "exact_species_match",
        "binomial_species_match",
        "genus_fallback_exploratory",
        "include_for_primary_burden",
        "no_primary_species_matches",
        "warning_burden_computation_must_not_proceed",
        "phase2f_integrated_us_species_name_diagnostics.tsv",
        "n_exact_species_matches",
        "n_binomial_species_matches",
        "n_genus_fallback_exploratory",
        "n_primary_species_unique",
        "n_duplicate_primary_rows_collapsed",
    ]
    for token in must_have:
        assert token in t


def test_statistics_script_has_real_execute_logic_and_bh_and_models():
    t = (ROOT / "scripts/phase2f_03_run_integrated_us_replication_statistics.py").read_text(
        encoding="utf-8", errors="ignore"
    ).lower()

    must_have = [
        "mannwhitneyu",
        "smf.logit",
        "_apply_bh_qvalues",
        "curlicarrierburden_log1p",
        "curlicarrierpresence",
        "host_age + c(sex) + host_body_mass_index",
        "pd_vs_pc_only",
        "pd_vs_hc_only",
        "status = dry_run_only",
    ]
    for token in must_have:
        assert token in t

    # Execute mode must explicitly guard against scaffold-only rows.
    assert 'out["status"] == "not_run_in_scaffold"' in t


def test_no_download_rawread_or_phase3_or_heavy_tooling_commands_in_scripts():
    txt = "\n".join(p.read_text(encoding="utf-8", errors="ignore").lower() for p in SCRIPTS)
    forbidden = [
        "wget",
        "prefetch",
        "fasterq-dump",
        "fastq",
        "sra",
        "qiita",
        "phase 3",
        "hmmsearch",
        "hmmbuild",
        "kraken",
        "humann",
        "metaphlan cli",
        "bowtie2",
        "blast",
        "diamond",
        "mafft",
    ]
    for token in forbidden:
        assert token not in txt
