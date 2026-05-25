from __future__ import annotations

import py_compile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_scripts_exist_and_compile():
    s1 = ROOT / "scripts/phase2e_05_extract_integrated_us_standardized_tables.py"
    s2 = ROOT / "scripts/phase2e_06_validate_integrated_us_standardized_tables.py"
    assert s1.exists()
    assert s2.exists()
    py_compile.compile(str(s1), doraise=True)
    py_compile.compile(str(s2), doraise=True)


def test_extraction_script_has_execute_gate_and_rscript_override_and_phenotype_mapping():
    s1 = ROOT / "scripts/phase2e_05_extract_integrated_us_standardized_tables.py"
    t = s1.read_text(encoding="utf-8").lower()

    must_have = [
        "--execute",
        "phase2e_rscript",
        "which -a rscript",
        "command -v rscript",
        "if not args.execute",
        "status=\"dry_run_only\"",
        "species_abundance.tsv",
        "mandatory_files",
        "missing mandatory files in stage",
        "phylo_objects/species",
        "phase2e_integrated_us_r_extraction_log.txt",
        "phase2e_integrated_us_phenotype_mapping_report.tsv",
        "\"yes\": \"pd\", \"no\": \"control\"",
        "\"yes\": 1, \"no\": 0",
        "if \"original_sample_id\" not in md.columns",
        "if \"id\" in md.columns",
        "elif \"sample_name\" in md.columns",
        "cannot create original_sample_id: neither id nor sample_name is available",
        "required_for_mapping = [\"pd\", \"donor_group\"]",
        "donor_group",
        "study_group",
        "donor_group/pd consistency conflicts",
        "os.replace",
    ]
    for m in must_have:
        assert m in t


def test_r_helper_uses_phylo_objects_species_and_extracts_otu_tax_and_aggregates_duplicates():
    rf = ROOT / "scripts/_tmp_extract_integrated_us_species_phyloseq.R"
    t = rf.read_text(encoding="utf-8")

    must_have = [
        "library(phyloseq)",
        "ps <- Phylo_Objects$Species",
        "stopifnot(inherits(ps, \"phyloseq\"))",
        "otu <- as(otu_table(ps), \"matrix\")",
        "tax <- as(tax_table(ps), \"matrix\")",
        "meta <- as(sample_data(ps), \"data.frame\")",
        "if (!taxa_are_rows(otu_table(ps)))",
        "duplicate clade_name count",
        "aggregate(. ~ clade_name, data=ab, FUN=sum)",
        "species_abundance.tsv",
    ]
    for m in must_have:
        assert m in t


def test_validation_requires_species_abundance_and_phenotype_checks_and_fails():
    s2 = ROOT / "scripts/phase2e_06_validate_integrated_us_standardized_tables.py"
    t = s2.read_text(encoding="utf-8").lower()
    assert "species_abundance_exists" in t
    assert "original_sample_id" in t
    assert "case_status_is_pd_control" in t
    assert "case_status_pd_binary_agree" in t
    assert "donor_group_pd_consistency" in t
    assert "if status != \"pass\":" in t
    assert "raise systemexit(1)" in t


def test_no_disallowed_or_heavy_commands_present():
    files = [
        ROOT / "scripts/phase2e_05_extract_integrated_us_standardized_tables.py",
        ROOT / "scripts/_tmp_extract_integrated_us_species_phyloseq.R",
        ROOT / "scripts/phase2e_06_validate_integrated_us_standardized_tables.py",
    ]
    text = "\n".join(f.read_text(encoding="utf-8", errors="ignore").lower() for f in files)

    forbidden = [
        "wget",
        "curl",
        "prefetch",
        "fasterq-dump",
        "fastq",
        "sra",
        "raw metagenomic processing",
        "qiita",
        "curli carrier burden",
        "replication statistics",
        "phase 3",
    ]
    for term in forbidden:
        assert term not in text
