from pathlib import Path
import py_compile

ROOT = Path(__file__).resolve().parents[1]

REQUIRED_DIRS = [
    "config",
    "scripts",
    "data/phase1/raw",
    "data/phase1/curated",
    "data/phase1/alignments",
    "data/phase1/hmms",
    "data/phase1/validation",
    "data/phase1/reports",
    "logs",
    "tests",
    "results/phase1_reference_package",
]

REQUIRED_FILES = [
    "README.md",
    "environment.yml",
    "Makefile",
    "RUNBOOK_PHASE1_MANUAL.md",
    "config/phase1_config.yml",
    "config/curli_gene_roles.yml",
    "config/curli_query_terms.yml",
    "scripts/phase1_00_check_environment.py",
    "scripts/phase1_01_fetch_curli_sequences.py",
    "scripts/phase1_02_filter_and_curate.py",
    "scripts/phase1_03_build_alignments.py",
    "scripts/phase1_04_build_hmms.py",
    "scripts/phase1_05_validate_hmms.py",
    "scripts/phase1_06_make_report.py",
    "tests/test_phase1_structure.py",
]

SCRIPT_FILES = [
    "scripts/phase1_00_check_environment.py",
    "scripts/phase1_01_fetch_curli_sequences.py",
    "scripts/phase1_02_filter_and_curate.py",
    "scripts/phase1_03_build_alignments.py",
    "scripts/phase1_04_build_hmms.py",
    "scripts/phase1_05_validate_hmms.py",
    "scripts/phase1_06_make_report.py",
]

EXPECTED_GENES = ["csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG"]


def test_required_directories_exist():
    missing = [d for d in REQUIRED_DIRS if not (ROOT / d).is_dir()]
    assert not missing, f"Missing required directories: {missing}"


def test_required_files_exist_and_nonempty():
    missing = []
    empty = []
    for f in REQUIRED_FILES:
        path = ROOT / f
        if not path.is_file():
            missing.append(f)
        elif path.stat().st_size == 0:
            empty.append(f)
    assert not missing, f"Missing required files: {missing}"
    assert not empty, f"Empty required files: {empty}"


def test_phase1_scripts_compile():
    failures = []
    for script in SCRIPT_FILES:
        path = ROOT / script
        try:
            py_compile.compile(str(path), doraise=True)
        except Exception as exc:
            failures.append(f"{script}: {exc}")
    assert not failures, "Script syntax failures: " + "; ".join(failures)


def test_curli_gene_roles_mentions_all_genes():
    text = (ROOT / "config/curli_gene_roles.yml").read_text()
    missing = [gene for gene in EXPECTED_GENES if f"{gene}:" not in text]
    assert not missing, f"curli_gene_roles.yml missing genes: {missing}"


def test_query_terms_mentions_all_genes():
    text = (ROOT / "config/curli_query_terms.yml").read_text()
    missing = [gene for gene in EXPECTED_GENES if f"{gene}:" not in text]
    assert not missing, f"curli_query_terms.yml missing genes: {missing}"


def test_makefile_has_expected_phase1_targets():
    text = (ROOT / "Makefile").read_text()
    expected_targets = [
        "env-check",
        "phase1-fetch-dry-run",
        "phase1-fetch",
        "phase1-curate",
        "phase1-align",
        "phase1-hmms",
        "phase1-validate",
        "phase1-report",
        "phase1-all",
        "phase1-clean-temp",
    ]
    missing = [target for target in expected_targets if f"{target}:" not in text]
    assert not missing, f"Makefile missing targets: {missing}"


def test_runbook_contains_manual_checkpoints():
    text = (ROOT / "RUNBOOK_PHASE1_MANUAL.md").read_text()
    expected_checkpoints = [
        "Checkpoint A",
        "Checkpoint B",
        "Checkpoint C",
        "Checkpoint D",
        "Checkpoint E",
        "Checkpoint F",
        "Checkpoint G",
        "Checkpoint H",
    ]
    missing = [checkpoint for checkpoint in expected_checkpoints if checkpoint not in text]
    assert not missing, f"Runbook missing checkpoints: {missing}"
