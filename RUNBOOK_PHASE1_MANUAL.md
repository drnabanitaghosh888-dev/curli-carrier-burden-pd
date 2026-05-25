# RUNBOOK: PHASE 1 (MANUAL, STEP-GATED)

Project: `curli_vag_pd`  
Goal: Build a reproducible curli reference database workflow for `csgA`-`csgG`.

Important constraints:
- Do not run large downloads.
- Do not run metagenomic/SRA workflows.
- Run each checkpoint in order and review outputs before proceeding.

---

## Checkpoint A: Environment Check

Run:

```bash
python scripts/phase1_00_check_environment.py
```

Review:
- `data/phase1/reports/environment_check.md`

Proceed criteria:
- Required tools found (`mafft`, `hmmbuild`, `hmmpress`, `hmmsearch`, `seqkit`, `cd-hit`, `diamond`, `blastp`)
- Python imports pass
- UniProt connectivity check passes

---

## Checkpoint B: Fetch Dry Run

Run:

```bash
python scripts/phase1_01_fetch_curli_sequences.py --dry-run
```

Review:
- Planned UniProt queries shown in console
- Planned output paths under `data/phase1/raw/`

Proceed criteria:
- Queries are specific to `csgA`-`csgG`
- Taxonomic scope and max sequence settings are acceptable

---

## Checkpoint C: Fetch Real Sequences

Run:

```bash
python scripts/phase1_01_fetch_curli_sequences.py
```

Review:
- `data/phase1/raw/uniprot_download_manifest.tsv`
- If present: `data/phase1/raw/uniprot_failed_queries.tsv`

Proceed criteria:
- All genes have output FASTA/TSV files
- No unexplained fetch failures

---

## Checkpoint D: Curate Sequences

Run:

```bash
python scripts/phase1_02_filter_and_curate.py
```

Review:
- `data/phase1/curated/curli_filtering_summary.tsv`

Proceed criteria:
- Curated outputs exist for all genes
- Low-count genes are identified and noted

---

## Checkpoint E: Build Alignments

Run:

```bash
python scripts/phase1_03_build_alignments.py
```

Review:
- `data/phase1/alignments/alignment_summary.tsv`

Proceed criteria:
- Alignment files generated for genes with valid curated input
- No unexplained MAFFT failures

---

## Checkpoint F: Build HMMs

Run:

```bash
python scripts/phase1_04_build_hmms.py
```

Review:
- `data/phase1/hmms/hmm_build_summary.tsv`

Proceed criteria:
- Per-gene HMMs built
- Combined HMM exists
- `hmmpress` index files exist (`.h3f`, `.h3i`, `.h3m`, `.h3p`)

---

## Checkpoint G: Validate HMMs

Run:

```bash
python scripts/phase1_05_validate_hmms.py
```

Review:
- `data/phase1/validation/hmm_validation_summary.tsv`

Proceed criteria:
- Genes classified (`reliable`, `borderline`, `weak`)
- Weak/borderline genes documented for interpretation

---

## Checkpoint H: Final Report

Run:

```bash
python scripts/phase1_06_make_report.py
```

Review:
- `data/phase1/reports/phase1_status.txt`
- `data/phase1/reports/phase1_curli_reference_report.md`

Proceed criteria:
- Final status is `PASS` or `PASS_WITH_WARNINGS`
- If `FAIL`, resolve missing artifacts and rerun report

---

## Manual Execution Order (Quick Reference)

1. `python scripts/phase1_00_check_environment.py`
2. `python scripts/phase1_01_fetch_curli_sequences.py --dry-run`
3. `python scripts/phase1_01_fetch_curli_sequences.py`
4. `python scripts/phase1_02_filter_and_curate.py`
5. `python scripts/phase1_03_build_alignments.py`
6. `python scripts/phase1_04_build_hmms.py`
7. `python scripts/phase1_05_validate_hmms.py`
8. `python scripts/phase1_06_make_report.py`
