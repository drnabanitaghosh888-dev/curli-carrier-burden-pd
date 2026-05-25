# PHASE 2C-W2 Wallen Source Download and Inspection Runbook

## 1) Scope

This phase prepares manual acquisition and inspection of small Wallen processed/source files only.
No raw-read analysis is performed.

## 2) Safety Rules

- Download only small processed/source files from Zenodo.
- Do not download FASTQ.
- Do not download SRA.
- Do not run raw-read processing.
- Do not run HMMER/Kraken/HUMAnN/MetaPhlAn/Bowtie2.

## 3) Manual Download Instructions

From Wallen Zenodo record (`Zenodo_record_7246185`), manually download:
- `Source_Data_24Oct2022.xlsx`
- `Supplementary_Code_24Oct2022.zip`

## 4) Expected Local Folder

Use:

```bash
mkdir -p data/phase2c/manual_sources/wallen_prjna834801
```

After manual download:

```bash
cp -v ~/Downloads/Source_Data_24Oct2022.xlsx data/phase2c/manual_sources/wallen_prjna834801/
cp -v ~/Downloads/Supplementary_Code_24Oct2022.zip data/phase2c/manual_sources/wallen_prjna834801/
```

## 5) MD5 Verification Instructions

Run local checker:

```bash
python scripts/phase2c_02_check_wallen_source_files.py
```

It compares computed MD5 against:
- `4672da6a00ae951281441dd0a1620fdd` (workbook)
- `732c6273349d7a2efc0bc7d13d0c7d69` (zip)

## 6) Workbook Sheet-Inspection Instructions

Run:

```bash
python scripts/phase2c_03_inspect_wallen_source_workbook.py
```

This only lists sheets and basic dimensions/previews.  
It does not extract full data or perform analysis.

## 7) Supplementary-Code Zip Inspection Instructions

Run:

```bash
python scripts/phase2c_04_inspect_wallen_code_zip.py
```

This only inventories zip entries and likely file categories.  
It does not extract or execute code.

## 8) How to Interpret Sheet Inventory

Use candidate table-type tags to prepare mapping to:
- `sample_metadata.tsv`
- `species_abundance.tsv`
- `genus_abundance.tsv`
- `gene_family_abundance.tsv`
- `pathway_abundance.tsv`
- `run_accession_mapping.tsv`

Then run:

```bash
python scripts/phase2c_05_prepare_wallen_sheet_extraction_plan.py
```

## 9) Readiness Gates Before Next Phases

Wallen cannot be marked ready for Phase 2D until:
- species abundance verified
- sample metadata verified
- PD/control labels verified

Wallen cannot be marked ready for Phase 3 until:
- run-accession mapping verified

## 10) What Must Not Be Done in This Phase

- No FASTQ/SRA download
- No raw-read processing
- No metagenome profiling workflows
- No HMM screening on metagenomes
