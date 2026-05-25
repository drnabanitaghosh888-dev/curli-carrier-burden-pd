from __future__ import annotations

import argparse
import csv
import os
import shutil
import statistics
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "data/phase2e/reports"
REPORTS.mkdir(parents=True, exist_ok=True)
PROCESSED = ROOT / "data/phase2e/processed_tables/integrated_us_multicohort_pd"

RDATA = ROOT / "data/phase2e/manual_sources/integrated_us_multicohort_pd/github_jboktor/inspection_extract/PD_Metagenomic_Analysis-master/EDA_App/PhyloseqObj.RData"
README_SRC = ROOT / "data/phase2e/manual_sources/integrated_us_multicohort_pd/github_jboktor/inspection_extract/PD_Metagenomic_Analysis-master/README.md"
LICENSE_SRC = ROOT / "data/phase2e/manual_sources/integrated_us_multicohort_pd/github_jboktor/inspection_extract/PD_Metagenomic_Analysis-master/LICENSE"
METADATA_KEYS = ROOT / "data/phase2e/manual_sources/integrated_us_multicohort_pd/github_jboktor/inspection_extract/PD_Metagenomic_Analysis-master/files/metadata_keys.csv"

PLAN_TSV = REPORTS / "phase2e_integrated_us_standardized_extraction_plan.tsv"
STATUS_TXT = REPORTS / "phase2e_integrated_us_standardized_extraction_status.txt"
R_ENV_TSV = REPORTS / "phase2e_integrated_us_r_environment.tsv"
SUMMARY_TSV = REPORTS / "phase2e_integrated_us_standardized_extraction_summary.tsv"
R_LOG_TXT = REPORTS / "phase2e_integrated_us_r_extraction_log.txt"
PHENO_REPORT_TSV = REPORTS / "phase2e_integrated_us_phenotype_mapping_report.tsv"

TARGET_COMPONENT = "Phylo_Objects/Species"
MANDATORY_FILES = [
    "sample_metadata.tsv",
    "species_abundance.tsv",
    "README_or_data_dictionary.txt",
    "license_or_access_terms.txt",
]


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def _r_probe(rscript_path: str) -> dict[str, str]:
    info = {
        "candidate_rscript": rscript_path,
        "exists": "yes" if Path(rscript_path).exists() else "no",
        "phyloseq_available": "false",
        "phyloseq_version": "",
        "phyloseq_dir": "",
        "r_home": "",
        "lib_paths": "",
        "selected": "no",
        "notes": "",
    }

    probe = (
        'cat("Rscript=", normalizePath(Sys.which("Rscript")), "\\n", sep=""); '
        'cat("R_HOME=", R.home(), "\\n", sep=""); '
        'cat(".libPaths=", paste(.libPaths(), collapse=" | "), "\\n", sep=""); '
        'cat("phyloseq_available=", requireNamespace("phyloseq", quietly=TRUE), "\\n", sep=""); '
        'cat("phyloseq_dir=", system.file(package="phyloseq"), "\\n", sep=""); '
        'if(requireNamespace("phyloseq", quietly=TRUE)){cat("phyloseq_version=", as.character(packageVersion("phyloseq")), "\\n", sep="")}'
    )
    if not Path(rscript_path).exists():
        info["notes"] = "candidate_not_found"
        return info

    p = _run([rscript_path, "-e", probe])
    out = p.stdout.strip().splitlines() if p.returncode == 0 else []
    if p.returncode != 0:
        info["notes"] = (p.stderr or "probe_failed").strip()[:400]
        return info

    kv = {}
    for line in out:
        if "=" in line:
            k, v = line.split("=", 1)
            kv[k.strip()] = v.strip()

    info["phyloseq_available"] = kv.get("phyloseq_available", "false").lower()
    info["phyloseq_version"] = kv.get("phyloseq_version", "")
    info["phyloseq_dir"] = kv.get("phyloseq_dir", "")
    info["r_home"] = kv.get("R_HOME", "")
    info["lib_paths"] = kv.get(".libPaths", "")
    info["notes"] = "ok"
    return info


def _candidate_rscripts() -> list[str]:
    seen = set()
    cands: list[str] = []

    env = os.environ.get("PHASE2E_RSCRIPT", "").strip()
    if env:
        cands.append(env)

    cfg = ROOT / "config/phase2e_r_environment.tsv"
    if cfg.exists():
        for line in cfg.read_text(encoding="utf-8").splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            for part in parts:
                p = part.strip()
                if p.endswith("Rscript"):
                    cands.append(p)

    p = _run(["/bin/zsh", "-lc", "which -a Rscript 2>/dev/null || true"])
    for line in p.stdout.splitlines():
        if line.strip():
            cands.append(line.strip())

    p2 = _run(["/bin/zsh", "-lc", "command -v Rscript 2>/dev/null || true"])
    if p2.stdout.strip():
        cands.append(p2.stdout.strip())

    conda_prefix = os.environ.get("CONDA_PREFIX", "").strip()
    if conda_prefix:
        cands.append(str(Path(conda_prefix) / "bin" / "Rscript"))

    cands.extend([
        "/Users/krishnendu/miniconda3/envs/molevo/bin/Rscript",
        "/Users/krishnendu/miniforge3/envs/molevo/bin/Rscript",
        "/Users/krishnendu/anaconda3/envs/molevo/bin/Rscript",
        "/opt/homebrew/Caskroom/miniconda/base/envs/molevo/bin/Rscript",
        "/opt/homebrew/bin/Rscript",
        "/usr/local/bin/Rscript",
        "/Library/Frameworks/R.framework/Resources/bin/Rscript",
    ])

    out = []
    for c in cands:
        if c and c not in seen:
            seen.add(c)
            out.append(c)
    return out


def _write_r_env_report(rows: list[dict[str, str]]) -> None:
    cols = [
        "candidate_rscript",
        "exists",
        "phyloseq_available",
        "phyloseq_version",
        "phyloseq_dir",
        "r_home",
        "lib_paths",
        "selected",
        "notes",
    ]
    with R_ENV_TSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def _select_rscript() -> tuple[Optional[str], list[dict[str, str]]]:
    rows = []
    selected = None
    for cand in _candidate_rscripts():
        info = _r_probe(cand)
        if selected is None and info["phyloseq_available"] == "true":
            info["selected"] = "yes"
            selected = cand
        rows.append(info)
    return selected, rows


def _write_plan(selected_rscript: Optional[str], phyloseq_visible: bool) -> None:
    rows = [
        [
            "sample_metadata.tsv",
            "PhyloseqObj.RData",
            TARGET_COMPONENT,
            str(PROCESSED / "sample_metadata.tsv"),
            "yes",
            "high" if phyloseq_visible else "low",
            "yes",
            "execute mode only",
        ],
        [
            "species_abundance.tsv",
            "PhyloseqObj.RData",
            TARGET_COMPONENT,
            str(PROCESSED / "species_abundance.tsv"),
            "yes",
            "high" if phyloseq_visible else "low",
            "yes",
            "execute mode only",
        ],
        [
            "README_or_data_dictionary.txt",
            "README.md + metadata_keys.csv",
            "docs",
            str(PROCESSED / "README_or_data_dictionary.txt"),
            "yes",
            "high",
            "yes",
            "execute mode only",
        ],
        [
            "license_or_access_terms.txt",
            "LICENSE",
            "license",
            str(PROCESSED / "license_or_access_terms.txt"),
            "yes",
            "high",
            "yes",
            "execute mode only",
        ],
        [
            "run_accession_mapping.NOT_AVAILABLE_YET.txt",
            "metadata",
            "run mapping",
            str(PROCESSED / "run_accession_mapping.NOT_AVAILABLE_YET.txt"),
            "no",
            "low",
            "yes",
            "execute mode only",
        ],
    ]
    with PLAN_TSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "target_standard_file",
            "source_object_or_file",
            "source_component",
            "would_write_path",
            "required_for_replication",
            "extraction_feasibility",
            "execute_required",
            "notes",
        ])
        w.writerows(rows)


def _write_status(selected: Optional[str], env_rows: list[dict[str, str]], execute: bool, ok: bool) -> None:
    phy = "FALSE"
    ver = ""
    phy_dir = ""
    if selected:
        for r in env_rows:
            if r["candidate_rscript"] == selected:
                phy = r["phyloseq_available"].upper()
                ver = r["phyloseq_version"]
                phy_dir = r["phyloseq_dir"]
                break

    mode = "EXECUTE_MODE_REQUESTED" if execute else "DRY_RUN_ONLY"
    if not execute:
        status = "DRY_RUN_READY" if selected else "PHYLOSEQ_NOT_VISIBLE_TO_CODEX"
    else:
        status = "EXECUTE_SUCCESS" if ok else "EXECUTE_FAIL"

    STATUS_TXT.write_text(
        f"MODE = {mode}\n"
        f"STATUS = {status}\n"
        f"SELECTED_RSCRIPT = {selected or ''}\n"
        f"PHYLOSEQ_AVAILABLE = {phy}\n"
        f"PHYLOSEQ_VERSION = {ver}\n"
        f"PHYLOSEQ_DIR = {phy_dir}\n"
        f"RDATA_EXISTS = {RDATA.exists()}\n"
        f"TARGET_COMPONENT = {TARGET_COMPONENT}\n"
        f"FINAL_TABLES_WRITTEN = {'yes' if execute and ok else 'no'}\n",
        encoding="utf-8",
    )


def _write_tmp_r_extractor(path: Path) -> None:
    txt = r'''args <- commandArgs(trailingOnly=TRUE)
rdata_path <- args[1]
out_dir <- args[2]
log_path <- args[3]

log_lines <- character()
add_log <- function(k, v) {
  log_lines <<- c(log_lines, paste0(k, " = ", v))
}
flush_log <- function() {
  writeLines(log_lines, con=log_path)
}

suppressPackageStartupMessages(library(phyloseq))
add_log("selected Rscript", normalizePath(Sys.which("Rscript")))
add_log("phyloseq version", as.character(packageVersion("phyloseq")))

obj_names <- load(rdata_path)
if (!"Phylo_Objects" %in% obj_names) stop("Phylo_Objects not found")
if (is.null(Phylo_Objects$Species)) stop("Phylo_Objects$Species not found")
ps <- Phylo_Objects$Species
stopifnot(inherits(ps, "phyloseq"))

otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(otu_table(ps))) {
  otu <- t(otu)
}

tax <- as(tax_table(ps), "matrix")
meta <- as(sample_data(ps), "data.frame")
meta$sample_name <- rownames(meta)

add_log("nsamples", nsamples(ps))
add_log("ntaxa", ntaxa(ps))
add_log("taxa_are_rows", taxa_are_rows(otu_table(ps)))
add_log("otu_table dimensions", paste(dim(otu), collapse="x"))
add_log("tax_table dimensions", paste(dim(tax), collapse="x"))
add_log("sample_data dimensions", paste(dim(meta), collapse="x"))

rank_or_blank <- function(m, rank_name) {
  if (rank_name %in% colnames(m)) {
    out <- as.character(m[, rank_name])
    out[is.na(out)] <- ""
    out
  } else {
    rep("", nrow(m))
  }
}

rn <- rownames(otu)
tax <- tax[rn, , drop=FALSE]
sp <- rank_or_blank(tax, "Species")
if (!"Species" %in% colnames(tax) || any(sp == "")) {
  missing_idx <- which(sp == "")
  sp[missing_idx] <- rn[missing_idx]
  warning("Species rank missing/blank for one or more taxa; using taxa_names(ps)")
}

clade_name <- paste0(
  "k__", rank_or_blank(tax, "Kingdom"), "|",
  "p__", rank_or_blank(tax, "Phylum"), "|",
  "c__", rank_or_blank(tax, "Class"), "|",
  "o__", rank_or_blank(tax, "Order"), "|",
  "f__", rank_or_blank(tax, "Family"), "|",
  "g__", rank_or_blank(tax, "Genus"), "|",
  "s__", sp
)

ab <- data.frame(clade_name=clade_name, otu, check.names=FALSE, stringsAsFactors=FALSE)
dup_count <- sum(duplicated(ab$clade_name))
if (dup_count > 0) {
  ab <- aggregate(. ~ clade_name, data=ab, FUN=sum)
}

sm_path <- file.path(out_dir, "sample_metadata.tsv")
sp_path <- file.path(out_dir, "species_abundance.tsv")
write.table(meta, file=sm_path, sep="\t", quote=FALSE, row.names=FALSE)
write.table(ab, file=sp_path, sep="\t", quote=FALSE, row.names=FALSE)

add_log("species_abundance output path", sp_path)
add_log("duplicate clade_name count", dup_count)
add_log("final write status", ifelse(file.exists(sp_path), "SUCCESS", "FAIL"))
flush_log()
'''
    path.write_text(txt, encoding="utf-8")


def _write_summary(
    status: str,
    n_samples: str = "",
    n_species_rows: str = "",
    case_counts: str = "",
    pd_counts: str = "",
    abundance_cols: str = "",
    sample_rows: str = "",
    sample_match: str = "",
    col_min: str = "",
    col_median: str = "",
    col_max: str = "",
    dup_count: str = "",
    observed_donor_group_values: str = "",
    observed_pd_values: str = "",
    phenotype_mapping_status: str = "",
    donor_group_conflict_count: str = "",
    notes: str = "",
) -> None:
    with SUMMARY_TSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "status",
            "n_samples",
            "n_species_rows",
            "Case_status_counts",
            "PD_binary_counts",
            "abundance_sample_column_count",
            "sample_metadata_row_count",
            "sample_id_match_status",
            "abundance_column_sums_min",
            "abundance_column_sums_median",
            "abundance_column_sums_max",
            "duplicate_clade_name_count",
            "observed_donor_group_values",
            "observed_PD_values",
            "phenotype_mapping_status",
            "donor_group_conflict_count",
            "notes",
        ])
        w.writerow([
            status,
            n_samples,
            n_species_rows,
            case_counts,
            pd_counts,
            abundance_cols,
            sample_rows,
            sample_match,
            col_min,
            col_median,
            col_max,
            dup_count,
            observed_donor_group_values,
            observed_pd_values,
            phenotype_mapping_status,
            donor_group_conflict_count,
            notes,
        ])


def _counts_to_str(series: pd.Series) -> str:
    vc = series.fillna("NA").astype(str).value_counts(dropna=False).sort_index()
    return ";".join(f"{k}:{int(v)}" for k, v in vc.items())


def _parse_dup_count(log_path: Path) -> str:
    if not log_path.exists():
        return ""
    for line in log_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.lower().startswith("duplicate clade_name count") and "=" in line:
            return line.split("=", 1)[1].strip()
    return ""


def _apply_pd_mapping_and_write_report(metadata_path: Path) -> tuple[str, str, str, str, str]:
    md = pd.read_csv(metadata_path, sep="\t")
    required_for_mapping = ["PD", "donor_group"]
    missing_mapping = [c for c in required_for_mapping if c not in md.columns]
    if missing_mapping:
        raise RuntimeError(f"Missing required phenotype mapping columns: {', '.join(missing_mapping)}")

    if "original_sample_id" not in md.columns:
        if "id" in md.columns:
            md["original_sample_id"] = md["id"].astype(str)
        elif "sample_name" in md.columns:
            md["original_sample_id"] = md["sample_name"].astype(str)
        else:
            raise RuntimeError("Cannot create original_sample_id: neither id nor sample_name is available")
    if "study_group" in md.columns:
        # primary mapping must not be study_group
        pass

    pd_series = md["PD"].astype(str).str.strip()
    donor_group_series = md["donor_group"].astype(str).str.strip()
    case_status = pd_series.map({"Yes": "PD", "No": "Control"})
    pd_binary = pd_series.map({"Yes": 1, "No": 0})
    n_unmapped = int(case_status.isna().sum())
    if n_unmapped > 0:
        raise RuntimeError(f"Unmapped PD values found: {n_unmapped}")

    dg_conflict = ((donor_group_series == "PD") & (pd_series != "Yes")) | (
        donor_group_series.isin(["PC", "HC"]) & (pd_series != "No")
    )
    conflict_count = int(dg_conflict.sum())
    if conflict_count > 0:
        raise RuntimeError(f"donor_group/PD consistency conflicts: {conflict_count}")

    md["Case_status"] = case_status
    md["PD_binary"] = pd_binary.astype(int)
    required_output_cols = [
        "sample_name",
        "original_sample_id",
        "donor_id",
        "donor_group",
        "PD",
        "Case_status",
        "PD_binary",
        "host_age",
        "sex",
        "host_body_mass_index",
        "race",
    ]
    missing_output = [c for c in required_output_cols if c not in md.columns]
    if missing_output:
        raise RuntimeError(f"Missing required phenotype output columns: {', '.join(missing_output)}")
    ordered_cols = required_output_cols + [c for c in md.columns if c not in required_output_cols]
    md = md.loc[:, ordered_cols]
    md.to_csv(metadata_path, sep="\t", index=False)

    observed_donor_group_values = ";".join(sorted(donor_group_series.dropna().unique().tolist()))
    observed_pd_values = ";".join(sorted(pd_series.dropna().unique().tolist()))
    n_pd = int((pd_series == "Yes").sum())
    n_control = int((pd_series == "No").sum())
    status = "PASS"

    with PHENO_REPORT_TSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "source_column",
            "observed_values",
            "mapping_used",
            "pd_values",
            "control_values",
            "n_pd",
            "n_control",
            "n_unmapped",
            "donor_group_consistency_status",
            "donor_group_conflict_count",
            "status",
            "notes",
        ])
        w.writerow([
            "PD",
            observed_pd_values,
            "Yes->PD/1; No->Control/0",
            "Yes",
            "No",
            str(n_pd),
            str(n_control),
            str(n_unmapped),
            "PASS",
            str(conflict_count),
            status,
            "",
        ])
    return observed_donor_group_values, observed_pd_values, str(conflict_count), str(n_pd), str(n_control)


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase2E-D6 Integrated-US standardized extraction scaffold")
    parser.add_argument("--execute", action="store_true", help="Write standardized tables (manual use only)")
    args = parser.parse_args()

    selected, rows = _select_rscript()
    _write_r_env_report(rows)
    phyloseq_visible = selected is not None
    _write_plan(selected, phyloseq_visible)
    _write_tmp_r_extractor(ROOT / "scripts/_tmp_extract_integrated_us_species_phyloseq.R")

    if not args.execute:
        _write_summary(status="DRY_RUN_ONLY", notes="No tables written in dry-run mode")
        _write_status(selected, rows, execute=False, ok=False)
        return

    if selected is None:
        _write_summary(status="FAIL", notes="No Rscript with phyloseq available")
        _write_status(selected, rows, execute=True, ok=False)
        raise SystemExit(1)

    stage_dir = Path(tempfile.mkdtemp(prefix=".phase2e_stage_", dir=str(PROCESSED.parent)))
    ok = False
    try:
        r_helper = ROOT / "scripts/_tmp_extract_integrated_us_species_phyloseq.R"
        run = _run([selected, str(r_helper), str(RDATA), str(stage_dir), str(R_LOG_TXT)])
        if run.returncode != 0:
            _write_summary(status="FAIL", notes=(run.stderr or "R extraction failed")[:500])
            _write_status(selected, rows, execute=True, ok=False)
            raise SystemExit(1)

        if README_SRC.exists() and METADATA_KEYS.exists():
            text = "# README\n" + README_SRC.read_text(encoding="utf-8", errors="ignore")
            text += "\n\n# metadata_keys.csv\n" + METADATA_KEYS.read_text(encoding="utf-8", errors="ignore")
            (stage_dir / "README_or_data_dictionary.txt").write_text(text, encoding="utf-8")
        if LICENSE_SRC.exists():
            shutil.copy2(LICENSE_SRC, stage_dir / "license_or_access_terms.txt")

        missing = [name for name in MANDATORY_FILES if not (stage_dir / name).exists()]
        if missing:
            _write_summary(status="FAIL", notes=f"Missing mandatory files in stage: {', '.join(missing)}")
            _write_status(selected, rows, execute=True, ok=False)
            raise SystemExit(1)

        (stage_dir / "run_accession_mapping.NOT_AVAILABLE_YET.txt").write_text(
            "Run accession mapping is not explicitly available in current extraction scaffold.\n",
            encoding="utf-8",
        )

        PROCESSED.mkdir(parents=True, exist_ok=True)
        for fname in MANDATORY_FILES + ["run_accession_mapping.NOT_AVAILABLE_YET.txt"]:
            src = stage_dir / fname
            dst = PROCESSED / fname
            os.replace(src, dst)

        smdf = pd.read_csv(PROCESSED / "sample_metadata.tsv", sep="\t")
        spdf = pd.read_csv(PROCESSED / "species_abundance.tsv", sep="\t")
        observed_donor_group_values, observed_pd_values, conflict_count, _n_pd, _n_control = _apply_pd_mapping_and_write_report(
            PROCESSED / "sample_metadata.tsv"
        )
        smdf = pd.read_csv(PROCESSED / "sample_metadata.tsv", sep="\t")

        abundance_samples = [c for c in spdf.columns if c != "clade_name"]
        sample_rows = int(smdf.shape[0])
        n_species_rows = int(spdf.shape[0])
        n_samples = len(abundance_samples)

        md_samples = set(smdf["sample_name"].astype(str).tolist()) if "sample_name" in smdf.columns else set()
        sample_match = "PASS" if md_samples == set(abundance_samples) else "FAIL"

        col_sums = pd.to_numeric(spdf[abundance_samples].sum(axis=0), errors="coerce") if abundance_samples else pd.Series(dtype=float)
        if len(col_sums) > 0:
            col_min = str(float(col_sums.min()))
            col_median = str(float(statistics.median(col_sums.tolist())))
            col_max = str(float(col_sums.max()))
        else:
            col_min = ""
            col_median = ""
            col_max = ""

        case_counts = _counts_to_str(smdf["Case_status"]) if "Case_status" in smdf.columns else "MISSING_COLUMN"
        pd_counts = _counts_to_str(smdf["PD_binary"]) if "PD_binary" in smdf.columns else "MISSING_COLUMN"
        dup_count = _parse_dup_count(R_LOG_TXT)

        _write_summary(
            status="SUCCESS",
            n_samples=str(n_samples),
            n_species_rows=str(n_species_rows),
            case_counts=case_counts,
            pd_counts=pd_counts,
            abundance_cols=str(n_samples),
            sample_rows=str(sample_rows),
            sample_match=sample_match,
            col_min=col_min,
            col_median=col_median,
            col_max=col_max,
            dup_count=dup_count,
            observed_donor_group_values=observed_donor_group_values,
            observed_pd_values=observed_pd_values,
            phenotype_mapping_status="PASS",
            donor_group_conflict_count=conflict_count,
            notes="",
        )

        _write_status(selected, rows, execute=True, ok=True)
        ok = True
    finally:
        shutil.rmtree(stage_dir, ignore_errors=True)
        if not ok and args.execute and not SUMMARY_TSV.exists():
            _write_summary(status="FAIL", notes="Extraction failed before summary generation")


if __name__ == "__main__":
    main()
