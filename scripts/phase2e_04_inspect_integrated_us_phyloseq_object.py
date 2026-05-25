from __future__ import annotations

import csv
import json
import re
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "data/phase2e/reports"
REPORTS.mkdir(parents=True, exist_ok=True)

PHYLOSEQ_RDATA = ROOT / "data/phase2e/manual_sources/integrated_us_multicohort_pd/github_jboktor/inspection_extract/PD_Metagenomic_Analysis-master/EDA_App/PhyloseqObj.RData"
META_FILES = [
    ROOT / "data/phase2e/manual_sources/integrated_us_multicohort_pd/github_jboktor/inspection_extract/PD_Metagenomic_Analysis-master/files/metadata_phyloseq.csv",
    ROOT / "data/phase2e/manual_sources/integrated_us_multicohort_pd/github_jboktor/inspection_extract/PD_Metagenomic_Analysis-master/files/metadata_phyloseq_TBC.csv",
    ROOT / "data/phase2e/manual_sources/integrated_us_multicohort_pd/github_jboktor/inspection_extract/PD_Metagenomic_Analysis-master/files/metadata_phyloseq_RUSH.csv",
    ROOT / "data/phase2e/manual_sources/integrated_us_multicohort_pd/github_jboktor/inspection_extract/PD_Metagenomic_Analysis-master/files/metadata_phyloseq_Bonn.csv",
    ROOT / "data/phase2e/manual_sources/integrated_us_multicohort_pd/github_jboktor/inspection_extract/PD_Metagenomic_Analysis-master/files/metadata_phyloseq_SHANGHAI.csv",
    ROOT / "data/phase2e/manual_sources/integrated_us_multicohort_pd/github_jboktor/inspection_extract/PD_Metagenomic_Analysis-master/files/metadata_keys.csv",
]
TABLE_S10 = ROOT / "data/phase2e/manual_sources/integrated_us_multicohort_pd/zenodo/supplementary_tables_inspection/Supplementary-Tables/Table_S10 Sample-Metadata.xlsx"

PHYLOSEQ_OBJECT_INVENTORY = REPORTS / "phase2e_integrated_us_phyloseq_object_inventory.tsv"
PHYLOSEQ_SAMPLE_METADATA_INVENTORY = REPORTS / "phase2e_integrated_us_phyloseq_sample_metadata_inventory.tsv"
PHYLOSEQ_TAXA_INVENTORY = REPORTS / "phase2e_integrated_us_phyloseq_taxa_inventory.tsv"
PHYLOSEQ_COMPONENT_INVENTORY = REPORTS / "phase2e_integrated_us_phyloseq_component_inventory.tsv"
CANDIDATE_COMPONENT_SELECTION = REPORTS / "phase2e_integrated_us_candidate_component_selection.tsv"
LABEL_MAPPING_CANDIDATES = REPORTS / "phase2e_integrated_us_label_mapping_candidates.tsv"
TAXONOMY_RANK_READINESS = REPORTS / "phase2e_integrated_us_taxonomy_rank_readiness.tsv"
METADATA_COLUMN_INVENTORY = REPORTS / "phase2e_integrated_us_metadata_column_inventory.tsv"
EXTRACTION_PLAN = REPORTS / "phase2e_integrated_us_extraction_plan.tsv"
D5_STATUS = REPORTS / "phase2e_integrated_us_d5_status.txt"
D5B_STATUS = REPORTS / "phase2e_integrated_us_d5b_status.txt"
D5B_REPORT = REPORTS / "phase2e_integrated_us_d5b_report.md"

TMP_R_SCRIPT = ROOT / "scripts/_tmp_inspect_integrated_us_phyloseq.R"
TMP_JSON = REPORTS / "_tmp_phase2e_integrated_us_phyloseq_inspection.json"

R_SCRIPT_TEXT = r'''
args <- commandArgs(trailingOnly=TRUE)
rdata_path <- args[1]
out_json <- args[2]

safe_len <- function(x){
  if (is.null(x)) return(NA_integer_)
  as.integer(length(x))
}

safe_head_str <- function(x, n=10){
  if (is.null(x)) return(list())
  as.list(as.character(utils::head(x, n)))
}

safe_colnames <- function(x, n=100){
  cn <- colnames(x)
  if (is.null(cn)) return(list())
  as.list(as.character(utils::head(cn, n)))
}

inspect_one <- function(obj, nm, path, depth, phyloseq_available){
  cls_vec <- class(obj)
  cls <- paste(cls_vec, collapse=';')
  is_phy <- any(grepl('phyloseq', cls_vec, ignore.case=TRUE))
  has_sample_data <- FALSE
  has_otu <- FALSE
  has_tax <- FALSE
  sample_count <- NA
  taxa_count <- NA
  taxa_are_rows <- NA
  sample_cols <- list()
  sample_ids <- list()
  tax_ranks <- list()
  taxa_names <- list()
  notes <- c()

  if (is_phy && phyloseq_available) {
    tryCatch({
      sample_count <- phyloseq::nsamples(obj)
      taxa_count <- phyloseq::ntaxa(obj)
      sm <- phyloseq::sample_data(obj)
      ot <- phyloseq::otu_table(obj)
      tt <- phyloseq::tax_table(obj)
      has_sample_data <- !is.null(sm)
      has_otu <- !is.null(ot)
      has_tax <- !is.null(tt)
      sample_cols <- safe_colnames(sm, 200)
      sample_ids <- safe_head_str(rownames(sm), 10)
      tax_ranks <- safe_colnames(tt, 50)
      taxa_names <- safe_head_str(phyloseq::taxa_names(obj), 10)
      taxa_are_rows <- phyloseq::taxa_are_rows(ot)
    }, error=function(e){
      notes <<- c(notes, paste('phyloseq inspection error:', e$message))
    })
  } else if (is_phy && !phyloseq_available) {
    notes <- c(notes, 'phyloseq package unavailable; deep slot inspection not executed')
  }

  row <- list(
    component_path = path,
    component_name = nm,
    component_class = cls,
    is_phyloseq = is_phy,
    sample_count = sample_count,
    taxa_count = taxa_count,
    taxa_are_rows = taxa_are_rows,
    sample_data_n_columns = safe_len(sample_cols),
    sample_data_columns_preview = sample_cols,
    tax_rank_names = tax_ranks,
    taxa_names_preview = taxa_names,
    has_sample_data = has_sample_data,
    has_otu_table = has_otu,
    has_tax_table = has_tax,
    inspection_status = 'ok',
    notes = notes,
    depth = depth
  )
  list(row=row, children=list())
}

walk_obj <- function(obj, nm, path, depth, max_depth, phyloseq_available){
  out <- list()
  cur <- inspect_one(obj, nm, path, depth, phyloseq_available)
  out[[length(out)+1]] <- cur$row

  if (depth >= max_depth) return(out)

  if (is.list(obj)) {
    nms <- names(obj)
    for (i in seq_along(obj)) {
      child <- obj[[i]]
      child_nm <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else paste0('idx_', i)
      child_path <- paste0(path, '/', child_nm)
      child_rows <- walk_obj(child, child_nm, child_path, depth+1, max_depth, phyloseq_available)
      for (r in child_rows) out[[length(out)+1]] <- r
    }
  }
  out
}

result <- list(status='OK', rscript_available=TRUE, phyloseq_available=FALSE, errors=list(), objects=list(), components=list())
phyloseq_available <- requireNamespace('phyloseq', quietly=TRUE)
result$phyloseq_available <- phyloseq_available

obj_names <- c()
tryCatch({
  obj_names <- load(rdata_path)
}, error=function(e){
  result$status <<- 'LOAD_FAILED'
  result$errors <<- c(result$errors, as.character(e$message))
})

for (nm in obj_names) {
  obj <- get(nm)
  top <- inspect_one(obj, nm, nm, 0, phyloseq_available)
  result$objects[[length(result$objects)+1]] <- top$row
  rows <- walk_obj(obj, nm, nm, 0, 3, phyloseq_available)
  for (r in rows) result$components[[length(result$components)+1]] <- r
}

if (length(result$components) == 0) {
  result$status <- 'EMPTY_COMPONENT_INVENTORY'
}

if (!phyloseq_available) {
  result$status <- ifelse(result$status == 'OK', 'NEEDS_PHYLOSEQ_PACKAGE', result$status)
}

jsonlite::write_json(result, out_json, auto_unbox=TRUE, pretty=TRUE)
'''


def write_tsv(path: Path, columns: list[str], rows: list[list[str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(columns)
        w.writerows(rows)


def run_r_inspection() -> dict:
    TMP_R_SCRIPT.write_text(R_SCRIPT_TEXT, encoding="utf-8")
    rscript = subprocess.run(
        ["Rscript", str(TMP_R_SCRIPT), str(PHYLOSEQ_RDATA), str(TMP_JSON)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if rscript.returncode != 0:
        return {
            "status": "LOAD_FAILED",
            "rscript_available": True,
            "phyloseq_available": False,
            "errors": [rscript.stderr.strip() or "Rscript inspection failed"],
            "objects": [],
            "components": [],
        }
    if not TMP_JSON.exists():
        return {
            "status": "LOAD_FAILED",
            "rscript_available": True,
            "phyloseq_available": False,
            "errors": ["R script did not emit JSON"],
            "objects": [],
            "components": [],
        }
    return json.loads(TMP_JSON.read_text(encoding="utf-8"))


def likely_cols(columns: list[str]) -> tuple[list[str], list[str], list[str], list[str]]:
    c = [x.strip() for x in columns]
    sid = [x for x in c if re.search(r"sample|subject|id$|sample_id|subject_id|biosample|run", x, re.I)]
    disease = [x for x in c if re.search(r"pd|parkinson|disease|control|case|status|diagnosis|group|condition|phenotype", x, re.I)]
    cohort = [x for x in c if re.search(r"cohort|study|site|center|country|dataset|source|rush|tbc|bonn|shanghai", x, re.I)]
    cov = [x for x in c if re.search(r"age|sex|bmi|med|drug|updrs|duration|constipation|hny|score", x, re.I)]
    return sid[:10], disease[:10], cohort[:10], cov[:10]


def inspect_metadata_files() -> list[list[str]]:
    rows: list[list[str]] = []
    for p in META_FILES + [TABLE_S10]:
        if not p.exists():
            rows.append([p.name, str(p), "missing", "", "", "", "", "", "", "", ""])
            continue
        try:
            if p.suffix.lower() == ".csv":
                df = pd.read_csv(p)
            else:
                df = pd.read_excel(p)
            cols = [str(x) for x in df.columns]
            sid, disease, cohort, cov = likely_cols(cols)
            miss = ""
            if disease:
                miss_parts = []
                for col in disease[:5]:
                    frac = float(df[col].isna().mean()) if col in df.columns else 1.0
                    miss_parts.append(f"{col}:{frac:.3f}")
                miss = "; ".join(miss_parts)
            rows.append([
                p.name,
                str(p.relative_to(ROOT)),
                "ok",
                str(df.shape[0]),
                str(df.shape[1]),
                "; ".join(cols[:30]),
                "; ".join(sid),
                "; ".join(disease),
                "; ".join(cohort),
                "; ".join(cov),
                miss,
            ])
        except Exception as e:
            rows.append([
                p.name,
                str(p.relative_to(ROOT)),
                f"unreadable:{type(e).__name__}",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ])
    return rows


def lower_join(values: list[str]) -> str:
    return " ".join([str(v).lower() for v in values if v is not None])


def component_rows_to_tables(components: list[dict], meta_rows: list[list[str]]) -> tuple[list[list[str]], list[list[str]], list[list[str]], list[list[str]]]:
    comp_inv_rows: list[list[str]] = []
    cand_rows: list[list[str]] = []
    label_rows: list[list[str]] = []
    tax_rows: list[list[str]] = []

    for c in components:
        sample_cols = [str(x) for x in c.get("sample_data_columns_preview", []) if x is not None]
        tax_ranks = [str(x) for x in c.get("tax_rank_names", []) if x is not None]
        taxa_preview = [str(x) for x in c.get("taxa_names_preview", []) if x is not None]
        notes = "; ".join([str(x) for x in c.get("notes", [])]) if isinstance(c.get("notes"), list) else str(c.get("notes", ""))

        comp_inv_rows.append([
            str(c.get("component_path", "")),
            str(c.get("component_name", "")),
            str(c.get("component_class", "")),
            str(bool(c.get("is_phyloseq", False))).lower(),
            str(c.get("sample_count", "")),
            str(c.get("taxa_count", "")),
            str(c.get("taxa_are_rows", "")),
            str(c.get("sample_data_n_columns", "")),
            "; ".join(sample_cols[:30]),
            "; ".join(tax_ranks[:20]),
            "; ".join(taxa_preview[:10]),
            str(c.get("inspection_status", "")),
            notes,
        ])

        is_phy = bool(c.get("is_phyloseq", False))
        has_sample = bool(c.get("has_sample_data", False))
        has_otu = bool(c.get("has_otu_table", False))
        has_tax = bool(c.get("has_tax_table", False))
        has_species_rank = any(re.search(r"species|\bs\b", r, re.I) for r in tax_ranks)
        has_genus_rank = any(re.search(r"genus|\bg\b", r, re.I) for r in tax_ranks)
        metaphlan_like_species = any("s__" in t for t in taxa_preview)

        joined = lower_join([c.get("component_path", ""), c.get("component_name", ""), ";".join(sample_cols)])
        us_hint = ("tbc" in joined) or ("rush" in joined)
        non_us_hint = ("bonn" in joined) or ("shanghai" in joined)

        candidate = "no"
        reason = ""
        block = ""
        if is_phy and has_sample and has_otu and has_tax:
            if us_hint and not non_us_hint:
                candidate = "yes"
                reason = "phyloseq component has sample_data+otu_table+tax_table and US-cohort hint (TBC/RUSH)"
                block = "none"
            elif non_us_hint and not us_hint:
                candidate = "no"
                reason = "component appears cohort-specific non-US (Bonn/Shanghai)"
                block = "component_not_selected"
            else:
                candidate = "needs_manual_selection"
                reason = "phyloseq-capable component but integrated-US scope is not uniquely identifiable"
                block = "component_not_selected"
        else:
            reason = "component lacks required phyloseq structures"
            block = "component_not_selected"

        cand_rows.append([
            str(c.get("component_path", "")),
            str(c.get("component_name", "")),
            str(is_phy).lower(),
            str(c.get("sample_count", "")),
            str(c.get("taxa_count", "")),
            str(has_sample).lower(),
            str(has_otu).lower(),
            str(has_tax).lower(),
            str(has_species_rank).lower(),
            str(has_genus_rank).lower(),
            candidate,
            reason,
            block,
            "",
        ])

        # label mapping from component sample_data columns
        if candidate in {"yes", "needs_manual_selection"}:
            disease_cols = [x for x in sample_cols if re.search(r"pd|disease|diagnosis|status|group|condition|phenotype|control", x, re.I)]
            sid_cols = [x for x in sample_cols if re.search(r"sample|subject|id|biosample|run", x, re.I)]
            cohort_cols = [x for x in sample_cols if re.search(r"cohort|study|site|center|country|rush|tbc|bonn|shanghai|source", x, re.I)]
            label_rows.append([
                "phyloseq_component",
                str(c.get("component_path", "")),
                str(c.get("sample_count", "")),
                sid_cols[0] if sid_cols else "",
                disease_cols[0] if disease_cols else "",
                "unknown_from_structure_only",
                "PD|Parkinson|case (candidate)",
                "Control|HC (candidate)",
                cohort_cols[0] if cohort_cols else "",
                "moderate" if (sid_cols and disease_cols) else "low",
                "observed label values need explicit extraction step",
            ])

        species_rank_column = ""
        genus_rank_column = ""
        for r in tax_ranks:
            if not species_rank_column and re.search(r"species|\bs\b", r, re.I):
                species_rank_column = r
            if not genus_rank_column and re.search(r"genus|\bg\b", r, re.I):
                genus_rank_column = r

        ready_species = "yes" if has_species_rank else ("no" if not metaphlan_like_species else "no")
        tax_block = "none" if has_species_rank else ("species_rank_unclear")
        tax_note = ""
        if not has_species_rank and metaphlan_like_species:
            tax_note = "possible_metaphlan_species_names_needing_parser"

        tax_rows.append([
            str(c.get("component_path", "")),
            "; ".join(tax_ranks[:20]),
            str(has_species_rank).lower(),
            species_rank_column,
            str(has_genus_rank).lower(),
            genus_rank_column,
            "; ".join(taxa_preview[:8]),
            ready_species,
            tax_block,
            tax_note,
        ])

    # Metadata file label candidates
    for r in meta_rows:
        if r[2] != "ok":
            continue
        label_rows.append([
            f"metadata_file:{r[0]}",
            "",
            r[3],
            (r[6].split("; ")[0] if r[6] else ""),
            (r[7].split("; ")[0] if r[7] else ""),
            "unknown_from_header_only",
            "PD|Parkinson|case (candidate)",
            "Control|HC (candidate)",
            (r[8].split("; ")[0] if r[8] else ""),
            "moderate" if (r[6] and r[7]) else "low",
            "column-level candidate from metadata header inspection",
        ])

    return comp_inv_rows, cand_rows, label_rows, tax_rows


def build_extraction_plan(cand_rows: list[list[str]], tax_rows: list[list[str]], meta_rows: list[list[str]]) -> tuple[list[list[str]], str]:
    yes_components = [r for r in cand_rows if r[10] == "yes"]
    needs_select = [r for r in cand_rows if r[10] == "needs_manual_selection"]

    species_ready = any(r[7] == "yes" for r in tax_rows if any(y[0] == r[0] for y in yes_components + needs_select))
    species_parser_possible = any("possible_metaphlan_species_names_needing_parser" in r[9] for r in tax_rows)

    has_sample_cols = any((r[6] and r[7]) for r in meta_rows if r[2] == "ok")

    comp_block = "none"
    status = "PARTIAL_READY_NEEDS_LIST_COMPONENT_INSPECTION"
    if not yes_components and not needs_select:
        status = "NOT_READY_NO_COMPONENT_PHYLOSEQ_OBJECT"
        comp_block = "component_not_selected"
    elif needs_select and not yes_components:
        status = "PARTIAL_READY_MULTIPLE_COMPONENTS_NEED_SELECTION"
        comp_block = "component_not_selected"
    elif yes_components and not has_sample_cols:
        status = "PARTIAL_READY_NEEDS_MANUAL_COLUMN_MAPPING"
        comp_block = "sample_id_unclear"
    elif yes_components and not species_ready:
        status = "PARTIAL_READY_SPECIES_RANK_UNCLEAR"
        comp_block = "species_rank_unclear"
    elif yes_components and has_sample_cols and species_ready:
        status = "READY_FOR_STANDARDIZED_EXTRACTION"
        comp_block = "none"

    plan = [
        [
            "sample_metadata.tsv",
            "candidate phyloseq component + metadata_phyloseq*.csv + Table_S10",
            "sample_data and metadata columns",
            "yes",
            "high" if yes_components and has_sample_cols else ("medium" if needs_select else "low"),
            "select Integrated-US component, export sample_data, reconcile sample ID and disease labels with metadata files",
            "none" if yes_components and has_sample_cols else ("component_not_selected" if not yes_components else "sample_id_unclear"),
            "do not extract in D5b; plan only",
        ],
        [
            "species_abundance.tsv",
            "candidate phyloseq component",
            "otu_table joined to species taxonomy",
            "yes",
            "high" if yes_components and species_ready else ("medium" if yes_components or needs_select else "low"),
            "export otu_table and map species rank or parse metaphlan species labels",
            "none" if yes_components and species_ready else ("species_rank_unclear" if (yes_components or needs_select) else "component_not_selected"),
            "do not mark high until species rank/component are unambiguous",
        ],
        [
            "genus_abundance.tsv",
            "candidate phyloseq component",
            "otu_table joined to genus taxonomy",
            "no",
            "medium" if (yes_components or needs_select) else "low",
            "export genus rank if present",
            "component_not_selected" if not (yes_components or needs_select) else "none",
            "optional",
        ],
        [
            "README_or_data_dictionary.txt",
            "README.md + metadata_keys.csv",
            "documentation",
            "yes",
            "high",
            "copy/normalize docs",
            "none",
            "already available",
        ],
        [
            "license_or_access_terms.txt",
            "LICENSE",
            "license text",
            "yes",
            "high",
            "copy license text",
            "none",
            "already available",
        ],
        [
            "run_accession_mapping.tsv",
            "metadata files",
            "run/accession columns if explicit",
            "no",
            "low",
            "extract run/BioSample columns if present",
            "run_accession_mapping_missing",
            "pending unless explicit mapping needed",
        ],
    ]

    # hard guard: never allow READY if top-level list-only unresolved evidence
    top_level_only_list = False
    if cand_rows:
        # no component with is_phyloseq true means unresolved list-only for practical extraction
        any_phy = any(r[2] == "true" for r in cand_rows)
        top_level_only_list = not any_phy
    if top_level_only_list:
        if status == "READY_FOR_STANDARDIZED_EXTRACTION":
            status = "PARTIAL_READY_NEEDS_LIST_COMPONENT_INSPECTION"

    return plan, status


def main() -> None:
    try:
        r_result = run_r_inspection()
        objects = r_result.get("objects", [])
        components = r_result.get("components", [])

        obj_rows = []
        sm_rows = []
        tx_rows = []
        for obj in objects:
            obj_rows.append([
                obj.get("component_name", obj.get("object_name", "")),
                obj.get("component_class", obj.get("object_class", "")),
                str(bool(obj.get("is_phyloseq", obj.get("is_phyloseq_like", False)))).lower(),
                str(obj.get("sample_count", "")),
                str(obj.get("taxa_count", "")),
                str(obj.get("taxa_are_rows", "")),
                r_result.get("status", ""),
                "; ".join(r_result.get("errors", [])),
            ])
            for col in obj.get("sample_data_columns_preview", obj.get("sample_columns", []))[:300]:
                sm_rows.append([
                    obj.get("component_name", obj.get("object_name", "")),
                    str(col),
                    "; ".join([str(x) for x in obj.get("sample_ids_preview", [])[:10]]),
                ])
            for rank in obj.get("tax_rank_names", [])[:100]:
                tx_rows.append([
                    obj.get("component_name", obj.get("object_name", "")),
                    str(rank),
                    "; ".join([str(x) for x in obj.get("taxa_names_preview", [])[:10]]),
                ])

        write_tsv(
            PHYLOSEQ_OBJECT_INVENTORY,
            [
                "object_name",
                "object_class",
                "is_phyloseq_like",
                "sample_count",
                "taxa_count",
                "taxa_are_rows",
                "inspection_status",
                "errors",
            ],
            obj_rows,
        )
        write_tsv(PHYLOSEQ_SAMPLE_METADATA_INVENTORY, ["object_name", "sample_data_column", "sample_ids_preview"], sm_rows)
        write_tsv(PHYLOSEQ_TAXA_INVENTORY, ["object_name", "tax_rank_name", "taxa_names_preview"], tx_rows)

        meta_rows = inspect_metadata_files()
        write_tsv(
            METADATA_COLUMN_INVENTORY,
            [
                "file_name",
                "relative_path",
                "read_status",
                "row_count",
                "column_count",
                "first_30_columns",
                "likely_sample_id_columns",
                "likely_disease_status_columns",
                "likely_cohort_source_columns",
                "likely_covariate_columns",
                "missingness_summary_for_disease_status_columns",
            ],
            meta_rows,
        )

        comp_inv_rows, cand_rows, label_rows, tax_rows = component_rows_to_tables(components, meta_rows)

        write_tsv(
            PHYLOSEQ_COMPONENT_INVENTORY,
            [
                "component_path",
                "component_name",
                "component_class",
                "is_phyloseq",
                "sample_count",
                "taxa_count",
                "taxa_are_rows",
                "sample_data_n_columns",
                "sample_data_columns_preview",
                "tax_rank_names",
                "taxa_names_preview",
                "inspection_status",
                "notes",
            ],
            comp_inv_rows,
        )

        write_tsv(
            CANDIDATE_COMPONENT_SELECTION,
            [
                "component_path",
                "component_name",
                "is_phyloseq",
                "sample_count",
                "taxa_count",
                "has_sample_data",
                "has_otu_table",
                "has_tax_table",
                "has_species_rank",
                "has_genus_rank",
                "candidate_for_integrated_us",
                "reason",
                "blocking_issue",
                "notes",
            ],
            cand_rows,
        )

        write_tsv(
            LABEL_MAPPING_CANDIDATES,
            [
                "source_object_or_file",
                "component_path",
                "n_rows",
                "sample_id_column_candidate",
                "disease_label_column_candidate",
                "observed_label_values",
                "pd_label_candidate",
                "control_label_candidate",
                "cohort_column_candidate",
                "evidence_quality",
                "notes",
            ],
            label_rows,
        )

        write_tsv(
            TAXONOMY_RANK_READINESS,
            [
                "component_path",
                "tax_rank_names",
                "species_rank_present",
                "species_rank_column",
                "genus_rank_present",
                "genus_rank_column",
                "taxonomy_preview",
                "ready_for_species_export",
                "blocking_issue",
                "notes",
            ],
            tax_rows,
        )

        plan_rows, status = build_extraction_plan(cand_rows, tax_rows, meta_rows)
        write_tsv(
            EXTRACTION_PLAN,
            [
                "target_standard_file",
                "source_object_or_file",
                "source_component",
                "required_for_replication",
                "extraction_feasibility",
                "proposed_method",
                "blocking_issue",
                "notes",
            ],
            plan_rows,
        )

        # D5 conservative status (guard against list-only unresolved)
        sample_inv_empty = len(sm_rows) == 0
        taxa_inv_empty = len(tx_rows) == 0
        has_phy_component = any(r[2] == "true" for r in cand_rows)
        species_unclear = not any(r[2] == "true" and r[7] == "yes" for r in tax_rows)

        d5_status = status
        if (not has_phy_component) or sample_inv_empty or taxa_inv_empty or species_unclear:
            d5_status = "PARTIAL_READY_NEEDS_LIST_COMPONENT_INSPECTION"

        D5_STATUS.write_text(
            f"D5_STATUS = {d5_status}\n"
            "READY_FOR_PROCESSED_REPLICATION_ANALYSIS = NO\n"
            "READY_FOR_RAW_READ_VALIDATION = NO\n"
            f"PHYLOSEQ_INSPECTION_STATUS = {r_result.get('status', 'UNKNOWN')}\n"
            f"PHYLOSEQ_PACKAGE_AVAILABLE = {bool(r_result.get('phyloseq_available', False))}\n",
            encoding="utf-8",
        )

        D5B_STATUS.write_text(
            f"D5B_STATUS = {status}\n"
            "READY_FOR_PROCESSED_REPLICATION_ANALYSIS = NO\n"
            "READY_FOR_RAW_READ_VALIDATION = NO\n"
            f"PHYLOSEQ_COMPONENTS_TOTAL = {len(comp_inv_rows)}\n",
            encoding="utf-8",
        )

        D5B_REPORT.write_text(
            "# PHASE 2E-D5b Integrated-US Component Inspection Report\n\n"
            f"- RData top-level objects inspected: {len(obj_rows)}\n"
            f"- Recursive components inspected (depth<=3): {len(comp_inv_rows)}\n"
            f"- Candidate components (yes): {sum(1 for r in cand_rows if r[10]=='yes')}\n"
            f"- Candidate components (needs_manual_selection): {sum(1 for r in cand_rows if r[10]=='needs_manual_selection')}\n"
            f"- D5b status: {status}\n\n"
            "No FASTQ/SRA/raw-read inspection was performed.\n"
            "No Curli Carrier Burden or replication statistics were run.\n",
            encoding="utf-8",
        )
    except Exception as exc:
        D5_STATUS.write_text(f"D5_STATUS = INSPECTION_FAILED\nERROR = {type(exc).__name__}: {exc}\n", encoding="utf-8")
        D5B_STATUS.write_text(f"D5B_STATUS = INSPECTION_FAILED\nERROR = {type(exc).__name__}: {exc}\n", encoding="utf-8")
        raise


if __name__ == "__main__":
    main()
