
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
