args <- commandArgs(trailingOnly=TRUE)
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
