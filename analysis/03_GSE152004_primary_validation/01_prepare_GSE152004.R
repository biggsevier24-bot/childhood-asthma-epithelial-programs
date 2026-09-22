options(stringsAsFactors = FALSE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
source(file.path(repo, "R/common/public_geo_inputs.R"))
use_codex_library_if_available()
source(file.path(repo, "R/common/module_scoring.R"))

expr_path <- Sys.getenv("GSE152004_EXPR", file.path(repo, "data/public/GSE152004/GSE152004_VST_gene_expression.csv"))
meta_path <- Sys.getenv("GSE152004_META", file.path(repo, "data/public/GSE152004/GSE152004_metadata.csv"))
raw_dir <- file.path(repo, "data/public/GSE152004/raw")
counts_path <- file.path(raw_dir, "GSE152004_695_raw_counts.txt.gz")
series_path <- file.path(raw_dir, "GSE152004_series_matrix.txt.gz")

if (!file.exists(expr_path) || !file.exists(meta_path)) {
  if (!file.exists(counts_path) || !file.exists(series_path)) source(file.path(repo, "analysis/03_GSE152004_primary_validation/00_download_GSE152004.R"))
  sm <- read_series_sample_fields(series_path)
  metadata <- data.frame(
    sample_id = sm$title,
    gsm_accession = sm$geo_accession,
    subject_id = sm$title,
    tissue = "nasal",
    asthma = as.integer(sub("^asthma status: ", "", sm$characteristics_ch1.1) == "asthmatic"),
    atopy = NA_integer_, age = NA_real_, sex = NA_character_, batch = NA_character_,
    stringsAsFactors = FALSE
  )
  x <- read.delim(gzfile(counts_path), check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
  counts <- as.matrix(x)
  storage.mode(counts) <- "numeric"
  vst <- prepare_vst(counts, rownames(counts), "GSE152004")
  stopifnot(setequal(colnames(vst), metadata$sample_id))
  vst <- vst[, metadata$sample_id, drop = FALSE]
  dir.create(dirname(expr_path), recursive = TRUE, showWarnings = FALSE)
  write_expression_csv(vst, expr_path)
  write.csv(metadata, meta_path, row.names = FALSE)
}

a <- align_expression_metadata(read_gene_expression_csv(expr_path), read.csv(meta_path, check.names = FALSE))
stopifnot(nrow(a$metadata) == 695L)
dir.create(file.path(repo, "outputs/03_GSE152004"), recursive = TRUE, showWarnings = FALSE)
write.csv(a$metadata, file.path(repo, "outputs/03_GSE152004/GSE152004_processed_metadata.csv"), row.names = FALSE)
saveRDS(a, file.path(repo, "outputs/03_GSE152004/GSE152004_aligned_input.rds"))
