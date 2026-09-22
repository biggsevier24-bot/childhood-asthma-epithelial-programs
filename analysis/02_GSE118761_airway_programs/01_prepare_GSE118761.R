options(stringsAsFactors = FALSE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
source(file.path(repo, "R/common/public_geo_inputs.R"))
use_codex_library_if_available()
source(file.path(repo, "R/common/module_scoring.R"))

expr_path <- Sys.getenv("GSE118761_EXPR", file.path(repo, "data/public/GSE118761/GSE118761_VST_gene_expression.csv"))
meta_path <- Sys.getenv("GSE118761_META", file.path(repo, "data/public/GSE118761/GSE118761_metadata.csv"))
raw_dir <- file.path(repo, "data/public/GSE118761/raw")
counts_path <- file.path(raw_dir, "GSE118761_genecounts.csv.gz")
series_path <- file.path(raw_dir, "GSE118761_series_matrix.txt.gz")

if (!file.exists(expr_path) || !file.exists(meta_path)) {
  if (!file.exists(counts_path) || !file.exists(series_path)) source(file.path(repo, "analysis/02_GSE118761_airway_programs/00_download_GSE118761.R"))
  sm <- read_series_sample_fields(series_path)
  metadata <- data.frame(
    sample_id = sub("^Sample name: ", "", sm$description),
    gsm_accession = sm$geo_accession,
    sample_title = sm$title,
    diagnosis = sub("^diagnosis: ", "", sm$characteristics_ch1),
    subject_id = sub("^subject_id: ", "", sm$characteristics_ch1.2),
    tissue = ifelse(grepl("^tissue: Nasal", sm$characteristics_ch1.1), "nasal", "tracheal"),
    baseline = TRUE,
    asthma = as.integer(grepl("asthmatic$", sub("^diagnosis: ", "", sm$characteristics_ch1)) & !grepl("non-asthmatic", sm$characteristics_ch1, ignore.case = TRUE)),
    atopy = as.integer(grepl("^Atopic", sub("^diagnosis: ", "", sm$characteristics_ch1))),
    wheeze = as.integer(grepl("wheezer", sm$title, ignore.case = TRUE) & !grepl("non-wheezer", sm$title, ignore.case = TRUE)),
    age = as.numeric(sub("^age \\(yrs\\): ", "", sm$characteristics_ch1.4)),
    sex = sub("^gender: ", "", sm$characteristics_ch1.3),
    batch = NA_character_, stringsAsFactors = FALSE
  )
  x <- read.csv(gzfile(counts_path), check.names = FALSE, stringsAsFactors = FALSE)
  sample_cols <- setdiff(names(x), c("EnsemblID", "Symbol", "Description"))
  counts <- as.matrix(x[, sample_cols, drop = FALSE])
  storage.mode(counts) <- "numeric"
  vst <- prepare_vst(counts, x$Symbol, "GSE118761")
  stopifnot(setequal(colnames(vst), metadata$sample_id))
  metadata <- metadata[match(colnames(vst), metadata$sample_id), , drop = FALSE]
  vst <- vst[, metadata$sample_id, drop = FALSE]
  dir.create(dirname(expr_path), recursive = TRUE, showWarnings = FALSE)
  write_expression_csv(vst, expr_path)
  write.csv(metadata, meta_path, row.names = FALSE)
}

a <- align_expression_metadata(read_gene_expression_csv(expr_path), read.csv(meta_path, check.names = FALSE))
stopifnot(nrow(a$metadata) == 104L, sum(a$metadata$tissue == "nasal") == 55L, sum(a$metadata$tissue == "tracheal") == 49L)
dir.create(file.path(repo, "outputs/02_GSE118761"), recursive = TRUE, showWarnings = FALSE)
write.csv(a$metadata, file.path(repo, "outputs/02_GSE118761/GSE118761_processed_metadata.csv"), row.names = FALSE)
saveRDS(a, file.path(repo, "outputs/02_GSE118761/GSE118761_aligned_input.rds"))
