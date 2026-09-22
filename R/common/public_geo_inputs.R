options(stringsAsFactors = FALSE)

use_codex_library_if_available <- function() {
  helper <- Sys.getenv("CODEX_R_LIB_HELPER", "")
  if (nzchar(helper) && file.exists(helper)) source(helper)
}

download_or_copy_geo <- function(url, destination, cache_env = "") {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(destination) && file.info(destination)$size > 0) return(destination)
  cached <- if (nzchar(cache_env)) Sys.getenv(cache_env, "") else ""
  if (nzchar(cached)) {
    if (!file.exists(cached)) stop(cache_env, " points to a missing file: ", cached)
    if (!file.copy(cached, destination, overwrite = TRUE)) stop("Could not copy cached GEO file")
  } else {
    download.file(url, destination, mode = "wb", quiet = FALSE)
  }
  if (!file.exists(destination) || file.info(destination)$size == 0) stop("Empty GEO input: ", destination)
  destination
}

read_series_sample_fields <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[grepl("^!Sample_", lines)]
  keys <- sub("\t.*$", "", lines)
  vals <- strsplit(sub("^[^\t]+\t", "", lines), "\t", fixed = TRUE)
  vals <- lapply(vals, function(x) gsub('^"|"$', "", x))
  n <- max(lengths(vals))
  out <- as.data.frame(matrix(NA_character_, n, length(vals)), stringsAsFactors = FALSE)
  names(out) <- make.unique(sub("^!Sample_", "", keys))
  for (i in seq_along(vals)) out[seq_along(vals[[i]]), i] <- vals[[i]]
  out
}

prepare_vst <- function(counts, symbols, label, min_nonzero_fraction = 0.20) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Bioconductor package DESeq2 is required")
  symbols <- toupper(trimws(as.character(symbols)))
  keep_symbol <- !is.na(symbols) & nzchar(symbols) & symbols != "NA" & !grepl("///", symbols, fixed = TRUE)
  counts <- counts[keep_symbol, , drop = FALSE]
  symbols <- symbols[keep_symbol]
  keep_expression <- rowSums(counts > 0, na.rm = TRUE) >= ceiling(min_nonzero_fraction * ncol(counts))
  counts <- counts[keep_expression, , drop = FALSE]
  symbols <- symbols[keep_expression]
  counts <- rowsum(counts, group = symbols, reorder = FALSE)
  counts <- round(as.matrix(counts))
  storage.mode(counts) <- "integer"
  coldata <- data.frame(sample_id = colnames(counts), row.names = colnames(counts))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~1)
  dds <- DESeq2::estimateSizeFactors(dds)
  vst <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))
  if (anyNA(vst)) stop("NA values in VST matrix for ", label)
  vst
}

write_expression_csv <- function(matrix, path) {
  out <- data.frame(gene_symbol = rownames(matrix), matrix, check.names = FALSE)
  write.csv(out, path, row.names = FALSE)
}
