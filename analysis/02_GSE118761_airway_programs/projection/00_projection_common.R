codex_r_helper <- Sys.getenv("CODEX_R_LIB_HELPER", "")
if (nzchar(codex_r_helper) && file.exists(codex_r_helper)) source(codex_r_helper)
options(stringsAsFactors = FALSE)

projection_script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
projection_repo <- Sys.getenv("PROJECTION_REPO_ROOT", "")
if (!nzchar(projection_repo)) {
  projection_repo <- file.path(dirname(projection_script), "../../..")
}
projection_repo <- normalizePath(projection_repo, winslash = "/", mustWork = TRUE)

if (!requireNamespace("yaml", quietly = TRUE)) stop("Package yaml is required")
if (!requireNamespace("digest", quietly = TRUE)) stop("Package digest is required")

projection_config_path <- file.path(projection_repo, "config/GSE118761_projection_final.yaml")
projection_config <- yaml::read_yaml(projection_config_path)
projection_out <- file.path(projection_repo, projection_config$outputs$projection_dir)
projection_figure_dir <- file.path(projection_repo, projection_config$outputs$figure_dir)
projection_table_dir <- file.path(projection_repo, projection_config$outputs$table_dir)
projection_audit_dir <- file.path(projection_repo, projection_config$outputs$audit_dir)
projection_docs_dir <- file.path(projection_repo, projection_config$outputs$docs_dir)
projection_log_path <- file.path(projection_repo, projection_config$outputs$log_file)

for (d in c(projection_out, file.path(projection_out, "provenance"), projection_figure_dir,
            projection_table_dir, projection_audit_dir, projection_docs_dir,
            dirname(projection_log_path))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

projection_log <- function(...) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), paste(..., collapse = " "))
  cat(line, "\n", file = projection_log_path, append = TRUE)
  message(line)
}

projection_write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

projection_read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

projection_sha256 <- function(path) {
  if (!file.exists(path)) stop("Missing file for SHA256: ", path)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

projection_input_path <- function(config_key, env_name) {
  env_value <- Sys.getenv(env_name, "")
  path <- if (nzchar(env_value)) env_value else file.path(projection_repo, projection_config$inputs[[config_key]])
  if (!file.exists(path)) stop("Missing required input ", config_key, ": ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

projection_short_state <- function(x) {
  out <- sub("_.*$", "", as.character(x))
  if (any(!out %in% c("E1", "E2", "E3"))) stop("Unexpected state label")
  out
}

projection_read_expression <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  genes <- toupper(trimws(as.character(x[[1]])))
  x[[1]] <- NULL
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  rownames(m) <- genes
  if (anyDuplicated(genes)) {
    mu <- rowMeans(m, na.rm = TRUE)
    keep <- unlist(tapply(seq_along(genes), genes, function(i) i[which.max(mu[i])]), use.names = FALSE)
    m <- m[sort(keep), , drop = FALSE]
  }
  m
}

projection_read_gene_set <- function(path) {
  x <- readLines(path, warn = FALSE)
  unique(toupper(trimws(x[nzchar(trimws(x))])))
}

projection_mean_gene_z <- function(expr, genes) {
  present <- intersect(unique(toupper(genes)), rownames(expr))
  if (!length(present)) stop("No requested genes are present")
  m <- expr[present, , drop = FALSE]
  keep <- apply(m, 1, sd, na.rm = TRUE) > 0
  m <- m[keep, , drop = FALSE]
  z <- t(scale(t(m)))
  ans <- colMeans(z, na.rm = TRUE)
  attr(ans, "genes_present") <- rownames(z)
  ans
}

projection_expected_count <- function(tissue, state) {
  as.integer(projection_config$expected_projection_counts[[tissue]][[state]])
}
