# External ECM sensitivity entry point; it is independent of fixed-module derivation.
options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "R/09_external_ECM_sensitivity.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "data/external_sensitivity")
expression_file <- file.path(input_dir, "standardized_gene_expression.tsv")
gmt_file <- file.path(input_dir, "external_ECM_gene_sets.gmt")
if (!all(file.exists(c(expression_file, gmt_file)))) {
  stop("External sensitivity inputs are intentionally separate. See data/external_sensitivity/README.md.")
}
expression <- as.matrix(read.delim(expression_file, row.names = 1, check.names = FALSE))
gmt <- readLines(gmt_file, warn = FALSE)
gene_sets <- lapply(strsplit(gmt, "\t", fixed = TRUE), function(x) x[-c(1, 2)])
names(gene_sets) <- vapply(strsplit(gmt, "\t", fixed = TRUE), `[`, character(1), 1)
z <- t(scale(t(expression)))
z[!is.finite(z)] <- 0
scores <- t(vapply(gene_sets, function(genes) {
  covered <- intersect(genes, rownames(z))
  if (!length(covered)) return(rep(NA_real_, ncol(z)))
  colMeans(z[covered, , drop = FALSE])
}, numeric(ncol(z))))
colnames(scores) <- colnames(z)
out_dir <- file.path(repo_root, "outputs/external_ECM_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(scores, file.path(out_dir, "mean_z_scores.tsv"), sep = "\t", quote = FALSE, col.names = NA)
cat("External ECM mean-z sensitivity scores written. These gene sets do not select the fixed 32 genes.\n")
