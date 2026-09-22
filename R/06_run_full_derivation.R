options(stringsAsFactors = FALSE)
suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"))
library_helper <- Sys.getenv("CODEX_R_LIB_HELPER", unset = "")
if (nzchar(library_helper) && file.exists(library_helper)) source(library_helper)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "R/06_run_full_derivation.R"
default_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(Sys.getenv("REPAIR_ECM_REPO_ROOT", unset = default_root), winslash = "/", mustWork = TRUE)
for (script in sprintf("R/%02d_%s.R", 0:5, c(
  "prepare_GSE18965", "map_GPL96_probes", "build_candidate_pool",
  "generate_expression_features", "generate_network_features", "repair_ECM_final_selection_rule"
))) source(file.path(repo_root, script), local = FALSE)

output_dir <- file.path(repo_root, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
series_file <- file.path(repo_root, "data/GSE18965/GSE18965_series_matrix.txt.gz")
gpl_file <- file.path(repo_root, "data/GSE18965/GPL96-57554.txt")
membership_file <- file.path(repo_root, "data/annotation_sources/annotation_gene_membership.tsv")
stopifnot(file.exists(series_file), file.exists(gpl_file), file.exists(membership_file))

log_file <- file.path(output_dir, "derivation_run_log.txt")
if (file.exists(log_file)) file.remove(log_file)
log_message <- function(...) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = ""), sep = " | ")
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_message("Reading deposited GSE18965 Series Matrix")
series <- read_gse18965_series(series_file)
log_message("Parsing GPL96 and retaining highest-mean probe per gene")
mapped <- map_gpl96_highest_mean_probe(gpl_file, series$probe_expression)
log_message("Constructing 140-gene biologically annotated candidate pool")
candidate <- build_repair_ecm_candidates(mapped$gene_expression, series$metadata, membership_file)
log_message("Generating expression and Pearson features")
expression_result <- generate_expression_features(
  candidate$candidates, mapped$gene_expression, candidate$annotation, series$metadata,
  candidate$go_columns, candidate$family_columns
)
log_message("Applying frozen network prefilter and generating rule features")
network_result <- generate_network_features(
  expression_result, candidate$candidates, mapped$gene_expression,
  candidate$annotation, series$metadata, candidate$go_columns, candidate$family_columns
)
log_message("Applying frozen explicit if/else selection rule")
selected_mask <- select_repair_ECM_gene(network_result$feature)
final_genes <- sort(network_result$feature$gene[selected_mask])

write.table(series$metadata, file.path(output_dir, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mapped$selected_map, file.path(output_dir, "selected_highest_mean_probe_per_gene.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mapped$gene_expression, file.path(output_dir, "GSE18965_gene_expression.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(candidate$annotation[candidate$candidates, , drop = FALSE], file.path(output_dir, "candidate_pool.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(network_result$feature, file.path(output_dir, "gene_feature_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(network_result$signed_correlation, file.path(output_dir, "Pearson_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(network_result$anchor$audit, file.path(output_dir, "network_prefilter_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(expression_result$feature_cluster$metrics, file.path(output_dir, "feature_cluster_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(network_result$anchor$metrics, file.path(output_dir, "anchor_cluster_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene = final_genes), file.path(output_dir, "final_derived_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

stopifnot(length(final_genes) == 32L)
log_message("Derivation complete: FINAL_N=32")
cat("FINAL_N = ", length(final_genes), "\n", sep = "")
