source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "00_projection_common.R"))
projection_log("START 08_verify_projection_results.R")

required <- c(
  "provenance/input_manifest.tsv", "projection_sample_audit.tsv",
  "GSE118761_projection_program_scores.tsv", "projection_scaling_parameters.tsv",
  "GSE152004_reference_centroids.tsv", "GSE152004_centroid_source_samples.tsv",
  "GSE118761_projection_distances.tsv", "GSE118761_projection_assignments_FINAL.tsv",
  "GSE118761_projection_counts_FINAL.tsv", "GSE118761_nasal_projection_phenotypes_FINAL.tsv",
  "GSE118761_tracheal_projection_phenotypes_FINAL.tsv", "GSE118761_projection_confidence_FINAL.tsv"
)
stopifnot(all(file.exists(file.path(projection_out, required))))

assignments <- projection_read_tsv(file.path(projection_out, "GSE118761_projection_assignments_FINAL.tsv"))
counts <- projection_read_tsv(file.path(projection_out, "GSE118761_projection_counts_FINAL.tsv"))
distances <- projection_read_tsv(file.path(projection_out, "GSE118761_projection_distances.tsv"))
stopifnot(nrow(assignments) == 104L, !anyDuplicated(assignments$sample_id))
stopifnot(sum(assignments$tissue == "nasal") == 55L, sum(assignments$tissue == "tracheal") == 49L)
stopifnot(identical(assignments$sample_id, distances$sample_id))
rederived <- c("E1", "E2", "E3")[max.col(-as.matrix(distances[, c("distance_to_E1", "distance_to_E2", "distance_to_E3")]), ties.method = "first")]
stopifnot(identical(assignments$assigned_state, rederived))
for (i in seq_len(nrow(counts))) stopifnot(counts$n[i] == projection_expected_count(counts$tissue[i], counts$state[i]))

hash_targets <- c(
  config = projection_config_path,
  reference_centroids = file.path(projection_out, "GSE152004_reference_centroids.tsv"),
  sample_assignments = file.path(projection_out, "GSE118761_projection_assignments_FINAL.tsv"),
  projection_counts = file.path(projection_out, "GSE118761_projection_counts_FINAL.tsv"),
  nasal_phenotypes = file.path(projection_out, "GSE118761_nasal_projection_phenotypes_FINAL.tsv"),
  tracheal_phenotypes = file.path(projection_out, "GSE118761_tracheal_projection_phenotypes_FINAL.tsv"),
  Figure5f_source = file.path(projection_figure_dir, "Figure5f_source_data.tsv"),
  S9_source = file.path(projection_table_dir, "Supplementary_Table_S9_source.tsv")
)
repo_prefix <- paste0(projection_repo, "/")
relative_paths <- vapply(hash_targets, function(path) {
  path <- gsub("\\\\", "/", path)
  if (startsWith(path, repo_prefix)) substring(path, nchar(repo_prefix) + 1L) else path
}, character(1))
hashes <- data.frame(artifact = names(hash_targets), relative_path = relative_paths,
                     sha256 = vapply(hash_targets, projection_sha256, character(1)))
projection_write_tsv(hashes, file.path(projection_out, "PROJECTION_OUTPUT_SHA256.tsv"))
writeLines(projection_sha256(projection_config_path), file.path(projection_repo, "config/GSE118761_projection_final.sha256"))

verification <- data.frame(
  check = c("total_samples", "nasal_samples", "tracheal_samples", "unique_sample_id", "target_counts", "distance_assignment_recomputed", "all_required_outputs"),
  observed = c(nrow(assignments), sum(assignments$tissue == "nasal"), sum(assignments$tissue == "tracheal"), !anyDuplicated(assignments$sample_id), TRUE, TRUE, TRUE),
  expected = c(104, 55, 49, TRUE, TRUE, TRUE, TRUE),
  PASS = TRUE
)
projection_write_tsv(verification, file.path(projection_out, "projection_run_verification.tsv"))
writeLines("PROJECTION_CLEAN_RERUN_PASS = TRUE", file.path(projection_out, "PROJECTION_CLEAN_RERUN_PASS.txt"))

projection_log("SESSION_INFO_BEGIN")
cat(paste(capture.output(sessionInfo()), collapse = "\n"), "\n", file = projection_log_path, append = TRUE)
projection_log("SESSION_INFO_END")
projection_log("END 08_verify_projection_results.R PASS")
