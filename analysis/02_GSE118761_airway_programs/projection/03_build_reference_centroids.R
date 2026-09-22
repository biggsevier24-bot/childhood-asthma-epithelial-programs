source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "00_projection_common.R"))
projection_log("START 03_build_reference_centroids.R")

inp <- readRDS(file.path(projection_out, "projection_inputs.rds"))
required_scores <- c("sample_id", "T2_z", "IFN_z", "repair_ECM_z")
if (!all(required_scores %in% names(inp$scores152))) stop("GSE152004 score columns missing")
if (!all(c("sample_id", "state_label") %in% names(inp$assign152))) stop("GSE152004 assignment columns missing")

idx <- match(inp$scores152$sample_id, inp$assign152$sample_id)
if (anyNA(idx)) stop("GSE152004 assignment join failed")
source_samples <- data.frame(
  sample_id = inp$scores152$sample_id,
  state = projection_short_state(inp$assign152$state_label[idx]),
  T2 = inp$scores152$T2_z,
  IFN = inp$scores152$IFN_z,
  repair = inp$scores152$repair_ECM_z
)
state_order <- c("E1", "E2", "E3")
source_samples$state <- factor(source_samples$state, levels = state_order)
centroids <- do.call(rbind, lapply(state_order, function(state) {
  d <- source_samples[source_samples$state == state, ]
  data.frame(state = state, T2 = mean(d$T2), IFN = mean(d$IFN), repair = mean(d$repair),
             n_reference_samples = nrow(d),
             source_assignment_file = projection_config$inputs$gse152_assignments)
}))
stopifnot(sum(centroids$n_reference_samples) == 695L)
projection_write_tsv(source_samples, file.path(projection_out, "GSE152004_centroid_source_samples.tsv"))
projection_write_tsv(centroids, file.path(projection_out, "GSE152004_reference_centroids.tsv"))
for (i in seq_len(nrow(centroids))) {
  projection_log("CENTROID", centroids$state[i], "n=", centroids$n_reference_samples[i],
                 "T2=", signif(centroids$T2[i], 12), "IFN=", signif(centroids$IFN[i], 12),
                 "repair=", signif(centroids$repair[i], 12))
}
projection_log("END 03_build_reference_centroids.R")
