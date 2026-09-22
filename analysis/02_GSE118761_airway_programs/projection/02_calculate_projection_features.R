source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "00_projection_common.R"))
projection_log("START 02_calculate_projection_features.R")

inp <- readRDS(file.path(projection_out, "projection_inputs.rds"))
sets <- list(
  T2 = projection_read_gene_set(inp$paths$t2_gene_set),
  IFN = projection_read_gene_set(inp$paths$ifn_gene_set),
  repair = projection_read_gene_set(inp$paths$repair_gene_set)
)
raw_scores <- lapply(sets, function(g) projection_mean_gene_z(inp$expr118, g))

scores <- inp$meta118
scores$T2_score_raw <- raw_scores$T2[scores$sample_id]
scores$IFN_score_raw <- raw_scores$IFN[scores$sample_id]
scores$repair_score_raw <- raw_scores$repair[scores$sample_id]
scores$T2_projection <- NA_real_
scores$IFN_projection <- NA_real_
scores$repair_projection <- NA_real_

scaling_rows <- list()
for (tissue in c("nasal", "tracheal")) {
  idx <- which(scores$tissue == tissue)
  for (feature in c("T2", "IFN", "repair")) {
    raw_col <- paste0(feature, "_score_raw")
    projection_col <- paste0(feature, "_projection")
    mu <- mean(scores[[raw_col]][idx])
    sigma <- sd(scores[[raw_col]][idx])
    if (!is.finite(sigma) || sigma <= 0) stop("Non-positive projection SD")
    scores[[projection_col]][idx] <- (scores[[raw_col]][idx] - mu) / sigma
    scaling_rows[[length(scaling_rows) + 1L]] <- data.frame(
      feature = feature,
      reference_population = paste0("GSE118761_", tissue, "_samples"),
      mean = mu,
      sd = sigma,
      transformation = "(raw_score - tissue_mean) / tissue_sample_sd"
    )
  }
}
if (anyNA(scores[, c("T2_projection", "IFN_projection", "repair_projection")])) stop("Missing projection coordinate")
projection_write_tsv(scores, file.path(projection_out, "GSE118761_projection_program_scores.tsv"))
projection_write_tsv(do.call(rbind, scaling_rows), file.path(projection_out, "projection_scaling_parameters.tsv"))

coverage <- do.call(rbind, lapply(names(sets), function(nm) data.frame(
  module = nm,
  requested_n = length(sets[[nm]]),
  detected_n = length(attr(raw_scores[[nm]], "genes_present")),
  detected_genes = paste(attr(raw_scores[[nm]], "genes_present"), collapse = ";")
)))
projection_write_tsv(coverage, file.path(projection_out, "projection_gene_coverage.tsv"))
projection_log("SCORING gene-wise z across all 104 samples; projection z within tissue")
projection_log("END 02_calculate_projection_features.R")
