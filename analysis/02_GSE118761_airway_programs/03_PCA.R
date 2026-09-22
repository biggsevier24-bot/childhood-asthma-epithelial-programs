options(stringsAsFactors = FALSE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)

a <- readRDS(file.path(repo, "outputs/02_GSE118761/GSE118761_aligned_input.rds"))
modules <- read.csv(file.path(repo, "outputs/02_GSE118761/GSE118761_recomputed_module_scores.csv"), check.names = FALSE)
source_dir <- file.path(repo, "figure_source_data/Figure4")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

variables <- c(
  T2 = "T2_z",
  IFN = "IFN_z",
  repair_ECM = "repair_ECM_z",
  atopy = "atopy",
  asthma = "asthma",
  wheeze = "wheeze"
)

for (tissue in c("nasal", "tracheal")) {
  idx <- which(tolower(a$metadata$tissue) == tissue)
  stopifnot(length(idx) > 10L)

  # Historical Figure 4 implementation: PCA was performed separately by tissue
  # using the top 2,000 genes ranked by within-tissue variance, with both
  # centering and per-gene scaling enabled.
  tissue_expr <- a$expr[, idx, drop = FALSE]
  gene_var <- apply(tissue_expr, 1L, var, na.rm = TRUE)
  eligible <- which(is.finite(gene_var) & gene_var > 0)
  stopifnot(length(eligible) >= 2000L)
  ord <- eligible[order(gene_var[eligible], decreasing = TRUE)]
  top <- ord[seq_len(2000L)]

  p <- prcomp(t(tissue_expr[top, , drop = FALSE]), center = TRUE, scale. = TRUE, rank. = 10L)
  npc <- min(10L, ncol(p$x))
  scores <- data.frame(
    sample_id = a$metadata$sample_id[idx],
    tissue = tissue,
    p$x[, seq_len(npc), drop = FALSE],
    check.names = FALSE
  )

  matched <- modules[match(scores$sample_id, modules$sample_id), , drop = FALSE]
  stopifnot(!anyNA(matched$sample_id), identical(as.character(scores$sample_id), as.character(matched$sample_id)))

  rows <- do.call(rbind, lapply(seq_len(npc), function(j) {
    do.call(rbind, lapply(names(variables), function(v) {
      y <- matched[[variables[[v]]]]
      x <- p$x[, j]
      ok <- is.finite(x) & is.finite(y)
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      data.frame(
        PC = paste0("PC", j),
        variable = v,
        rho = unname(ct$estimate),
        P = unname(ct$p.value),
        n = sum(ok),
        stringsAsFactors = FALSE
      )
    }))
  }))

  top_genes <- data.frame(
    rank = seq_along(top),
    gene = rownames(a$expr)[top],
    variance = gene_var[top],
    stringsAsFactors = FALSE
  )

  out_dir <- file.path(repo, "outputs/02_GSE118761")
  write.csv(scores, file.path(out_dir, paste0("GSE118761_", tissue, "_PCA_scores_recomputed.csv")), row.names = FALSE)
  write.csv(rows, file.path(out_dir, paste0("GSE118761_", tissue, "_PCA_correlations.csv")), row.names = FALSE)
  write.csv(top_genes, file.path(out_dir, paste0("GSE118761_", tissue, "_PCA_top2000_genes.csv")), row.names = FALSE)

  # Figure 4 source data are generated directly from the canonical analysis.
  write.csv(scores, file.path(source_dir, paste0("GSE118761_", tissue, "_PCA_scores.csv")), row.names = FALSE)
  write.csv(rows, file.path(source_dir, paste0("GSE118761_", tissue, "_PCA_correlations.csv")), row.names = FALSE)
}
