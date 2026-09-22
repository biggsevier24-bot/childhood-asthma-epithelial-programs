z_rows <- function(x) {
  z <- t(scale(t(x)))
  z[!is.finite(z)] <- 0
  z
}

cohen_d <- function(values, groups) {
  x <- values[groups == "AA"]
  y <- values[groups == "HN"]
  pooled <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) /
                   (length(x) + length(y) - 2))
  if (!is.finite(pooled) || pooled == 0) return(NA_real_)
  (mean(x) - mean(y)) / pooled
}

select_fixed_feature_cluster <- function(candidates, gene_expression, annotation, metadata,
                                         signed_correlation, absolute_correlation,
                                         go_columns, family_columns) {
  hc <- hclust(as.dist(1 - absolute_correlation[candidates, candidates, drop = FALSE]), method = "average")
  labels <- cutree(hc, h = 0.65)
  metrics <- lapply(sort(unique(labels)), function(cluster_id) {
    genes <- names(labels)[labels == cluster_id]
    block <- absolute_correlation[genes, genes, drop = FALSE]
    within <- if (length(genes) > 1L) block[upper.tri(block)] else NA_real_
    score <- colMeans(z_rows(gene_expression[genes, , drop = FALSE]))
    a <- annotation[genes, , drop = FALSE]
    data.frame(
      cluster_id = cluster_id,
      cluster_n = length(genes),
      module_score_AA_minus_HN = mean(score[metadata$group == "AA"]) - mean(score[metadata$group == "HN"]),
      module_score_cohen_d = cohen_d(score, metadata$group),
      fraction_AA_gt_HN = mean(a$AA_gt_HN),
      median_within_association = if (length(genes) > 1L) median(within, na.rm = TRUE) else NA_real_,
      curated_family_fraction = mean(a$family_any),
      category_breadth = sum(colSums(a[, c(go_columns, family_columns), drop = FALSE]) > 0),
      stringsAsFactors = FALSE
    )
  })
  metrics <- do.call(rbind, metrics)
  metrics$eligible <- metrics$cluster_n >= 3L & metrics$module_score_AA_minus_HN > 0 &
    metrics$fraction_AA_gt_HN >= 0.50 & metrics$median_within_association >= 0.60
  eligible <- metrics[metrics$eligible, , drop = FALSE]
  eligible <- eligible[order(-eligible$category_breadth, -eligible$curated_family_fraction,
                             -eligible$module_score_cohen_d, -eligible$median_within_association,
                             -eligible$cluster_n, eligible$cluster_id), , drop = FALSE]
  stopifnot(nrow(eligible) > 0L)
  selected_id <- eligible$cluster_id[1]
  selected <- sort(names(labels)[labels == selected_id])
  stopifnot(length(selected) == 97L)
  list(labels = labels, metrics = metrics, selected_genes = selected)
}

generate_expression_features <- function(candidates, gene_expression, annotation, metadata,
                                         go_columns, family_columns) {
  signed <- cor(t(gene_expression[candidates, , drop = FALSE]), method = "pearson", use = "pairwise.complete.obs")
  absolute <- abs(signed)
  selected <- select_fixed_feature_cluster(candidates, gene_expression, annotation, metadata,
                                            signed, absolute, go_columns, family_columns)
  feature <- annotation[candidates, c("gene", go_columns, family_columns, "annotation_count",
                                      "mean_AA", "mean_HN", "mean_AA_minus_HN"), drop = FALSE]
  names(feature)[match(go_columns, names(feature))] <-
    c("GO_ECM", "GO_collagen_fibril", "GO_wound_healing", "GO_cell_adhesion")
  names(feature)[match(family_columns, names(feature))] <- paste0("family_", family_columns)
  names(feature)[names(feature) == "mean_AA_minus_HN"] <- "AA_minus_HN"
  feature$mean_expression <- rowMeans(gene_expression[feature$gene, , drop = FALSE])
  feature$expression_variance <- apply(gene_expression[feature$gene, , drop = FALSE], 1L, var)
  selected_genes <- selected$selected_genes
  feature$mean_signed_correlation_to_selected_cluster <-
    rowMeans(signed[feature$gene, selected_genes, drop = FALSE])
  for (threshold in c(0.30, 0.40, 0.50, 0.60)) {
    suffix <- sprintf("%03d", round(threshold * 100))
    feature[[paste0("signed_degree_", suffix)]] <- vapply(feature$gene, function(gene) {
      peers <- setdiff(selected_genes, gene)
      mean(signed[gene, peers] >= threshold)
    }, numeric(1))
  }
  list(feature = feature, signed_correlation = signed, absolute_correlation = absolute,
       feature_cluster = selected)
}
