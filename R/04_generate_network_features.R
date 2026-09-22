select_b1_anchor_candidates <- function(candidates, gene_expression, annotation, metadata,
                                        absolute_correlation, go_columns, family_columns) {
  seeds <- sort(intersect(
    annotation$gene[annotation$family_any & annotation$GO_any],
    annotation$gene[annotation$AA_gt_HN]
  ))
  stopifnot(length(seeds) == 73L)
  hc <- hclust(as.dist(1 - absolute_correlation[candidates, candidates, drop = FALSE]), method = "average")
  labels <- cutree(hc, k = 2L)
  metrics <- lapply(sort(unique(labels)), function(cluster_id) {
    genes <- names(labels)[labels == cluster_id]
    block <- absolute_correlation[genes, genes, drop = FALSE]
    within <- block[upper.tri(block)]
    score <- colMeans(z_rows(gene_expression[genes, , drop = FALSE]))
    a <- annotation[genes, , drop = FALSE]
    data.frame(
      cluster_id = cluster_id, cluster_n = length(genes),
      annotation_density = mean(a$annotation_count / length(c(go_columns, family_columns))),
      mean_within_Pearson = if (length(within)) mean(within) else 1,
      module_AA_minus_HN = mean(score[metadata$group == "AA"]) - mean(score[metadata$group == "HN"]),
      stringsAsFactors = FALSE
    )
  })
  metrics <- do.call(rbind, metrics)
  eligible <- metrics[metrics$cluster_n >= 10L & metrics$module_AA_minus_HN > 0, , drop = FALSE]
  eligible <- eligible[order(-eligible$annotation_density, -eligible$mean_within_Pearson,
                             -eligible$module_AA_minus_HN, -eligible$cluster_n, eligible$cluster_id), , drop = FALSE]
  stopifnot(nrow(eligible) > 0L)
  selected_id <- eligible$cluster_id[1]
  initial <- sort(names(labels)[labels == selected_id])
  stopifnot(length(initial) == 138L)
  seed_block <- absolute_correlation[initial, seeds, drop = FALSE]
  for (i in seq_along(initial)) {
    if (initial[i] %in% seeds) seed_block[i, initial[i]] <- NA_real_
  }
  seed_degree <- rowMeans(seed_block >= 0.30, na.rm = TRUE)
  retained <- sort(initial[is.finite(seed_degree) & seed_degree >= 0.30])
  stopifnot(length(retained) == 115L)
  audit <- data.frame(gene = initial, seed_degree_abs_r030 = seed_degree,
                      retained = initial %in% retained, stringsAsFactors = FALSE)
  list(seeds = seeds, labels = labels, metrics = metrics, initial = initial,
       retained = retained, audit = audit)
}

add_seed_feature <- function(feature, signed, absolute, annotation, candidates,
                             seed_name, seed_genes, metric, correlation_type) {
  matrix_used <- if (correlation_type == "signed") signed else absolute
  values <- vapply(feature$gene, function(gene) {
    peers <- setdiff(seed_genes, gene)
    x <- matrix_used[gene, peers]
    if (metric == "mean") mean(x) else if (metric == "median") median(x) else max(x)
  }, numeric(1))
  feature[[paste(tolower(seed_name), correlation_type, metric, sep = "_")]] <- values
  feature
}

generate_network_features <- function(expression_result, candidates, gene_expression,
                                      annotation, metadata, go_columns, family_columns) {
  signed <- expression_result$signed_correlation
  absolute <- expression_result$absolute_correlation
  anchor <- select_b1_anchor_candidates(candidates, gene_expression, annotation, metadata,
                                        absolute, go_columns, family_columns)
  feature <- expression_result$feature
  seed_definitions <- list(
    AA_POSITIVE_ALL = candidates[annotation[candidates, "mean_AA_minus_HN"] > 0],
    AA_POSITIVE_GO_COLLAGEN = candidates[annotation[candidates, "mean_AA_minus_HN"] > 0 & annotation[candidates, "collagen_fibril_organization"]],
    AA_POSITIVE_INTEGRIN_FAMILY = candidates[annotation[candidates, "mean_AA_minus_HN"] > 0 & annotation[candidates, "integrin"]],
    AA_POSITIVE_MMP_FAMILY = candidates[annotation[candidates, "mean_AA_minus_HN"] > 0 & annotation[candidates, "MMP"]],
    AA_POSITIVE_MATRICELLULAR_FAMILY = candidates[annotation[candidates, "mean_AA_minus_HN"] > 0 & annotation[candidates, "matricellular"]]
  )
  feature <- add_seed_feature(feature, signed, absolute, annotation, candidates,
                              "AA_POSITIVE_INTEGRIN_FAMILY", seed_definitions$AA_POSITIVE_INTEGRIN_FAMILY,
                              "median", "absolute")
  feature <- add_seed_feature(feature, signed, absolute, annotation, candidates,
                              "AA_POSITIVE_MATRICELLULAR_FAMILY", seed_definitions$AA_POSITIVE_MATRICELLULAR_FAMILY,
                              "mean", "signed")
  feature <- add_seed_feature(feature, signed, absolute, annotation, candidates,
                              "AA_POSITIVE_MMP_FAMILY", seed_definitions$AA_POSITIVE_MMP_FAMILY,
                              "max", "signed")
  feature <- add_seed_feature(feature, signed, absolute, annotation, candidates,
                              "AA_POSITIVE_ALL", seed_definitions$AA_POSITIVE_ALL,
                              "median", "signed")
  collagen_seeds <- seed_definitions$AA_POSITIVE_GO_COLLAGEN
  feature$aa_positive_go_collagen_signed_degree_060 <- vapply(feature$gene, function(gene) {
    peers <- setdiff(collagen_seeds, gene)
    mean(signed[gene, peers] >= 0.60)
  }, numeric(1))
  feature <- feature[match(anchor$retained, feature$gene), , drop = FALSE]
  stopifnot(all(feature$gene == anchor$retained), nrow(feature) == 115L)
  list(feature = feature, anchor = anchor, signed_correlation = signed,
       absolute_correlation = absolute, seed_definitions = seed_definitions)
}
