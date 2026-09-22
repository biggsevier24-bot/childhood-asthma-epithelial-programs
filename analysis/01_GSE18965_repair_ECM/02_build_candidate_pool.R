build_repair_ecm_candidates <- function(gene_expression, metadata, membership_file) {
  membership <- read.delim(membership_file, check.names = FALSE, stringsAsFactors = FALSE)
  go_columns <- c("ECM_organization", "collagen_fibril_organization", "wound_healing", "cell_adhesion")
  family_columns <- c("collagen", "integrin", "laminin", "MMP", "matricellular")
  annotation <- data.frame(gene = rownames(gene_expression), stringsAsFactors = FALSE)
  for (column in c(go_columns, family_columns)) {
    annotation[[column]] <- annotation$gene %in% membership$gene_symbol[membership$source_name == column]
  }
  annotation$GO_any <- rowSums(annotation[, go_columns, drop = FALSE]) > 0
  annotation$family_any <- rowSums(annotation[, family_columns, drop = FALSE]) > 0
  annotation$annotation_count <- rowSums(annotation[, c(go_columns, family_columns), drop = FALSE])
  aa <- metadata$gsm[metadata$group == "AA"]
  hn <- metadata$gsm[metadata$group == "HN"]
  annotation$mean_AA <- rowMeans(gene_expression[annotation$gene, aa, drop = FALSE])
  annotation$mean_HN <- rowMeans(gene_expression[annotation$gene, hn, drop = FALSE])
  annotation$mean_AA_minus_HN <- annotation$mean_AA - annotation$mean_HN
  annotation$AA_gt_HN <- annotation$mean_AA_minus_HN > 0
  rownames(annotation) <- annotation$gene
  candidate <- sort(annotation$gene[
    annotation$GO_any & (annotation$family_any | annotation$collagen_fibril_organization)
  ])
  stopifnot(length(candidate) == 140L)
  list(annotation = annotation, candidates = candidate,
       go_columns = go_columns, family_columns = family_columns)
}
