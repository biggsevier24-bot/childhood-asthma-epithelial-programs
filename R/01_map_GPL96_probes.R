map_gpl96_highest_mean_probe <- function(gpl_file, probe_expression) {
  gpl <- read.delim(gpl_file, sep = "\t", comment.char = "#", quote = "", check.names = FALSE)
  probe_map <- data.frame(
    probe_id = gpl$ID,
    gene_symbol = toupper(trimws(gpl$`Gene Symbol`)),
    gene_title = gpl$`Gene Title`,
    entrez_id = gpl$ENTREZ_GENE_ID,
    stringsAsFactors = FALSE
  )
  probe_map$present_in_series <- probe_map$probe_id %in% rownames(probe_expression)
  probe_map$reliable_symbol <- !is.na(probe_map$gene_symbol) & probe_map$gene_symbol != "" &
    !grepl("///", probe_map$gene_symbol, fixed = TRUE)
  probe_map$included_before_collapse <- probe_map$present_in_series & probe_map$reliable_symbol
  probe_map$overall_mean_expression <- NA_real_
  means <- rowMeans(probe_expression, na.rm = TRUE)
  probe_map$overall_mean_expression[probe_map$present_in_series] <-
    means[probe_map$probe_id[probe_map$present_in_series]]
  eligible <- probe_map[probe_map$included_before_collapse, , drop = FALSE]
  eligible <- eligible[order(eligible$gene_symbol, -eligible$overall_mean_expression, eligible$probe_id), , drop = FALSE]
  selected <- eligible[!duplicated(eligible$gene_symbol), , drop = FALSE]
  selected$tie_rule <- "highest overall mean; lexical probe_id if exact mean tie"
  gene_expression <- probe_expression[selected$probe_id, , drop = FALSE]
  rownames(gene_expression) <- selected$gene_symbol
  list(probe_audit = probe_map, selected_map = selected, gene_expression = gene_expression)
}
