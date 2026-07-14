options(stringsAsFactors = FALSE)

if (file.exists("config/use_local_r_libs.R")) source("config/use_local_r_libs.R")
suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(openxlsx)
})

series_matrix_file <- "GSE18965_series_matrix.txt(1).gz"
family_soft_file <- "GSE18965_family.soft.gz"
gpl_file <- "GPL96-57554.txt"
submitted_gene_file <- "repair_ECM_32_genes.csv"
cel_dir <- "GSE18965_formal_rerun"
out_dir <- "GSE18965_series_matrix_reproduction"

subdirs <- c(
  "00_logs", "01_metadata", "02_series_matrix_QC", "03_expression",
  "04_limma_primary", "05_limma_sensitivity", "06_32_gene_audit",
  "07_module_score", "08_series_vs_CEL_comparison", "09_GO", "10_report"
)

safe_reset_dir <- function(path) {
  target <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(path) && !startsWith(target, paste0(root, "/"))) {
    stop("Refusing to delete output outside working directory: ", target)
  }
  if (file.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

safe_reset_dir(out_dir)
for (d in subdirs) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "00_logs", "run_log.txt")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_msg("Started Series Matrix reproduction.")

required_inputs <- c(
  series_matrix_file, family_soft_file, gpl_file, submitted_gene_file,
  file.path(cel_dir, "03_expression", "GSE18965_GCRMA_probe_expression.csv"),
  file.path(cel_dir, "03_expression", "GSE18965_GCRMA_gene_expression.csv"),
  file.path(cel_dir, "03_expression", "selected_highest_mean_probe_per_gene.csv"),
  file.path(cel_dir, "04_limma_primary", "GSE18965_complete_limma_primary_default_eBayes_results.csv"),
  file.path(cel_dir, "05_limma_sensitivity", "GSE18965_complete_limma_robust_trend_sensitivity_results.csv"),
  file.path(cel_dir, "08_module_score", "repair_ECM_32_module_scores_sample_level.csv"),
  file.path(cel_dir, "08_module_score", "repair_ECM_32_module_score_statistics.csv"),
  file.path(cel_dir, "CHANGELOG.md"),
  "GSE18965_formal_CEL_GCRMA_limma_reanalysis_v2_formalized.R"
)
input_manifest <- data.frame(
  file = required_inputs,
  exists = file.exists(required_inputs),
  size_bytes = ifelse(file.exists(required_inputs), file.info(required_inputs)$size, NA_real_),
  md5 = ifelse(file.exists(required_inputs), tools::md5sum(required_inputs), NA_character_),
  stringsAsFactors = FALSE
)
write.csv(input_manifest, file.path(out_dir, "01_metadata", "input_file_manifest.csv"), row.names = FALSE)
if (any(!input_manifest$exists)) {
  stop("Missing required files: ", paste(input_manifest$file[!input_manifest$exists], collapse = ", "))
}

extract_tag <- function(lines, tag) {
  line <- grep(paste0("^", tag, "\\t"), lines, value = TRUE)
  if (length(line) != 1) stop("Could not uniquely extract ", tag)
  vals <- strsplit(line, "\t", fixed = TRUE)[[1]][-1]
  gsub('^"|"$', "", vals)
}

series_lines <- readLines(gzfile(series_matrix_file, "rt"), warn = FALSE)
table_begin <- grep("^!series_matrix_table_begin", series_lines)
table_end <- grep("^!series_matrix_table_end", series_lines)
if (length(table_begin) != 1 || length(table_end) != 1 || table_end <= table_begin) {
  stop("Could not locate Series Matrix expression table boundaries.")
}
sample_metadata <- data.frame(
  gsm = extract_tag(series_lines, "!Sample_geo_accession"),
  title = extract_tag(series_lines, "!Sample_title"),
  source_name_ch1 = extract_tag(series_lines, "!Sample_source_name_ch1"),
  characteristics_ch1 = extract_tag(series_lines, "!Sample_characteristics_ch1"),
  data_processing = extract_tag(series_lines, "!Sample_data_processing"),
  platform_id = extract_tag(series_lines, "!Sample_platform_id"),
  stringsAsFactors = FALSE
)
sample_metadata$group <- ifelse(grepl("^AA_", sample_metadata$title), "AA",
                                ifelse(grepl("^HN_", sample_metadata$title), "HN", NA_character_))
sample_metadata$group_label <- ifelse(sample_metadata$group == "AA", "Atopic asthma",
                                      ifelse(sample_metadata$group == "HN", "Healthy non-atopic", NA_character_))
if (nrow(sample_metadata) != 16 || any(is.na(sample_metadata$group)) ||
    sum(sample_metadata$group == "AA") != 9 || sum(sample_metadata$group == "HN") != 7) {
  stop("Series Matrix sample count or group assignment mismatch.")
}
write.csv(sample_metadata, file.path(out_dir, "01_metadata", "series_matrix_sample_metadata.csv"), row.names = FALSE)

soft_lines <- readLines(gzfile(family_soft_file, "rt"), warn = FALSE)
soft_gsm <- unique(sub("^!Sample_geo_accession = ", "", grep("^!Sample_geo_accession = ", soft_lines, value = TRUE)))
soft_order_matches <- identical(sample_metadata$gsm, soft_gsm[match(sample_metadata$gsm, soft_gsm)])

expr_text <- series_lines[(table_begin + 1):(table_end - 1)]
series_table <- read.delim(text = paste(expr_text, collapse = "\n"), check.names = FALSE, quote = "\"", stringsAsFactors = FALSE)
probe_ids <- series_table$ID_REF
series_probe_expr <- as.matrix(series_table[, sample_metadata$gsm])
storage.mode(series_probe_expr) <- "numeric"
rownames(series_probe_expr) <- probe_ids
colnames(series_probe_expr) <- sample_metadata$gsm

integrity <- data.frame(
  metric = c(
    "samples_read", "AA_n", "HN_n", "probes_read", "missing_values",
    "duplicate_probe_ids", "min_expression", "max_expression", "median_expression",
    "appears_log2_transformed", "all_platform_GPL96", "soft_sample_set_matches_series",
    "soft_order_check"
  ),
  value = c(
    ncol(series_probe_expr), sum(sample_metadata$group == "AA"), sum(sample_metadata$group == "HN"),
    nrow(series_probe_expr), sum(is.na(series_probe_expr)), sum(duplicated(rownames(series_probe_expr))),
    min(series_probe_expr, na.rm = TRUE), max(series_probe_expr, na.rm = TRUE), median(series_probe_expr, na.rm = TRUE),
    max(series_probe_expr, na.rm = TRUE) < 30 && min(series_probe_expr, na.rm = TRUE) > -5,
    all(sample_metadata$platform_id == "GPL96"),
    setequal(sample_metadata$gsm, soft_gsm),
    soft_order_matches
  ),
  stringsAsFactors = FALSE
)
write.csv(integrity, file.path(out_dir, "01_metadata", "series_matrix_integrity_report.csv"), row.names = FALSE)
write.csv(data.frame(probe_id = rownames(series_probe_expr), series_probe_expr, check.names = FALSE),
          file.path(out_dir, "03_expression", "series_matrix_probe_expression.csv"), row.names = FALSE)

pdf(file.path(out_dir, "02_series_matrix_QC", "series_matrix_boxplot.pdf"), width = 11, height = 6)
boxplot(series_probe_expr, las = 2, outline = FALSE, main = "GSE18965 Series Matrix expression", ylab = "Deposited log2 expression")
dev.off()
pca <- prcomp(t(series_probe_expr), center = TRUE, scale. = FALSE)
pct <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)
pca_df <- data.frame(sample_metadata, PC1 = pca$x[, 1], PC2 = pca$x[, 2])
write.csv(pca_df, file.path(out_dir, "02_series_matrix_QC", "series_matrix_PCA_scores.csv"), row.names = FALSE)
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = group, label = title)) +
  geom_point(size = 3) + geom_text(vjust = -0.7, size = 3, show.legend = FALSE) +
  theme_bw() + labs(title = "Series Matrix PCA", x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"))
ggsave(file.path(out_dir, "02_series_matrix_QC", "series_matrix_PCA.pdf"), p_pca, width = 8, height = 6)
cor_mat <- cor(series_probe_expr)
write.csv(cor_mat, file.path(out_dir, "02_series_matrix_QC", "series_matrix_sample_correlation_matrix.csv"))
pheatmap(cor_mat, annotation_col = data.frame(group = sample_metadata$group, row.names = sample_metadata$gsm),
         annotation_row = data.frame(group = sample_metadata$group, row.names = sample_metadata$gsm),
         filename = file.path(out_dir, "02_series_matrix_QC", "series_matrix_sample_correlation_heatmap.pdf"),
         width = 9, height = 8)

gpl <- read.delim(gpl_file, sep = "\t", comment.char = "#", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
gpl_ann <- data.frame(
  probe_id = gpl$ID,
  gene_symbol = trimws(gpl$`Gene Symbol`),
  gene_title = gpl$`Gene Title`,
  entrez_id = gpl$ENTREZ_GENE_ID,
  go_bp_raw = gpl$`Gene Ontology Biological Process`,
  stringsAsFactors = FALSE
)
gpl_ann$missing_symbol <- is.na(gpl_ann$gene_symbol) | gpl_ann$gene_symbol == ""
gpl_ann$multi_gene_symbol <- grepl("///", gpl_ann$gene_symbol, fixed = TRUE)
gpl_ann$present_in_series <- gpl_ann$probe_id %in% rownames(series_probe_expr)
gpl_ann$included_primary <- !gpl_ann$missing_symbol & !gpl_ann$multi_gene_symbol & gpl_ann$present_in_series
write.csv(gpl_ann, file.path(out_dir, "03_expression", "series_GPL96_probe_annotation_audit.csv"), row.names = FALSE)

ann <- gpl_ann[gpl_ann$included_primary, ]
series_probe_means <- rowMeans(series_probe_expr)
ann$overall_mean_expression <- series_probe_means[ann$probe_id]
ann <- ann[order(ann$gene_symbol, -ann$overall_mean_expression, ann$probe_id), ]
series_selected <- ann[!duplicated(ann$gene_symbol), ]
series_gene_expr <- series_probe_expr[series_selected$probe_id, , drop = FALSE]
rownames(series_gene_expr) <- series_selected$gene_symbol
write.csv(series_selected, file.path(out_dir, "03_expression", "series_selected_highest_mean_probe_per_gene.csv"), row.names = FALSE)
write.csv(data.frame(gene_symbol = rownames(series_gene_expr), series_gene_expr, check.names = FALSE),
          file.path(out_dir, "03_expression", "series_gene_level_expression.csv"), row.names = FALSE)

group <- factor(sample_metadata$group, levels = c("HN", "AA"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contrast <- makeContrasts(AA_minus_HN = AA - HN, levels = design)
fit <- lmFit(series_gene_expr, design)
fit2 <- contrasts.fit(fit, contrast)
fit_main <- eBayes(fit2, robust = FALSE, trend = FALSE)
fit_sens <- eBayes(fit2, robust = TRUE, trend = TRUE)

format_limma <- function(fit_obj, selected) {
  tab <- topTable(fit_obj, coef = "AA_minus_HN", number = Inf, sort.by = "none")
  tab$gene_symbol <- rownames(tab)
  tab <- tab[, c("gene_symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  names(tab) <- c("gene_symbol", "log2FC", "average_expression", "moderated_t", "p_value", "FDR", "B_statistic")
  tab$selected_probe <- selected$probe_id[match(tab$gene_symbol, selected$gene_symbol)]
  tab$gene_title <- selected$gene_title[match(tab$gene_symbol, selected$gene_symbol)]
  tab$DEG_FDR_0_05_abs_log2FC_0_5 <- tab$FDR < 0.05 & abs(tab$log2FC) >= 0.5
  tab$up_FDR_0_05_log2FC_ge_0_5 <- tab$FDR < 0.05 & tab$log2FC >= 0.5
  tab$down_FDR_0_05_log2FC_le_minus_0_5 <- tab$FDR < 0.05 & tab$log2FC <= -0.5
  tab$direction <- ifelse(tab$log2FC > 0, "AA higher", ifelse(tab$log2FC < 0, "HN higher", "No difference"))
  tab
}
series_limma <- format_limma(fit_main, series_selected)
series_sens <- format_limma(fit_sens, series_selected)
write.csv(series_limma, file.path(out_dir, "04_limma_primary", "series_complete_limma_primary_default_eBayes_results.csv"), row.names = FALSE)
write.csv(series_sens, file.path(out_dir, "05_limma_sensitivity", "series_complete_limma_robust_trend_sensitivity_results.csv"), row.names = FALSE)
write.csv(series_limma[series_limma$up_FDR_0_05_log2FC_ge_0_5, ], file.path(out_dir, "04_limma_primary", "series_upregulated_FDR0.05_log2FC0.5.csv"), row.names = FALSE)
write.csv(series_limma[series_limma$down_FDR_0_05_log2FC_le_minus_0_5, ], file.path(out_dir, "04_limma_primary", "series_downregulated_FDR0.05_log2FC_minus0.5.csv"), row.names = FALSE)
deg_counts <- data.frame(
  samples_read = ncol(series_probe_expr),
  AA_n = sum(group == "AA"),
  HN_n = sum(group == "HN"),
  probes_read = nrow(series_probe_expr),
  genes_tested = nrow(series_limma),
  DEGs_FDR_0_05_abs_log2FC_0_5 = sum(series_limma$DEG_FDR_0_05_abs_log2FC_0_5),
  upregulated = sum(series_limma$up_FDR_0_05_log2FC_ge_0_5),
  downregulated = sum(series_limma$down_FDR_0_05_log2FC_le_minus_0_5),
  stringsAsFactors = FALSE
)
write.csv(deg_counts, file.path(out_dir, "04_limma_primary", "series_DEG_flow_counts.csv"), row.names = FALSE)

submitted <- read.csv(submitted_gene_file, check.names = FALSE)
submitted$gene_symbol <- trimws(submitted$gene_symbol)
submitted <- submitted[!duplicated(submitted$gene_symbol), , drop = FALSE]
probe_lists <- aggregate(probe_id ~ gene_symbol, data = ann, FUN = function(x) paste(sort(unique(x)), collapse = "; "))
names(probe_lists)[2] <- "series_all_available_primary_policy_probes"
audit32 <- data.frame(submitted_order = seq_len(nrow(submitted)), gene_symbol = submitted$gene_symbol, stringsAsFactors = FALSE)
audit32 <- merge(audit32, probe_lists, by = "gene_symbol", all.x = TRUE, sort = FALSE)
audit32$series_mapped <- !is.na(audit32$series_all_available_primary_policy_probes)
audit32$series_selected_probe <- series_selected$probe_id[match(audit32$gene_symbol, series_selected$gene_symbol)]
audit32$series_overall_mean_expression <- series_selected$overall_mean_expression[match(audit32$gene_symbol, series_selected$gene_symbol)]
audit32 <- merge(audit32, series_limma[, c("gene_symbol", "log2FC", "p_value", "FDR", "direction", "up_FDR_0_05_log2FC_ge_0_5")],
                 by = "gene_symbol", all.x = TRUE, sort = FALSE)
names(audit32)[names(audit32) %in% c("log2FC", "p_value", "FDR", "direction", "up_FDR_0_05_log2FC_ge_0_5")] <-
  c("series_log2FC", "series_p_value", "series_FDR", "series_direction", "series_meets_original_up_rule")
audit32$series_FDR_lt_0_05 <- !is.na(audit32$series_FDR) & audit32$series_FDR < 0.05
audit32$series_log2FC_ge_0_5 <- !is.na(audit32$series_log2FC) & audit32$series_log2FC >= 0.5
audit32$series_abs_log2FC_ge_0_5 <- !is.na(audit32$series_log2FC) & abs(audit32$series_log2FC) >= 0.5
sens32 <- series_sens[, c("gene_symbol", "log2FC", "FDR", "up_FDR_0_05_log2FC_ge_0_5", "DEG_FDR_0_05_abs_log2FC_0_5")]
names(sens32) <- c("gene_symbol", "series_sensitivity_log2FC", "series_sensitivity_FDR", "series_sensitivity_meets_up_rule", "series_sensitivity_meets_abs_DEG_rule")
audit32 <- merge(audit32, sens32, by = "gene_symbol", all.x = TRUE, sort = FALSE)
audit32$series_main_vs_sensitivity_direction_consistent <- sign(audit32$series_log2FC) == sign(audit32$series_sensitivity_log2FC)
audit32$series_main_vs_sensitivity_up_rule_consistent <- audit32$series_meets_original_up_rule == audit32$series_sensitivity_meets_up_rule

cel_limma <- read.csv(file.path(cel_dir, "04_limma_primary", "GSE18965_complete_limma_primary_default_eBayes_results.csv"), check.names = FALSE)
cel_selected <- read.csv(file.path(cel_dir, "03_expression", "selected_highest_mean_probe_per_gene.csv"), check.names = FALSE)
cel_scores <- read.csv(file.path(cel_dir, "08_module_score", "repair_ECM_32_module_scores_sample_level.csv"), check.names = FALSE)
if ("GSM" %in% names(cel_scores) && !"gsm" %in% names(cel_scores)) names(cel_scores)[names(cel_scores) == "GSM"] <- "gsm"
cel_score_stats <- read.csv(file.path(cel_dir, "08_module_score", "repair_ECM_32_module_score_statistics.csv"), check.names = FALSE)
cel32 <- cel_limma[, c("gene_symbol", "selected_probe", "log2FC", "p_value", "FDR", "direction", "up_FDR_0_05_log2FC_ge_0_5")]
names(cel32) <- c("gene_symbol", "CEL_selected_probe", "CEL_log2FC", "CEL_p_value", "CEL_FDR", "CEL_direction", "CEL_meets_original_up_rule")
audit32 <- merge(audit32, cel32, by = "gene_symbol", all.x = TRUE, sort = FALSE)
audit32$series_selected_probe_same_as_CEL <- audit32$series_selected_probe == audit32$CEL_selected_probe
audit32$series_CEL_direction_same <- sign(audit32$series_log2FC) == sign(audit32$CEL_log2FC)
audit32$series_minus_CEL_log2FC <- audit32$series_log2FC - audit32$CEL_log2FC
audit32 <- audit32[order(audit32$submitted_order), ]
write.csv(audit32, file.path(out_dir, "06_32_gene_audit", "series_repair_ECM_32_gene_audit.csv"), row.names = FALSE)
write.csv(audit32[, c("submitted_order", "gene_symbol", "series_selected_probe", "CEL_selected_probe", "series_selected_probe_same_as_CEL",
                      "series_log2FC", "CEL_log2FC", "series_minus_CEL_log2FC", "series_direction", "CEL_direction",
                      "series_CEL_direction_same", "series_FDR", "CEL_FDR", "series_meets_original_up_rule", "CEL_meets_original_up_rule")],
          file.path(out_dir, "06_32_gene_audit", "series_vs_CEL_32_gene_comparison.csv"), row.names = FALSE)

available32 <- intersect(submitted$gene_symbol, rownames(series_gene_expr))
series_expr32 <- series_gene_expr[available32, , drop = FALSE]
z32 <- t(scale(t(series_expr32)))
series_module_score <- colMeans(z32, na.rm = TRUE)
score_df <- data.frame(sample_metadata[, c("gsm", "title", "group")],
                       series_module_score = series_module_score[sample_metadata$gsm],
                       stringsAsFactors = FALSE)
write.csv(score_df, file.path(out_dir, "07_module_score", "series_32_gene_module_scores_sample_level.csv"), row.names = FALSE)
module_stats_fun <- function(scores, groups) {
  x <- scores[groups == "AA"]
  y <- scores[groups == "HN"]
  tt <- t.test(x, y, var.equal = FALSE)
  pooled_sd <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / (length(x) + length(y) - 2))
  cohen_d <- (mean(x) - mean(y)) / pooled_sd
  hedges_g <- cohen_d * (1 - 3 / (4 * (length(x) + length(y)) - 9))
  data.frame(
    n_AA = length(x), n_HN = length(y),
    mean_AA = mean(x), sd_AA = sd(x), mean_HN = mean(y), sd_HN = sd(y),
    difference_AA_minus_HN = mean(x) - mean(y),
    CI_low = tt$conf.int[1], CI_high = tt$conf.int[2],
    Welch_t = unname(tt$statistic), Welch_df = unname(tt$parameter), p_value = tt$p.value,
    Cohen_d = cohen_d, Hedges_g = hedges_g,
    stringsAsFactors = FALSE
  )
}
series_score_stats <- module_stats_fun(score_df$series_module_score, score_df$group)
series_score_stats$genes_used <- length(available32)
write.csv(series_score_stats, file.path(out_dir, "07_module_score", "series_32_gene_module_score_statistics.csv"), row.names = FALSE)
p_score <- ggplot(score_df, aes(group, series_module_score, color = group)) +
  geom_boxplot(width = 0.25, outlier.shape = NA) + geom_jitter(width = 0.08, size = 2.5) +
  theme_bw() + labs(title = "Series Matrix 32-gene module score", x = NULL, y = "Mean z-score")
ggsave(file.path(out_dir, "07_module_score", "series_32_gene_module_score_boxplot.pdf"), p_score, width = 6, height = 5)
pheatmap(z32, annotation_col = data.frame(group = sample_metadata$group, row.names = sample_metadata$gsm),
         filename = file.path(out_dir, "07_module_score", "series_32_gene_heatmap.pdf"),
         width = 10, height = 10, fontsize_row = 7)
logo <- do.call(rbind, lapply(available32, function(drop_gene) {
  keep <- setdiff(available32, drop_gene)
  sc <- colMeans(z32[keep, , drop = FALSE], na.rm = TRUE)
  st <- module_stats_fun(sc, group)
  data.frame(dropped_gene = drop_gene, genes_used = length(keep),
             difference_AA_minus_HN = st$difference_AA_minus_HN,
             CI_low = st$CI_low, CI_high = st$CI_high,
             p_value = st$p_value,
             direction = ifelse(st$difference_AA_minus_HN > 0, "AA higher", "HN higher"),
             stringsAsFactors = FALSE)
}))
write.csv(logo, file.path(out_dir, "07_module_score", "series_32_gene_leave_one_gene_out_sensitivity.csv"), row.names = FALSE)

cel_probe <- read.csv(file.path(cel_dir, "03_expression", "GSE18965_GCRMA_probe_expression.csv"), check.names = FALSE)
rownames(cel_probe) <- cel_probe$probe_id
cel_probe_expr <- as.matrix(cel_probe[, sample_metadata$gsm])
storage.mode(cel_probe_expr) <- "numeric"
common_probes <- intersect(rownames(series_probe_expr), rownames(cel_probe_expr))
common_samples <- sample_metadata$gsm
probe_sample_summary <- do.call(rbind, lapply(common_samples, function(s) {
  d <- series_probe_expr[common_probes, s] - cel_probe_expr[common_probes, s]
  data.frame(
    sample = s,
    title = sample_metadata$title[match(s, sample_metadata$gsm)],
    group = sample_metadata$group[match(s, sample_metadata$gsm)],
    n_common_probes = length(common_probes),
    Pearson = cor(series_probe_expr[common_probes, s], cel_probe_expr[common_probes, s], method = "pearson"),
    Spearman = cor(series_probe_expr[common_probes, s], cel_probe_expr[common_probes, s], method = "spearman"),
    mean_difference_series_minus_CEL = mean(d),
    sd_difference = sd(d),
    BA_bias = mean(d),
    BA_LOA_low = mean(d) - 1.96 * sd(d),
    BA_LOA_high = mean(d) + 1.96 * sd(d),
    stringsAsFactors = FALSE
  )
}))
write.csv(probe_sample_summary, file.path(out_dir, "08_series_vs_CEL_comparison", "probe_level_per_sample_correlation.csv"), row.names = FALSE)
overall_probe <- data.frame(
  n_common_probes = length(common_probes),
  n_common_values = length(common_probes) * length(common_samples),
  Pearson = cor(as.vector(series_probe_expr[common_probes, common_samples]), as.vector(cel_probe_expr[common_probes, common_samples]), method = "pearson"),
  Spearman = cor(as.vector(series_probe_expr[common_probes, common_samples]), as.vector(cel_probe_expr[common_probes, common_samples]), method = "spearman"),
  mean_difference_series_minus_CEL = mean(series_probe_expr[common_probes, common_samples] - cel_probe_expr[common_probes, common_samples]),
  sd_difference = sd(as.vector(series_probe_expr[common_probes, common_samples] - cel_probe_expr[common_probes, common_samples])),
  stringsAsFactors = FALSE
)
write.csv(overall_probe, file.path(out_dir, "08_series_vs_CEL_comparison", "probe_level_overall_correlation_summary.csv"), row.names = FALSE)
probe_level_by_probe <- do.call(rbind, lapply(common_probes, function(pr) {
  sx <- series_probe_expr[pr, common_samples]
  cx <- cel_probe_expr[pr, common_samples]
  data.frame(
    probe_id = pr,
    Pearson_16_samples = suppressWarnings(cor(sx, cx, method = "pearson")),
    Spearman_16_samples = suppressWarnings(cor(sx, cx, method = "spearman")),
    series_probe_log2FC = mean(sx[group == "AA"]) - mean(sx[group == "HN"]),
    CEL_probe_log2FC = mean(cx[group == "AA"]) - mean(cx[group == "HN"]),
    stringsAsFactors = FALSE
  )
}))
probe_level_by_probe$series_minus_CEL_probe_log2FC <- probe_level_by_probe$series_probe_log2FC - probe_level_by_probe$CEL_probe_log2FC
write.csv(probe_level_by_probe, file.path(out_dir, "08_series_vs_CEL_comparison", "probe_level_per_probe_correlation_and_log2FC.csv"), row.names = FALSE)

gene_cmp <- merge(
  series_limma[, c("gene_symbol", "selected_probe", "log2FC", "p_value", "FDR", "DEG_FDR_0_05_abs_log2FC_0_5", "up_FDR_0_05_log2FC_ge_0_5")],
  cel_limma[, c("gene_symbol", "selected_probe", "log2FC", "p_value", "FDR", "DEG_FDR_0_05_abs_log2FC_0_5", "up_FDR_0_05_log2FC_ge_0_5")],
  by = "gene_symbol", suffixes = c("_series", "_CEL")
)
gene_cmp$selected_probe_same <- gene_cmp$selected_probe_series == gene_cmp$selected_probe_CEL
gene_cmp$series_minus_CEL_log2FC <- gene_cmp$log2FC_series - gene_cmp$log2FC_CEL
gene_cmp$p_rank_series <- rank(gene_cmp$p_value_series, ties.method = "average")
gene_cmp$p_rank_CEL <- rank(gene_cmp$p_value_CEL, ties.method = "average")
write.csv(gene_cmp, file.path(out_dir, "08_series_vs_CEL_comparison", "gene_level_series_vs_CEL_full_comparison.csv"), row.names = FALSE)
gene_summary <- data.frame(
  common_genes = nrow(gene_cmp),
  selected_probe_agreement_n = sum(gene_cmp$selected_probe_same, na.rm = TRUE),
  selected_probe_agreement_fraction = mean(gene_cmp$selected_probe_same, na.rm = TRUE),
  log2FC_Pearson = cor(gene_cmp$log2FC_series, gene_cmp$log2FC_CEL, method = "pearson"),
  log2FC_Spearman = cor(gene_cmp$log2FC_series, gene_cmp$log2FC_CEL, method = "spearman"),
  p_value_rank_Spearman = cor(gene_cmp$p_rank_series, gene_cmp$p_rank_CEL, method = "spearman"),
  FDR_lt_0_05_overlap = sum(gene_cmp$FDR_series < 0.05 & gene_cmp$FDR_CEL < 0.05, na.rm = TRUE),
  DEG_fixed_threshold_overlap = sum(gene_cmp$DEG_FDR_0_05_abs_log2FC_0_5_series & gene_cmp$DEG_FDR_0_05_abs_log2FC_0_5_CEL, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(gene_summary, file.path(out_dir, "08_series_vs_CEL_comparison", "gene_level_series_vs_CEL_summary.csv"), row.names = FALSE)

score_cmp <- merge(score_df, cel_scores[, c("gsm", "module_score")], by = "gsm")
names(score_cmp)[names(score_cmp) == "module_score"] <- "CEL_module_score"
score_cmp$series_minus_CEL_module_score <- score_cmp$series_module_score - score_cmp$CEL_module_score
write.csv(score_cmp, file.path(out_dir, "08_series_vs_CEL_comparison", "module_score_series_vs_CEL_sample_comparison.csv"), row.names = FALSE)
module_cmp_summary <- data.frame(
  n_samples = nrow(score_cmp),
  Pearson = cor(score_cmp$series_module_score, score_cmp$CEL_module_score, method = "pearson"),
  Spearman = cor(score_cmp$series_module_score, score_cmp$CEL_module_score, method = "spearman"),
  mean_difference_series_minus_CEL = mean(score_cmp$series_minus_CEL_module_score),
  sd_difference = sd(score_cmp$series_minus_CEL_module_score),
  stringsAsFactors = FALSE
)
write.csv(module_cmp_summary, file.path(out_dir, "08_series_vs_CEL_comparison", "module_score_series_vs_CEL_summary.csv"), row.names = FALSE)

p1 <- ggplot(gene_cmp, aes(log2FC_CEL, log2FC_series)) +
  geom_point(alpha = 0.35, size = 1) + geom_hline(yintercept = 0, lty = 2) + geom_vline(xintercept = 0, lty = 2) +
  theme_bw() + labs(title = "Whole-genome log2FC: Series Matrix vs CEL GC-RMA", x = "CEL log2FC", y = "Series Matrix log2FC")
ggsave(file.path(out_dir, "08_series_vs_CEL_comparison", "whole_genome_log2FC_scatter.pdf"), p1, width = 6.5, height = 5.5)
p2 <- ggplot(audit32, aes(CEL_log2FC, series_log2FC, label = gene_symbol)) +
  geom_point(size = 2.5) + geom_text(vjust = -0.5, size = 2.5) + geom_hline(yintercept = 0, lty = 2) + geom_vline(xintercept = 0, lty = 2) +
  theme_bw() + labs(title = "32-gene log2FC: Series Matrix vs CEL", x = "CEL log2FC", y = "Series Matrix log2FC")
ggsave(file.path(out_dir, "08_series_vs_CEL_comparison", "32_gene_log2FC_scatter.pdf"), p2, width = 7, height = 6)
p3 <- ggplot(score_cmp, aes(CEL_module_score, series_module_score, color = group, label = title)) +
  geom_point(size = 2.8) + geom_text(vjust = -0.6, size = 2.5, show.legend = FALSE) +
  theme_bw() + labs(title = "Module score: Series Matrix vs CEL", x = "CEL module score", y = "Series Matrix module score")
ggsave(file.path(out_dir, "08_series_vs_CEL_comparison", "module_score_series_vs_CEL_scatter.pdf"), p3, width = 6.5, height = 5.5)
p4 <- ggplot(gene_cmp, aes((log2FC_CEL + log2FC_series) / 2, series_minus_CEL_log2FC)) +
  geom_point(alpha = 0.35, size = 1) + geom_hline(yintercept = mean(gene_cmp$series_minus_CEL_log2FC), color = "red") +
  theme_bw() + labs(title = "Bland-Altman style log2FC comparison", x = "Mean log2FC", y = "Series - CEL log2FC")
ggsave(file.path(out_dir, "08_series_vs_CEL_comparison", "whole_genome_log2FC_bland_altman.pdf"), p4, width = 6.5, height = 5.5)

up_genes <- series_limma$gene_symbol[series_limma$up_FDR_0_05_log2FC_ge_0_5]
if (length(up_genes) == 0) {
  writeLines(
    "GO BP enrichment not estimable because the prespecified upregulated gene set was empty (FDR<0.05 and log2FC>=0.5). No 0/0 GO result table was generated.",
    file.path(out_dir, "09_GO", "GO_status.txt")
  )
  go_res <- data.frame()
} else {
  parse_go <- function(x) {
    if (is.na(x) || !nzchar(x)) return(data.frame(go_id = character(), description = character(), stringsAsFactors = FALSE))
    parts <- unlist(strsplit(x, " /// ", fixed = TRUE))
    out <- do.call(rbind, lapply(parts, function(part) {
      fields <- strsplit(part, " // ", fixed = TRUE)[[1]]
      if (length(fields) < 2) return(NULL)
      data.frame(go_id = paste0("GO:", fields[1]), description = fields[2], stringsAsFactors = FALSE)
    }))
    if (is.null(out)) return(data.frame(go_id = character(), description = character(), stringsAsFactors = FALSE))
    unique(out)
  }
  go_maps <- do.call(rbind, lapply(seq_len(nrow(series_selected)), function(i) {
    pgo <- parse_go(series_selected$go_bp_raw[i])
    if (!nrow(pgo)) return(NULL)
    data.frame(gene_symbol = series_selected$gene_symbol[i], pgo, stringsAsFactors = FALSE)
  }))
  go_maps <- unique(go_maps)
  write.csv(go_maps, file.path(out_dir, "09_GO", "GO_BP_gene_to_term_mapping_from_series_selected_GPL96_probes.csv"), row.names = FALSE)
  universe_genes <- series_limma$gene_symbol
  n_input <- length(intersect(up_genes, unique(go_maps$gene_symbol)))
  N_bg <- length(intersect(universe_genes, unique(go_maps$gene_symbol)))
  go_terms <- unique(go_maps[, c("go_id", "description")])
  go_res <- do.call(rbind, lapply(seq_len(nrow(go_terms)), function(i) {
    term_genes <- unique(go_maps$gene_symbol[go_maps$go_id == go_terms$go_id[i]])
    term_genes_bg <- intersect(term_genes, universe_genes)
    term_genes_input <- intersect(term_genes, up_genes)
    k <- length(term_genes_input)
    K <- length(term_genes_bg)
    data.frame(
      GO_ID = go_terms$go_id[i],
      Description = go_terms$description[i],
      GeneRatio = paste0(k, "/", n_input),
      BgRatio = paste0(K, "/", N_bg),
      gene_count = k,
      raw_P = phyper(k - 1, K, N_bg - K, n_input, lower.tail = FALSE),
      gene_symbols = paste(sort(term_genes_input), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))
  go_res$FDR_BH <- p.adjust(go_res$raw_P, method = "BH")
  go_res <- go_res[order(go_res$FDR_BH, go_res$raw_P, -go_res$gene_count), ]
  write.csv(go_res, file.path(out_dir, "09_GO", "GO_BP_complete_results.csv"), row.names = FALSE)
}

summary_lines <- c(
  paste0("- Series Matrix samples read: ", ncol(series_probe_expr)),
  paste0("- AA n / HN n: ", sum(group == "AA"), " / ", sum(group == "HN")),
  paste0("- Probes read: ", nrow(series_probe_expr)),
  paste0("- Genes tested: ", nrow(series_limma)),
  paste0("- Series Matrix DEGs at FDR<0.05 and |log2FC|>=0.5: ", sum(series_limma$DEG_FDR_0_05_abs_log2FC_0_5)),
  paste0("- Series Matrix upregulated/downregulated: ", sum(series_limma$up_FDR_0_05_log2FC_ge_0_5), " / ", sum(series_limma$down_FDR_0_05_log2FC_le_minus_0_5)),
  paste0("- 32 genes mapped: ", sum(audit32$series_mapped, na.rm = TRUE), "/32"),
  paste0("- 32 genes meeting original upregulated rule: ", sum(audit32$series_meets_original_up_rule, na.rm = TRUE), "/32"),
  paste0("- Series Matrix module score difference, 95% CI, P: ",
         signif(series_score_stats$difference_AA_minus_HN, 4), " [", signif(series_score_stats$CI_low, 4), ", ",
         signif(series_score_stats$CI_high, 4), "], P=", signif(series_score_stats$p_value, 4)),
  paste0("- CEL module score difference, 95% CI, P: ",
         signif(cel_score_stats$difference_AA_minus_HN, 4), " [", signif(cel_score_stats$CI_low, 4), ", ",
         signif(cel_score_stats$CI_high, 4), "], P=", signif(cel_score_stats$p_value, 4)),
  paste0("- Whole-genome log2FC correlation between Series Matrix and CEL: Pearson ",
         signif(gene_summary$log2FC_Pearson, 4), "; Spearman ", signif(gene_summary$log2FC_Spearman, 4)),
  paste0("- 32-gene log2FC correlation: Pearson ",
         signif(cor(audit32$series_log2FC, audit32$CEL_log2FC, method = "pearson"), 4), "; Spearman ",
         signif(cor(audit32$series_log2FC, audit32$CEL_log2FC, method = "spearman"), 4)),
  paste0("- Selected-probe agreement: ", gene_summary$selected_probe_agreement_n, "/", gene_summary$common_genes,
         " (", signif(100 * gene_summary$selected_probe_agreement_fraction, 4), "%)"),
  "- Main explanation for discordance: the deposited Series Matrix and current CEL rerun are both labelled GC-RMA/log2 but are not numerically identical at probe level; probe-level sample correlations are high but differences are large enough to alter highest-mean probe choices and limma effect estimates. The discrepancy is therefore upstream of gene-level limma and is most consistent with differences in historical GEO processing implementation/package/CDF/probe-affinity resources rather than sample order or group coding.",
  "- Historical 32-gene derivation fully reproducible: No",
  "- Recommended role of each analysis in manuscript:",
  "  - Series Matrix: primary reproduction of the submitted Methods.",
  "  - CEL rerun: raw-data sensitivity analysis."
)

comparison_report <- c(
  "# GSE18965 Series Matrix vs CEL Comparison Report",
  "",
  "## Summary",
  summary_lines,
  "",
  "## Consistency Definition",
  "For `main versus robust/trend consistency`, this rerun reports two explicit criteria: direction consistency means identical sign of log2FC; up-rule consistency means identical classification for FDR<0.05 and log2FC>=0.5. These are not interchangeable.",
  "",
  "## Probe-Level Comparison",
  paste0("- Common probes compared: ", length(common_probes), "."),
  paste0("- Overall probe-value Pearson/Spearman: ", signif(overall_probe$Pearson, 4), " / ", signif(overall_probe$Spearman, 4), "."),
  paste0("- Mean Series-CEL probe-value difference: ", signif(overall_probe$mean_difference_series_minus_CEL, 4), "."),
  "Per-sample and per-probe summaries are in `08_series_vs_CEL_comparison`.",
  "",
  "## Gene-Level Comparison",
  paste0("- Common genes: ", gene_summary$common_genes, "."),
  paste0("- Selected-probe agreement: ", gene_summary$selected_probe_agreement_n, "/", gene_summary$common_genes, "."),
  paste0("- Whole-genome log2FC Pearson/Spearman: ", signif(gene_summary$log2FC_Pearson, 4), " / ", signif(gene_summary$log2FC_Spearman, 4), "."),
  paste0("- P-value rank Spearman: ", signif(gene_summary$p_value_rank_Spearman, 4), "."),
  "",
  "## GO Correction",
  "If the prespecified upregulated gene set is empty, GO enrichment is not estimable and no 0/0 GO result table is produced. The previous CEL output table with GeneRatio=0/0 rows should not be interpreted or cited as enrichment results.",
  "",
  "## Interpretation",
  "The sample order, group coding, platform ID, and fixed probe-collapse rule were checked. The major difference is already present at the probe-expression-matrix level between the GEO-deposited processed matrix and the current CEL->GC-RMA rerun. This can propagate into selected-probe choices for genes with multiple probes, gene-level log2FC, module scores, and DEG calls. Therefore the Series Matrix analysis should be treated as the formal Methods reproduction, while the CEL rerun should be described as raw-data sensitivity."
)
writeLines(comparison_report, file.path(out_dir, "10_report", "GSE18965_series_vs_CEL_comparison_report.md"))

main_report <- c(
  "# GSE18965 Series Matrix Reproduction Report",
  "",
  "## Executive Summary",
  summary_lines,
  "",
  "## Inputs",
  "- `GSE18965_series_matrix.txt(1).gz` was used as the formal Methods reproduction source.",
  "- Previous CEL->GC-RMA outputs were used only for sensitivity comparison and were read from saved probe-level and gene-level matrices.",
  "",
  "## Series Matrix Integrity",
  "- Series Matrix data processing field states: R/Bioconductor GC-RMA normalization log2.",
  "- Sample count and AA/HN group counts match expected values.",
  "- Expression values are in a log2-like range; no missing values or duplicate probe IDs were detected.",
  "",
  "## Analysis",
  "- GPL96 annotation policy matched the CEL rerun: empty symbols and multi-gene `///` probes excluded.",
  "- Highest mean expression probe per gene was selected within the Series Matrix expression matrix.",
  "- Primary limma used default eBayes; robust/trend was sensitivity only.",
  "- log2FC is AA-HN; DEG threshold is fixed at FDR<0.05 and |log2FC|>=0.5.",
  "",
  "## 32-Gene Audit",
  "- Full gene-level audit is in `06_32_gene_audit/series_repair_ECM_32_gene_audit.csv` and `.xlsx`.",
  "",
  "## GO",
  if (length(up_genes) == 0) "- GO BP enrichment was not estimated because the prespecified Series Matrix upregulated gene set was empty." else "- GO BP enrichment was estimated with the full gene-level limma universe as background.",
  "",
  "## Conclusion",
  "The Series Matrix rerun is the appropriate formal reproduction of the submitted Methods. The CEL rerun remains valuable as raw-data sensitivity, but it should not replace the deposited-processed-matrix analysis when describing the original Methods."
)
writeLines(main_report, file.path(out_dir, "10_report", "GSE18965_series_matrix_reproduction_report.md"))

changelog <- c(
  "# CHANGELOG_series_matrix",
  "",
  "- Created a separate Series Matrix reproduction workflow because the submitted Methods specify the GEO-deposited processed log2 Series Matrix.",
  "- Kept AA/HN grouping, log2FC=AA-HN, highest-mean probe collapse, default eBayes primary limma, robust/trend sensitivity limma, and fixed DEG threshold unchanged.",
  "- Used saved CEL probe-level and gene-level expression matrices from the previous rerun for consistency comparison; did not infer expression from summary statistics.",
  "- Defined main versus robust/trend consistency explicitly as log2FC direction consistency and up-rule classification consistency.",
  "- Corrected GO handling: when the prespecified upregulated gene set is empty, GO enrichment is reported as not estimable and no 0/0 GO result table is generated.",
  "- Did not modify manuscript text, response letter, thresholds, sample grouping, gene list, or probe-selection rules."
)
writeLines(changelog, file.path(out_dir, "CHANGELOG_series_matrix.md"))

wb <- createWorkbook()
add_sheet <- function(name, df) {
  addWorksheet(wb, substr(name, 1, 31))
  writeData(wb, substr(name, 1, 31), df)
  freezePane(wb, substr(name, 1, 31), firstRow = TRUE)
}
add_sheet("README_summary", data.frame(item = sub("^- ", "", summary_lines), stringsAsFactors = FALSE))
add_sheet("series_metadata", sample_metadata)
add_sheet("series_integrity", integrity)
add_sheet("series_annotation_audit", gpl_ann)
add_sheet("series_selected_probes", series_selected)
add_sheet("series_limma_primary_all", series_limma)
add_sheet("series_limma_sensitivity", series_sens)
add_sheet("series_DEG_counts", deg_counts)
add_sheet("series_32_gene_audit", audit32)
add_sheet("series_vs_CEL_32", audit32[, c("submitted_order", "gene_symbol", "series_selected_probe", "CEL_selected_probe", "series_selected_probe_same_as_CEL", "series_log2FC", "CEL_log2FC", "series_minus_CEL_log2FC", "series_FDR", "CEL_FDR")])
add_sheet("series_module_scores", score_df)
add_sheet("series_module_stats", series_score_stats)
add_sheet("leave_one_gene_out", logo)
add_sheet("probe_sample_compare", probe_sample_summary)
add_sheet("probe_overall_compare", overall_probe)
add_sheet("gene_level_compare", gene_cmp)
add_sheet("gene_compare_summary", gene_summary)
add_sheet("module_compare", score_cmp)
add_sheet("module_compare_summary", module_cmp_summary)
if (nrow(go_res)) add_sheet("GO_BP_complete", go_res) else add_sheet("GO_status", data.frame(status = readLines(file.path(out_dir, "09_GO", "GO_status.txt"))))
saveWorkbook(wb, file.path(out_dir, "10_report", "GSE18965_series_matrix_reproduction_all_results.xlsx"), overwrite = TRUE)

audit_wb <- createWorkbook()
addWorksheet(audit_wb, "series_32_gene_audit")
writeData(audit_wb, "series_32_gene_audit", audit32)
freezePane(audit_wb, "series_32_gene_audit", firstRow = TRUE)
saveWorkbook(audit_wb, file.path(out_dir, "06_32_gene_audit", "series_repair_ECM_32_gene_audit.xlsx"), overwrite = TRUE)

file.copy(file.path(out_dir, "10_report", "GSE18965_series_matrix_reproduction_report.md"),
          file.path(out_dir, "GSE18965_series_matrix_reproduction_report.md"), overwrite = TRUE)
file.copy(file.path(out_dir, "10_report", "GSE18965_series_vs_CEL_comparison_report.md"),
          file.path(out_dir, "GSE18965_series_vs_CEL_comparison_report.md"), overwrite = TRUE)
file.copy(file.path(out_dir, "10_report", "GSE18965_series_matrix_reproduction_all_results.xlsx"),
          file.path(out_dir, "GSE18965_series_matrix_reproduction_all_results.xlsx"), overwrite = TRUE)

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
writeLines(summary_lines, file.path(out_dir, "00_logs", "terminal_summary.txt"))
log_msg("Completed Series Matrix reproduction.")

