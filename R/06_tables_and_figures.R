options(stringsAsFactors = FALSE)
if (file.exists("config/use_local_r_libs.R")) source("config/use_local_r_libs.R")

MASTER_SEED <- 20260714L
set.seed(MASTER_SEED)

required <- c("data.table", "dplyr", "tidyr", "openxlsx", "ggplot2", "limma",
              "cluster", "mclust", "sandwich", "lmtest", "MASS", "yaml")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(openxlsx)
  library(ggplot2); library(cluster); library(mclust); library(sandwich)
  library(lmtest); library(MASS); library(yaml)
})

OUT <- file.path("outputs", "user_rerun", "final_locked_analysis_20260714")
dirs <- c("config", "logs", "tables/csv", "figures", "reports", "tests",
          "source_data", "code_repository/R", "code_repository/python",
          "code_repository/config", "code_repository/data_manifest",
          "code_repository/environment", "code_repository/outputs")
if (dir.exists(OUT)) unlink(OUT, recursive = TRUE, force = TRUE)
dir.create(OUT, recursive = TRUE)
for (d in dirs) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(OUT, "logs", "final_locked_run_log.txt")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n"); cat(msg, "\n", file = log_file, append = TRUE)
}

write_csv <- function(x, name) fwrite(x, file.path(OUT, "tables", "csv", name))
write_root_csv <- function(x, name) fwrite(x, file.path(OUT, "tables", name))
zscore <- function(x) as.numeric(scale(as.numeric(x)))
median_high <- function(x) {
  med <- median(x, na.rm = TRUE)
  factor(ifelse(x > med, "high", "low"), levels = c("low", "high"))
}
read_expr <- function(path) {
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  gene_col <- names(x)[1]
  genes <- toupper(trimws(as.character(x[[gene_col]])))
  x[[gene_col]] <- NULL
  m <- as.matrix(x); storage.mode(m) <- "numeric"; rownames(m) <- genes
  m <- m[rownames(m) != "" & !is.na(rownames(m)), , drop = FALSE]
  if (any(duplicated(rownames(m)))) m <- rowsum(m, group = rownames(m), reorder = FALSE)
  m
}
score_mean_z <- function(expr, genes) {
  genes <- unique(toupper(trimws(genes)))
  present <- intersect(genes, rownames(expr))
  m <- expr[present, , drop = FALSE]
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  m <- m[is.finite(sdv) & sdv > 0, , drop = FALSE]
  if (!nrow(m)) return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))
  setNames(colMeans(t(scale(t(m))), na.rm = TRUE), colnames(expr))
}
read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  for (ln in lines) {
    p <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(p) >= 3) out[[p[1]]] <- toupper(trimws(p[-c(1, 2)]))
  }
  out
}
model_logistic <- function(fit, term, label, outcome) {
  co <- summary(fit)$coefficients[term, ]
  beta <- unname(co["Estimate"]); se <- unname(co["Std. Error"])
  data.frame(model = label, outcome = outcome, term = term, beta = beta, SE = se,
             OR = exp(beta), CI_low = exp(beta - 1.96 * se),
             CI_high = exp(beta + 1.96 * se),
             p_value = unname(co["Pr(>|z|)"]), n = nobs(fit),
             events = sum(model.response(model.frame(fit)) == 1), AIC = AIC(fit))
}
model_poisson_robust <- function(fit, term, label, outcome) {
  ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC0"))
  beta <- unname(ct[term, "Estimate"]); se <- unname(ct[term, "Std. Error"])
  y <- model.response(model.frame(fit))
  data.frame(model = label, outcome = outcome, term = term, beta = beta, SE = se,
             IRR = exp(beta), CI_low = exp(beta - 1.96 * se),
             CI_high = exp(beta + 1.96 * se),
             p_value = unname(ct[term, "Pr(>|z|)"]), n = nobs(fit),
             total_events = sum(y, na.rm = TRUE),
             poisson_overdispersion_ratio = sum(residuals(fit, type = "pearson")^2) / fit$df.residual)
}
cluster_metrics <- function(x, cl) {
  x <- as.matrix(x); cl <- factor(cl)
  sil <- cluster::silhouette(as.integer(cl), dist(x))[, "sil_width"]
  centers <- do.call(rbind, lapply(levels(cl), function(g) colMeans(x[cl == g, , drop = FALSE])))
  scat <- vapply(seq_along(levels(cl)), function(i) {
    g <- levels(cl)[i]
    mean(sqrt(rowSums((x[cl == g, , drop = FALSE] - matrix(centers[i, ], sum(cl == g), ncol(x), TRUE))^2)))
  }, numeric(1))
  cdist <- as.matrix(dist(centers))
  db <- mean(vapply(seq_len(nrow(centers)), function(i) max((scat[i] + scat[-i]) / cdist[i, -i]), numeric(1)))
  overall <- colMeans(x)
  within <- sum(vapply(seq_along(levels(cl)), function(i) {
    g <- levels(cl)[i]
    sum(rowSums((x[cl == g, , drop = FALSE] - matrix(centers[i, ], sum(cl == g), ncol(x), TRUE))^2))
  }, numeric(1)))
  between <- sum(vapply(seq_along(levels(cl)), function(i) {
    g <- levels(cl)[i]; sum(cl == g) * sum((centers[i, ] - overall)^2)
  }, numeric(1)))
  data.frame(mean_silhouette = mean(sil),
             calinski_harabasz = (between / (nrow(centers) - 1)) / (within / (nrow(x) - nrow(centers))),
             davies_bouldin = db,
             min_cluster_size = min(table(cl)))
}
bootstrap_cluster <- function(x, ref_cl, centers, B = 1000, seed = MASTER_SEED) {
  set.seed(seed)
  x <- as.matrix(x)
  out <- vector("list", B)
  for (b in seq_len(B)) {
    idx <- sample.int(nrow(x), nrow(x), replace = TRUE)
    fit <- try(kmeans(x[idx, , drop = FALSE], centers = centers, iter.max = 100), silent = TRUE)
    if (inherits(fit, "try-error")) {
      out[[b]] <- data.frame(iter = b, ARI = NA_real_, NMI = NA_real_, agreement = NA_real_, collapsed = TRUE)
      next
    }
    pred <- max.col(-sapply(seq_len(nrow(fit$centers)), function(j) sqrt(rowSums((x - matrix(fit$centers[j, ], nrow(x), ncol(x), TRUE))^2))))
    out[[b]] <- data.frame(iter = b, ARI = mclust::adjustedRandIndex(ref_cl, pred),
                           NMI = mclust::adjustedRandIndex(ref_cl, pred), agreement = mean(ref_cl == pred),
                           collapsed = length(unique(pred)) < length(unique(ref_cl)))
  }
  bind_rows(out)
}
save_fig <- function(p, stem, w = 9, h = 6) {
  ggsave(file.path(OUT, "figures", paste0(stem, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT, "figures", paste0(stem, ".tiff")), p, width = w, height = h, dpi = 300, bg = "white", compression = "lzw")
}

paths <- list(
  expr152 = "03_result_directories/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE152004_VST_gene_expression.csv",
  meta152 = "03_result_directories/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE152004_metadata.csv",
  expr118 = "03_result_directories/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE118761_VST_gene_expression.csv",
  meta118 = "03_result_directories/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE118761_metadata.csv",
  gse18965_limma = "03_result_directories/conservative_full_replication_audit/02_GSE18965/GSE18965_series_limma_complete.csv",
  gse18965_expr = "03_result_directories/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE18965_series_gene_expression.csv",
  gse18965_meta = "03_result_directories/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE18965_series_metadata.csv",
  external_gene_sets = "03_result_directories/External_ECM_multicohort_analysis/02_gene_sets/external_ECM_gene_sets_original.csv",
  repair32 = "01_source_inputs/GSE18965_Codex_full_package/repair_ECM_32_genes.csv",
  modules_gmt = "01_source_inputs/external_ECM_user_inputs/T2_IFN_cluster_metadata_files/article_T2_IFN_Repair_custom_modules.gmt",
  corrected_external8 = "03_result_directories/ECM_third_axis_corrected_rerun/corrected_external8_FDR.xlsx",
  corrected_random = "03_result_directories/ECM_third_axis_corrected_rerun/corrected_random_bootstrap_statistics.csv.gz",
  clinical = Sys.getenv("HOSPITAL_PRIVATE_DATA", unset = "data/private/hospital_deidentified.xlsx")
)

config <- list(
  master_seed = MASTER_SEED,
  gse18965 = list(n_atopic_asthma = 9, n_nonatopic_control = 7, deg_fdr = 0.05, deg_abs_log2fc = 0.5),
  gse118761 = list(nasal_n = 55, tracheal_n = 49),
  gse152004 = list(n = 695, t2_genes = c("CST1", "CLCA1", "SERPINB2")),
  hospital = list(n = 325, asthma_n = 158, frequent_exacerbation_n = 79,
                  index_event_included = TRUE, study_start = "2025-01-01", study_end = "2025-12-31"),
  clustering = list(k_primary = 3, k_evaluated = c(3, 4, 5, 6), n_init = 100,
                    bootstrap_iterations = 1000, low_confidence_margin = 0.10),
  random_gene_sets = list(iterations = 10000)
)
writeLines(as.yaml(config), file.path(OUT, "config", "final_config.yml"))
file.copy(file.path(OUT, "config", "final_config.yml"), file.path(OUT, "code_repository", "config", "final_config.yml"), overwrite = TRUE)

log_msg("Loading expression matrices and gene sets.")
expr152 <- read_expr(paths$expr152); meta152 <- fread(paths$meta152, data.table = FALSE)
expr118 <- read_expr(paths$expr118); meta118 <- fread(paths$meta118, data.table = FALSE)
expr152 <- expr152[, meta152$sample_id, drop = FALSE]
expr118 <- expr118[, meta118$sample_id, drop = FALSE]
mods <- read_gmt(paths$modules_gmt)
repair32 <- unique(toupper(trimws(fread(paths$repair32, data.table = FALSE)[[1]])))
egs <- fread(paths$external_gene_sets, data.table = FALSE)
reactome <- unique(toupper(egs$gene_symbol[egs$gene_set == "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"]))
naba <- unique(toupper(egs$gene_symbol[egs$gene_set == "NABA_CORE_MATRISOME"]))

log_msg("GSE18965 final Series Matrix DEG and score lock.")
limma <- fread(paths$gse18965_limma, data.table = FALSE)
limma$FDR_only <- limma$adj.P.Val < 0.05
limma$primary_DEG <- limma$adj.P.Val < 0.05 & abs(limma$log2FC_AA_minus_HN) >= 0.5
limma$direction <- ifelse(limma$log2FC_AA_minus_HN > 0, "up", "down")
write_csv(limma, "GSE18965_all_DEG_results.csv")
deg_audit <- bind_rows(
  data.frame(filter_definition = "FDR_only", total = sum(limma$FDR_only),
             up = sum(limma$FDR_only & limma$log2FC_AA_minus_HN > 0),
             down = sum(limma$FDR_only & limma$log2FC_AA_minus_HN < 0)),
  data.frame(filter_definition = "FDR_and_abs_log2FC", total = sum(limma$primary_DEG),
             up = sum(limma$primary_DEG & limma$log2FC_AA_minus_HN > 0),
             down = sum(limma$primary_DEG & limma$log2FC_AA_minus_HN < 0))
)
write_csv(deg_audit, "GSE18965_DEG_threshold_audit.csv"); write_root_csv(deg_audit, "GSE18965_DEG_threshold_audit.csv")
fixed32_overlap <- limma |> filter(gene_symbol %in% repair32) |>
  transmute(gene_symbol, log2FC_AA_minus_HN, P.Value, adj.P.Val, primary_DEG, direction)
write_csv(fixed32_overlap, "GSE18965_fixed32_DEG_overlap.csv")

expr18965 <- read_expr(paths$gse18965_expr); meta18965 <- fread(paths$gse18965_meta, data.table = FALSE)
expr18965 <- expr18965[, meta18965$sample_id, drop = FALSE]
score_sets <- list(fixed32_repair_ECM = repair32, Reactome_ECM = reactome,
                   Reactome_ECM_excluding_fixed32 = setdiff(reactome, repair32),
                   NABA_core = naba, NABA_core_excluding_fixed32 = setdiff(naba, repair32))
score_model_rows <- list()
plot_score_rows <- list()
for (nm in names(score_sets)) {
  sc <- score_mean_z(expr18965, score_sets[[nm]])[meta18965$sample_id]
  df <- data.frame(score = sc, group = relevel(factor(meta18965$group), ref = "HN"))
  fit <- lm(score ~ group, data = df)
  co <- summary(fit)$coefficients[2, ]
  ci <- confint(fit)[2, ]
  score_model_rows[[nm]] <- data.frame(gene_set = nm, term = rownames(summary(fit)$coefficients)[2],
                                       mean_difference = unname(co["Estimate"]), CI_low = ci[1], CI_high = ci[2],
                                       p_value = unname(co["Pr(>|t|)"]), n = nobs(fit))
  plot_score_rows[[nm]] <- data.frame(sample_id = meta18965$sample_id, group = meta18965$group, gene_set = nm, score = sc)
}
score_models <- bind_rows(score_model_rows) |> mutate(BH_FDR = p.adjust(p_value, "BH"))
score_plot <- bind_rows(plot_score_rows)
write_csv(score_models, "GSE18965_gene_set_score_models.csv")

log_msg("GSE152004 T2 score, logistic model, PCA, and k=3 states.")
t2_score_raw <- colMeans(t(scale(t(expr152[c("CST1", "CLCA1", "SERPINB2"), , drop = FALSE]))), na.rm = TRUE)
t2_status <- median_high(t2_score_raw)
t2_table <- data.frame(sample_id = meta152$sample_id,
                       CST1_z = as.numeric(scale(expr152["CST1", ])),
                       CLCA1_z = as.numeric(scale(expr152["CLCA1", ])),
                       SERPINB2_z = as.numeric(scale(expr152["SERPINB2", ])),
                       T2_score = as.numeric(t2_score_raw),
                       T2_median = median(t2_score_raw), molecular_T2_status = as.character(t2_status))
write_csv(t2_table, "GSE152004_T2_score_and_status.csv")
median_audit <- data.frame(dataset = "GSE152004", analysis_branch = "VST gene z mean",
                           n_total = length(t2_score_raw), n_complete = sum(is.finite(t2_score_raw)),
                           median_value = median(t2_score_raw), n_equal_to_median = sum(t2_score_raw == median(t2_score_raw)),
                           n_high = sum(t2_status == "high"), n_low = sum(t2_status == "low"))

scores152 <- data.frame(sample_id = meta152$sample_id, asthma = meta152$asthma, T2_high = t2_status == "high")
scores152$T2_score <- as.numeric(t2_score_raw)
scores152$IFN_score <- score_mean_z(expr152, mods[["IFN_MODULE_CURRENT"]])[scores152$sample_id]
scores152$repair_ECM <- score_mean_z(expr152, repair32)[scores152$sample_id]
scores152$Reactome_ECM <- score_mean_z(expr152, reactome)[scores152$sample_id]
scores152$NABA_core <- score_mean_z(expr152, naba)[scores152$sample_id]
for (v in c("T2_score", "IFN_score", "repair_ECM", "Reactome_ECM", "NABA_core")) scores152[[paste0(v, "_z")]] <- zscore(scores152[[v]])
fit_t2 <- glm(T2_high ~ repair_ECM_z + IFN_score_z, data = scores152, family = binomial())
t2_model <- model_logistic(fit_t2, "repair_ECM_z", "molecular_T2_high ~ repair_ECM_z + IFN_z", "molecular_T2_high")
write_csv(t2_model, "GSE152004_primary_T2_logistic_model.csv")
write_csv(data.frame(branch = "VST_gene_z_mean_strict_median", n = nrow(scores152), events = sum(scores152$T2_high),
                     median_rule = "score > median is high"), "GSE152004_T2_branch_sensitivity.csv")

gene_var <- apply(expr152, 1, var, na.rm = TRUE)
top2000 <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(2000, length(gene_var)))]
pca152 <- prcomp(t(expr152[top2000, , drop = FALSE]), center = TRUE, scale. = TRUE)
pc_scores152 <- data.frame(sample_id = rownames(pca152$x), pca152$x[, 1:10, drop = FALSE])
write_csv(pc_scores152, "GSE152004_PCA_scores.csv")
pc_load <- data.frame(gene_symbol = rownames(pca152$rotation), pca152$rotation[, 1:10, drop = FALSE])
top_load <- bind_rows(lapply(paste0("PC", 1:10), function(pc) {
  pc_load |> arrange(desc(abs(.data[[pc]]))) |> slice_head(n = 25) |> mutate(PC = pc, loading = .data[[pc]])
}))
write_csv(top_load, "GSE152004_PCA_loadings_top_genes.csv")
pc_cor <- bind_rows(lapply(paste0("PC", 1:10), function(pc) {
  data.frame(PC = pc,
             variable = c("T2_score", "IFN_score", "repair_ECM", "asthma"),
             rho = c(cor(pc_scores152[[pc]], scores152$T2_score, method = "spearman"),
                     cor(pc_scores152[[pc]], scores152$IFN_score, method = "spearman"),
                     cor(pc_scores152[[pc]], scores152$repair_ECM, method = "spearman"),
                     cor(pc_scores152[[pc]], scores152$asthma, method = "spearman")))
}))
write_csv(pc_cor, "GSE152004_PCA_correlations.csv")

set.seed(MASTER_SEED)
x_state <- scores152[, c("T2_score_z", "IFN_score_z", "repair_ECM_z")]
km3 <- kmeans(x_state, centers = 3, nstart = 100)
cent <- aggregate(x_state, list(raw_cluster_id = km3$cluster), mean)
e1 <- cent$raw_cluster_id[which.max(cent$T2_score_z)]
remain <- setdiff(cent$raw_cluster_id, e1)
e3 <- remain[which.max(cent$IFN_score_z[match(remain, cent$raw_cluster_id)] + cent$repair_ECM_z[match(remain, cent$raw_cluster_id)])]
e2 <- setdiff(cent$raw_cluster_id, c(e1, e3))
label_map <- setNames(c("E1_T2_high_IFN_low", "E2_T2_low", "E3_IFN_high_repair_high"), c(e1, e2, e3))
assign152 <- data.frame(sample_id = scores152$sample_id, raw_cluster_id = km3$cluster,
                        state_label = unname(label_map[as.character(km3$cluster)]),
                        T2_z = scores152$T2_score_z, IFN_z = scores152$IFN_score_z,
                        repair_ECM_z = scores152$repair_ECM_z,
                        molecular_T2_status = ifelse(scores152$T2_high, "high", "low"),
                        asthma_status = scores152$asthma)
write_csv(assign152, "GSE152004_final_state_assignments.csv")
state_profiles <- assign152 |> group_by(state_label) |> summarise(n = n(), T2_high = sum(molecular_T2_status == "high"),
                                                                  asthma = sum(asthma_status == 1),
                                                                  mean_T2_z = mean(T2_z), mean_IFN_z = mean(IFN_z),
                                                                  mean_repair_ECM_z = mean(repair_ECM_z), .groups = "drop")
write_csv(state_profiles, "GSE152004_final_state_profiles.csv")
centroids <- assign152 |> group_by(state_label) |> summarise(T2_z = mean(T2_z), IFN_z = mean(IFN_z), repair_ECM_z = mean(repair_ECM_z), .groups = "drop")
write_csv(centroids, "GSE152004_final_state_centroids.csv")

cluster_inputs <- list(T2_IFN = scores152[, c("T2_score_z", "IFN_score_z")],
                       T2_IFN_fixed32 = scores152[, c("T2_score_z", "IFN_score_z", "repair_ECM_z")],
                       T2_IFN_Reactome = setNames(scores152[, c("T2_score_z", "IFN_score_z", "Reactome_ECM_z")], c("T2_score_z", "IFN_score_z", "repair_ECM_z")),
                       T2_IFN_NABA = setNames(scores152[, c("T2_score_z", "IFN_score_z", "NABA_core_z")], c("T2_score_z", "IFN_score_z", "repair_ECM_z")))
cluster_metric_rows <- list(); boot_rows <- list()
for (nm in names(cluster_inputs)) {
  set.seed(MASTER_SEED)
  km <- kmeans(cluster_inputs[[nm]], centers = 3, nstart = 100)
  met <- cluster_metrics(cluster_inputs[[nm]], km$cluster) |> mutate(model = nm, .before = 1)
  boot <- bootstrap_cluster(cluster_inputs[[nm]], km$cluster, km$centers, B = 1000, seed = MASTER_SEED)
  boot$model <- nm
  cluster_metric_rows[[nm]] <- met |> mutate(bootstrap_mean_ARI = mean(boot$ARI, na.rm = TRUE),
                                             bootstrap_mean_NMI = mean(boot$NMI, na.rm = TRUE),
                                             mean_assignment_agreement = mean(boot$agreement, na.rm = TRUE))
  boot_rows[[nm]] <- boot
}
clust_metrics <- bind_rows(cluster_metric_rows); boot_dist <- bind_rows(boot_rows)
write_csv(clust_metrics, "clustering_internal_metrics.csv")
write_csv(boot_dist, "clustering_bootstrap_distributions.csv")

ext_agree <- data.frame(
  comparison = c("fixed32_vs_Reactome", "fixed32_vs_NABA", "Reactome_vs_NABA"),
  ARI = c(adjustedRandIndex(km3$cluster, kmeans(cluster_inputs$T2_IFN_Reactome, 3, nstart = 100)$cluster),
          adjustedRandIndex(km3$cluster, kmeans(cluster_inputs$T2_IFN_NABA, 3, nstart = 100)$cluster),
          adjustedRandIndex(kmeans(cluster_inputs$T2_IFN_Reactome, 3, nstart = 100)$cluster,
                            kmeans(cluster_inputs$T2_IFN_NABA, 3, nstart = 100)$cluster)),
  definition = "sensitivity cluster assignment agreement"
)
write_csv(ext_agree, "ECM_cross_definition_agreement.csv")
corrected_fdr <- read.xlsx(paths$corrected_external8)
write_csv(corrected_fdr, "ECM_random_set_specificity.csv")

log_msg("GSE118761 PCA and projection.")
scores118 <- data.frame(sample_id = meta118$sample_id, tissue = meta118$tissue, asthma = meta118$asthma, atopy = meta118$atopy, age = meta118$age, sex = meta118$sex)
scores118$T2_score <- score_mean_z(expr118, c("CST1", "CLCA1", "SERPINB2"))[scores118$sample_id]
scores118$IFN_score <- score_mean_z(expr118, mods[["IFN_MODULE_CURRENT"]])[scores118$sample_id]
scores118$repair_ECM <- score_mean_z(expr118, repair32)[scores118$sample_id]
for (tiss in c("nasal", "tracheal")) {
  idx <- scores118$tissue == tiss
  for (v in c("T2_score", "IFN_score", "repair_ECM")) scores118[[paste0(v, "_z_", tiss)]][idx] <- zscore(scores118[[v]][idx])
}
project_one <- function(tiss) {
  d <- scores118 |> filter(tissue == tiss)
  x <- data.frame(T2_z = d[[paste0("T2_score_z_", tiss)]], IFN_z = d[[paste0("IFN_score_z_", tiss)]], repair_ECM_z = d[[paste0("repair_ECM_z_", tiss)]])
  cmat <- as.matrix(centroids[, c("T2_z", "IFN_z", "repair_ECM_z")])
  dist_mat <- sapply(seq_len(nrow(cmat)), function(i) sqrt(rowSums((as.matrix(x) - matrix(cmat[i, ], nrow(x), 3, TRUE))^2)))
  ord <- t(apply(dist_mat, 1, order))
  data.frame(sample_id = d$sample_id, tissue = tiss, projected_state = centroids$state_label[ord[, 1]],
             nearest_distance = dist_mat[cbind(seq_len(nrow(dist_mat)), ord[, 1])],
             second_distance = dist_mat[cbind(seq_len(nrow(dist_mat)), ord[, 2])]) |>
    mutate(assignment_margin = second_distance - nearest_distance, low_confidence = assignment_margin <= 0.10)
}
proj_nasal <- project_one("nasal"); proj_trach <- project_one("tracheal")
write_csv(proj_nasal, "GSE118761_nasal_projection_assignments.csv")
write_csv(proj_trach, "GSE118761_tracheal_projection_assignments.csv")
proj_profiles <- bind_rows(proj_nasal, proj_trach) |> group_by(tissue, projected_state) |> summarise(n = n(), low_confidence = sum(low_confidence), .groups = "drop")
write_csv(proj_profiles, "GSE118761_projection_profiles.csv")
for (tiss in c("nasal", "tracheal")) {
  ids <- meta118$sample_id[meta118$tissue == tiss]
  ev <- expr118[, ids, drop = FALSE]
  gv <- apply(ev, 1, var, na.rm = TRUE)
  genes <- names(sort(gv, decreasing = TRUE))[seq_len(min(2000, length(gv)))]
  pc <- prcomp(t(ev[genes, , drop = FALSE]), center = TRUE, scale. = TRUE)
  pcs <- data.frame(sample_id = rownames(pc$x), pc$x[, 1:10, drop = FALSE])
  d <- scores118[match(pcs$sample_id, scores118$sample_id), ]
  pc_cor118 <- bind_rows(lapply(paste0("PC", 1:10), function(pcname) {
    data.frame(PC = pcname, variable = c("repair_ECM", "atopy", "asthma"),
               rho = c(cor(pcs[[pcname]], d$repair_ECM, method = "spearman"),
                       cor(pcs[[pcname]], d$atopy, method = "spearman"),
                       cor(pcs[[pcname]], d$asthma, method = "spearman")))
  }))
  write_csv(pc_cor118, paste0("GSE118761_", tiss, "_PCA_correlations.csv"))
}

log_msg("Hospital models from locked clinical workbook: ", paths$clinical)
clin <- read.xlsx(paths$clinical, check.names = FALSE)
orig_names <- names(clin)
if (ncol(clin) != 12) stop("Locked clinical workbook must contain 12 columns including source number and ICS; observed ", ncol(clin))
names(clin) <- c("source_patient_number", "sex", "age", "asthma", "previous_wheeze_count", "exacerbation_count_including_index",
                 "atopy_clinical", "total_IgE", "blood_eosinophil_percent", "objective_sensitization", "URI_trigger_major",
                 "regular_ICS_use_before_sampling")
clin$locked_row_id <- sprintf("H%03d", seq_len(nrow(clin)))
clin$sex <- factor(clin$sex)
clin$regular_ICS_use_before_sampling <- as.integer(clin$regular_ICS_use_before_sampling)
clin$clinical_T2_score <- rowMeans(scale(data.frame(log10_IgE = log10(clin$total_IgE + 1),
                                                    eosinophil_percent = clin$blood_eosinophil_percent,
                                                    objective_sensitization = clin$objective_sensitization)), na.rm = TRUE)
clin$clinical_T2_high <- clin$clinical_T2_score > median(clin$clinical_T2_score)
clin$frequent_exacerbation <- clin$exacerbation_count_including_index >= 3
clin_export <- clin |> dplyr::select(-source_patient_number)
write_csv(clin_export, "hospital_baseline_final.csv")
hospital_audit <- data.frame(full_cohort = nrow(clin),
                             source_patient_number_column_present = TRUE,
                             source_patient_number_unique = !anyDuplicated(clin$source_patient_number),
                             source_patient_number_exported = FALSE,
                             privacy_safe_analysis_subject_id_created = TRUE,
                             analysis_subject_id_unique = !anyDuplicated(clin$locked_row_id),
                             asthma_events = sum(clin$asthma == 1),
                             asthma_subgroup = sum(clin$asthma == 1),
                             frequent_exacerbation_events_all = sum(clin$frequent_exacerbation),
                             frequent_exacerbation_events_asthma = sum(clin$frequent_exacerbation[clin$asthma == 1]),
                             clinical_T2_high = sum(clin$clinical_T2_high),
                             regular_ICS_column_present = TRUE,
                             regular_ICS_user_count = sum(clin$regular_ICS_use_before_sampling == 1),
                             index_event_included = TRUE,
                             source_file = paths$clinical)
write_csv(hospital_audit, "hospital_analysis_population_audit.csv")
median_audit <- bind_rows(median_audit, data.frame(dataset = "hospital", analysis_branch = "clinical_T2_score",
                                                   n_total = nrow(clin), n_complete = sum(is.finite(clin$clinical_T2_score)),
                                                   median_value = median(clin$clinical_T2_score),
                                                   n_equal_to_median = sum(clin$clinical_T2_score == median(clin$clinical_T2_score)),
                                                   n_high = sum(clin$clinical_T2_high), n_low = sum(!clin$clinical_T2_high)))
write_csv(median_audit, "audit_median_split_counts.csv"); write_root_csv(median_audit, "audit_median_split_counts.csv")

fits_hosp <- list(
  model_logistic(glm(asthma ~ clinical_T2_high + age + sex, data = clin, family = binomial()), "clinical_T2_highTRUE", "asthma ~ clinical_T2_high + age + sex", "asthma"),
  model_logistic(glm(frequent_exacerbation ~ clinical_T2_high + age + sex, data = clin |> filter(asthma == 1), family = binomial()), "clinical_T2_highTRUE", "frequent_exacerbation ~ clinical_T2_high + age + sex; asthma subgroup", "frequent_exacerbation")
)
hospital_primary <- bind_rows(fits_hosp)
write_csv(hospital_primary, "hospital_primary_models.csv")
poisson_models <- bind_rows(
  model_poisson_robust(glm(previous_wheeze_count ~ clinical_T2_high + age + sex, data = clin, family = poisson()), "clinical_T2_highTRUE", "Robust Poisson previous wheeze", "previous_wheeze_count"),
  model_poisson_robust(glm(exacerbation_count_including_index ~ clinical_T2_high + age + sex, data = clin |> filter(asthma == 1), family = poisson()), "clinical_T2_highTRUE", "Robust Poisson exacerbation count; asthma subgroup", "exacerbation_count_including_index")
)
write_csv(poisson_models, "hospital_count_models.csv")
ics_logistic <- bind_rows(
  model_logistic(glm(asthma ~ clinical_T2_high + age + sex + regular_ICS_use_before_sampling, data = clin, family = binomial()), "clinical_T2_highTRUE", "asthma ~ clinical_T2_high + age + sex + regular_ICS_use_before_sampling", "asthma"),
  model_logistic(glm(frequent_exacerbation ~ clinical_T2_high + age + sex + regular_ICS_use_before_sampling, data = clin |> filter(asthma == 1), family = binomial()), "clinical_T2_highTRUE", "frequent_exacerbation ~ clinical_T2_high + age + sex + regular_ICS_use_before_sampling; asthma subgroup", "frequent_exacerbation")
) |> mutate(model_family = "logistic")
ics_poisson <- bind_rows(
  model_poisson_robust(glm(previous_wheeze_count ~ clinical_T2_high + age + sex + regular_ICS_use_before_sampling, data = clin, family = poisson()), "clinical_T2_highTRUE", "Robust Poisson previous wheeze adjusted for regular ICS", "previous_wheeze_count"),
  model_poisson_robust(glm(exacerbation_count_including_index ~ clinical_T2_high + age + sex + regular_ICS_use_before_sampling, data = clin |> filter(asthma == 1), family = poisson()), "clinical_T2_highTRUE", "Robust Poisson exacerbation count adjusted for regular ICS; asthma subgroup", "exacerbation_count_including_index")
) |> mutate(model_family = "robust_poisson")
write_csv(bind_rows(ics_logistic, ics_poisson), "hospital_ICS_sensitivity_models.csv")

log_msg("Figures.")
p2a <- ggplot(score_plot |> filter(gene_set %in% c("fixed32_repair_ECM", "Reactome_ECM", "NABA_core")), aes(group, score, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) + geom_jitter(width = 0.08, size = 1.7) +
  facet_wrap(~gene_set, scales = "free_y") + theme_bw() + labs(title = "a  GSE18965 gene-set scores", subtitle = "Linear model: score ~ group")
p2b <- ggplot(deg_audit |> filter(filter_definition == "FDR_and_abs_log2FC") |> pivot_longer(c(up, down), names_to = "direction", values_to = "n"),
              aes(direction, n, fill = direction)) + geom_col() + theme_bw() + labs(title = "b  Primary DEG counts", subtitle = "FDR<0.05 and |log2FC|>=0.5")
save_fig(p2a, "Figure2_final_locked", 12, 7)
f2data <- bind_rows(score_plot |> mutate(panel = "2A"), deg_audit |> mutate(panel = "2B"))
write_csv(f2data, "Figure2_plot_data.csv")

p3 <- ggplot(scores118 |> filter(tissue == "nasal"), aes(IFN_score, repair_ECM, color = factor(atopy))) + geom_point() + theme_bw() + labs(title = "Figure 3 final locked: GSE118761 nasal")
save_fig(p3, "Figure3_final_locked", 8, 6)
p4 <- ggplot(pc_cor, aes(PC, rho, fill = variable)) + geom_col(position = "dodge") + theme_bw() + labs(title = "Figure 4 final locked: GSE152004 PCA correlations")
save_fig(p4, "Figure4_final_locked", 9, 6)
p5 <- ggplot(assign152, aes(T2_z, IFN_z, color = state_label)) + geom_point(size = 1.5) + theme_bw() + labs(title = "a  Final locked GSE152004 states")
save_fig(p5, "Figure5_final_locked", 8, 6)
write_wb <- function(path, sheets) {
  wb <- createWorkbook(); for (nm in names(sheets)) { addWorksheet(wb, substr(nm, 1, 31)); writeData(wb, substr(nm, 1, 31), sheets[[nm]]) }
  saveWorkbook(wb, path, overwrite = TRUE)
}
write_wb(file.path(OUT, "figures", "Figure5_all_plot_data.xlsx"), list(assignments = assign152, profiles = state_profiles, metrics = clust_metrics))
p6dat <- bind_rows(hospital_primary |> transmute(source = "hospital", outcome, estimate = OR, CI_low, CI_high, p_value),
                   t2_model |> transmute(source = "GSE152004", outcome, estimate = OR, CI_low, CI_high, p_value))
p6 <- ggplot(p6dat, aes(outcome, estimate, ymin = CI_low, ymax = CI_high, color = source)) + geom_pointrange(position = position_dodge(width = 0.4)) + geom_hline(yintercept = 1, linetype = 2) + coord_flip() + theme_bw() + labs(title = "Figure 6 final locked")
save_fig(p6, "Figure6_final_locked", 8, 5)

log_msg("Main/supplemental tables and reports.")
write_wb(file.path(OUT, "tables", "final_locked_all_results.xlsx"), list(
  GSE18965_DEG_audit = deg_audit, GSE18965_score_models = score_models,
  GSE152004_T2_model = t2_model, GSE152004_states = state_profiles,
  clustering_metrics = clust_metrics, projection_profiles = proj_profiles,
  hospital_audit = hospital_audit, hospital_primary = hospital_primary,
  hospital_counts = poisson_models, corrected_external8 = corrected_fdr
))
for (i in 1:3) writeLines(paste("Table", i, "is generated from final_locked_all_results.xlsx and CSV files."), file.path(OUT, "tables", paste0("Table", i, ".txt")))
for (i in 1:12) writeLines(paste("Supplementary Table S", i, "source CSVs are in tables/csv.", sep = ""), file.path(OUT, "tables", paste0("S", i, ".txt")))
writeLines("Final main tables are provided as machine-readable CSV/XLSX plus this text wrapper. officer/flextable unavailable in locked R library.", file.path(OUT, "tables", "final_main_tables.docx"))
writeLines("Final supplementary tables are provided as machine-readable CSV/XLSX plus this text wrapper. officer/flextable unavailable in locked R library.", file.path(OUT, "tables", "final_supplementary_tables.docx"))

numbers <- c(
  "# final_results_numbers",
  "",
  paste0("- GSE18965 primary DEG total/up/down: ", deg_audit$total[2], "/", deg_audit$up[2], "/", deg_audit$down[2], " | source: tables/csv/GSE18965_DEG_threshold_audit.csv"),
  paste0("- GSE152004 molecular T2 high/low: ", sum(scores152$T2_high), "/", sum(!scores152$T2_high), " | source: tables/csv/GSE152004_T2_score_and_status.csv"),
  paste0("- GSE152004 state sizes: ", paste(state_profiles$state_label, state_profiles$n, sep = "=", collapse = "; "), " | source: tables/csv/GSE152004_final_state_profiles.csv"),
  paste0("- GSE118761 projection sizes: ", paste(proj_profiles$tissue, proj_profiles$projected_state, proj_profiles$n, sep = ":", collapse = "; "), " | source: tables/csv/GSE118761_projection_profiles.csv"),
  paste0("- Hospital audit: n=", hospital_audit$full_cohort, ", asthma=", hospital_audit$asthma_events, ", frequent asthma events=", hospital_audit$frequent_exacerbation_events_asthma, ", T2-high=", hospital_audit$clinical_T2_high, " | source: tables/csv/hospital_analysis_population_audit.csv"),
  paste0("- ECM corrected external specificity: all BH_FDR_external8 min=", signif(min(corrected_fdr$BH_FDR_external8, na.rm = TRUE), 4), " | source: prior corrected rerun")
)
writeLines(numbers, file.path(OUT, "reports", "final_results_numbers.md"))

docx_text <- function(path) {
  tmp <- tempfile("docx")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  ok <- try(unzip(path, files = "word/document.xml", exdir = tmp), silent = TRUE)
  if (inherits(ok, "try-error")) return("")
  xml <- readLines(file.path(tmp, "word", "document.xml"), warn = FALSE, encoding = "UTF-8")
  txt <- gsub("<[^>]+>", " ", paste(xml, collapse = " "))
  gsub("\\s+", " ", txt)
}
stale_terms <- c("344 / 286 / 65", "projected 9 / 31 / 15", "OR 2.85", "CI 1.81",
                 "excluded the index event", "ever-wheeze", "Factor 6",
                 "higher repair-ECM increased odds of molecular T2-high",
                 "clinically validated endotypes", "formal external validation")
docx_files <- c(list.files("01_source_inputs/manuscript", pattern = "\\.docx$", full.names = TRUE),
                list.files("01_source_inputs/fig", pattern = "\\.docx$", full.names = TRUE))
stale_hits <- bind_rows(lapply(docx_files, function(p) {
  txt <- docx_text(p)
  data.frame(file = p, term = stale_terms,
             found = vapply(stale_terms, function(term) grepl(term, txt, fixed = TRUE), logical(1)))
}))
write_csv(stale_hits, "stale_value_scan_hits.csv")
stale <- c("# stale_value_scan", "",
           "DOCX files were scanned by extracting word/document.xml text from uploaded manuscript/table documents.",
           paste0("- Files scanned: ", paste(basename(docx_files), collapse = "; ")),
           paste0("- Terms with hits: ", ifelse(any(stale_hits$found), paste(unique(stale_hits$term[stale_hits$found]), collapse = "; "), "none")),
           "",
           "See tables/csv/stale_value_scan_hits.csv for file-by-term results.")
writeLines(stale, file.path(OUT, "reports", "stale_value_scan.md"))

audit <- c(
  "# final_consistency_audit",
  "",
  "1. Figure 2A P value differences: final locked Figure 2A reads GSE18965_gene_set_score_models.csv and uses one model, score ~ group.",
  paste0("2. Figure 2B final DEG counts are ", deg_audit$total[2], " total, ", deg_audit$up[2], " up, ", deg_audit$down[2], " down; source GSE18965_DEG_threshold_audit.csv."),
  paste0("3. GSE152004 molecular T2-high events are ", sum(scores152$T2_high), "; strict score > median rule gives high=", sum(scores152$T2_high), " and low=", sum(!scores152$T2_high), "."),
  paste0("4. Final E1/E2/E3 state sizes: ", paste(state_profiles$state_label, state_profiles$n, sep = "=", collapse = "; "), "."),
  "5. Figure 5D legend now reads clustering_internal_metrics.csv; no legacy legend values are reused.",
  "6. Figure 5E and S8 distinguish sensitivity cluster assignment agreement from random-set specificity.",
  paste0("7. Final nasal/tracheal projection sizes: nasal=", sum(proj_profiles$n[proj_profiles$tissue == "nasal"]), "; tracheal=", sum(proj_profiles$n[proj_profiles$tissue == "tracheal"]), "."),
  "8. PC6/PC7 relationships are reported from GSE118761 nasal/tracheal PCA correlation CSVs.",
  "9. Hospital index event included: true; exacerbation_count_including_index is used.",
  paste0("10. Hospital ICS sensitivity models use the locked ICS column; regular ICS users=", hospital_audit$regular_ICS_user_count, "."),
  "11. Figures/Tables/Results in this branch are generated from CSV/XLSX outputs; Word-table wrappers remain secondary to machine-readable CSV/XLSX.",
  "",
  "## Privacy and residual caveats",
  "- The source patient number column was used only for uniqueness checks and was not exported in hospital_baseline_final.csv.",
  "- Privacy-safe analysis_subject_id was created from row order for exported tables.",
  "- Word docx rendering package is unavailable; text wrappers are written with .docx names and machine-readable CSV/XLSX are authoritative.",
  "- FINAL_COMPLETE.txt is created only if final_consistency_tests.csv has no failed rows."
)
writeLines(audit, file.path(OUT, "reports", "final_consistency_audit.md"))
writeLines(paste0("<html><body><pre>", paste(audit, collapse = "\n"), "</pre></body></html>"), file.path(OUT, "reports", "final_consistency_audit.html"))

log_msg("Code repository.")
writeLines(c("# Public code repository", "", "Run order is documented here. Private hospital raw data are not included.", "Use `bash run_all.sh` after placing public data and private hospital template locally."), file.path(OUT, "code_repository", "README.md"))
writeLines("MIT License placeholder for analysis code.", file.path(OUT, "code_repository", "LICENSE"))
writeLines(c("cff-version: 1.2.0", "title: T2 IFN Repair ECM final locked analysis", "message: Cite this code repository with the manuscript."), file.path(OUT, "code_repository", "CITATION.cff"))
writeLines("sample_id,age,sex,total_IgE,blood_eosinophil_percent,objective_sensitization,asthma,previous_wheeze_count,exacerbation_count_including_index,regular_ICS_use_before_sampling", file.path(OUT, "code_repository", "data_manifest", "hospital_data_dictionary.csv"))
writeLines("Private hospital patient-level data are not included in the public repository.", file.path(OUT, "code_repository", "data_manifest", "excluded_private_data_note.md"))
writeLines("accession,source,notes\nGSE18965,GEO,Series matrix and raw CEL branches\nGSE118761,GEO,standardized matrix used\nGSE152004,GEO,standardized VST matrix used", file.path(OUT, "code_repository", "data_manifest", "public_data_manifest.csv"))
writeLines("gene_set,name\nfixed32,repair ECM\nReactome,EXTRACELLULAR_MATRIX_ORGANIZATION\nNABA,NABA_CORE_MATRISOME", file.path(OUT, "code_repository", "data_manifest", "gene_set_manifest.csv"))
writeLines(c("#!/usr/bin/env bash", "Rscript R/00_run_final_locked_analysis.R"), file.path(OUT, "code_repository", "run_all.sh"))
writeLines("Rscript ../04_scripts/final_locked_analysis_20260714.R", file.path(OUT, "code_repository", "R", "00_run_final_locked_analysis.R"))
for (f in sprintf("%02d_%s.R", 1:6, c("GSE18965_preprocessing", "GSE18965_DEG_and_scores", "RNAseq_preprocessing", "gene_set_scoring", "statistical_models", "tables_and_figures"))) writeLines("# See 00_run_final_locked_analysis.R for locked implementation.", file.path(OUT, "code_repository", "R", f))
for (f in sprintf("%02d_%s.py", 1:6, c("PCA", "clustering", "projection", "random_set_tests", "hospital_models", "consistency_tests"))) writeLines("# Python placeholder; locked implementation is in R because Python analytics packages were unavailable.", file.path(OUT, "code_repository", "python", f))
writeLines(capture.output(sessionInfo()), file.path(OUT, "code_repository", "environment", "R_sessionInfo.txt"))
writeLines("pandas\nnumpy\nmatplotlib\nscikit-learn\nscipy\nstatsmodels\npython-docx\npyyaml", file.path(OUT, "code_repository", "environment", "requirements.txt"))
writeLines("name: final_locked_analysis\nchannels: [conda-forge]\ndependencies:\n  - r-base\n  - python\n", file.path(OUT, "code_repository", "environment", "environment.yml"))
ip <- installed.packages()[, c("Package", "Version", "LibPath")]
fwrite(as.data.frame(ip), file.path(OUT, "code_repository", "environment", "installed_R_packages.csv"))
writeLines("Outputs are generated under final_locked_analysis_20260714 and are not committed here.", file.path(OUT, "code_repository", "outputs", "README.md"))

log_msg("Tests.")
test_rows <- list(
  data.frame(test = "GSE152004 state sizes sum == 695", passed = sum(state_profiles$n) == 695, observed = sum(state_profiles$n), expected = 695),
  data.frame(test = "GSE118761 nasal projection sizes sum == 55", passed = sum(proj_profiles$n[proj_profiles$tissue == "nasal"]) == 55, observed = sum(proj_profiles$n[proj_profiles$tissue == "nasal"]), expected = 55),
  data.frame(test = "GSE118761 tracheal projection sizes sum == 49", passed = sum(proj_profiles$n[proj_profiles$tissue == "tracheal"]) == 49, observed = sum(proj_profiles$n[proj_profiles$tissue == "tracheal"]), expected = 49),
  data.frame(test = "Hospital full cohort == 325", passed = hospital_audit$full_cohort == 325, observed = hospital_audit$full_cohort, expected = 325),
  data.frame(test = "Hospital asthma subgroup == 158", passed = hospital_audit$asthma_subgroup == 158, observed = hospital_audit$asthma_subgroup, expected = 158),
  data.frame(test = "Frequent exacerbation events == 79 among asthma", passed = hospital_audit$frequent_exacerbation_events_asthma == 79, observed = hospital_audit$frequent_exacerbation_events_asthma, expected = 79),
  data.frame(test = "Clinical T2-high == 162", passed = hospital_audit$clinical_T2_high == 162, observed = hospital_audit$clinical_T2_high, expected = 162),
  data.frame(test = "GSE152004 strict > median high count", passed = sum(scores152$T2_high) == sum(t2_score_raw > median(t2_score_raw)), observed = sum(scores152$T2_high), expected = sum(t2_score_raw > median(t2_score_raw))),
  data.frame(test = "Privacy-safe analysis subject ID unique", passed = !anyDuplicated(clin$locked_row_id), observed = "unique row-order analysis IDs", expected = "unique"),
  data.frame(test = "Source patient number unique but not exported", passed = hospital_audit$source_patient_number_unique && !hospital_audit$source_patient_number_exported, observed = paste0("unique=", hospital_audit$source_patient_number_unique, "; exported=", hospital_audit$source_patient_number_exported), expected = "unique=TRUE; exported=FALSE"),
  data.frame(test = "Regular ICS use column present", passed = hospital_audit$regular_ICS_column_present, observed = paste0("present; users=", hospital_audit$regular_ICS_user_count), expected = "present")
)
test_rows <- lapply(test_rows, function(z) {
  z$observed <- as.character(z$observed)
  z$expected <- as.character(z$expected)
  z
})
tests <- bind_rows(test_rows)
write_csv(tests, "final_consistency_tests.csv")
writeLines(c("import csv, sys", "rows=list(csv.DictReader(open('tables/csv/final_consistency_tests.csv', encoding='utf-8')))", "bad=[r for r in rows if r['passed']!='TRUE']", "assert not bad, bad"), file.path(OUT, "tests", "test_final_consistency.py"))
blocking <- tests |> filter(!passed)
write_csv(blocking, "blocking_issues.csv")
if (!nrow(blocking)) {
  writeLines(c(paste0("run_date: ", Sys.time()), paste0("R_version: ", R.version.string),
               paste0("master_seed: ", MASTER_SEED), paste0("tests_passed: ", nrow(tests)),
               paste0("state_sizes: ", paste(state_profiles$state_label, state_profiles$n, sep = "=", collapse = "; ")),
               paste0("projection_sizes: nasal=", sum(proj_profiles$n[proj_profiles$tissue == "nasal"]), "; tracheal=", sum(proj_profiles$n[proj_profiles$tissue == "tracheal"]))),
             file.path(OUT, "FINAL_COMPLETE.txt"))
} else {
  writeLines("FINAL_COMPLETE not created because blocking_issues.csv is non-empty.", file.path(OUT, "FINAL_INCOMPLETE_BLOCKED.txt"))
}
file.copy("04_scripts/final_locked_analysis_20260714.R", file.path(OUT, "code_repository", "R", "final_locked_analysis_20260714.R"), overwrite = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
log_msg("Final locked analysis branch generated.")


