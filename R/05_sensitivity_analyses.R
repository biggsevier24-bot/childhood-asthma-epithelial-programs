options(stringsAsFactors = FALSE)
if (file.exists("config/use_local_r_libs.R")) source("config/use_local_r_libs.R")
set.seed(20260626)

required <- c(
  "data.table", "dplyr", "tidyr", "mclust", "cluster", "openxlsx",
  "GSVA", "GSEABase"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(mclust)
  library(cluster)
  library(openxlsx)
  library(GSVA)
  library(GSEABase)
})

CONFIG <- list(
  seed = 20260626L,
  k = 3L,
  full_nstart = 500L,
  permutation_n = 10000L,
  permutation_boot_reps = 100L,
  descriptive_boot_reps = 1000L,
  checkpoint_every = 500L,
  low_confidence_margin = 0.10,
  out_dir = "03_result_directories/ECM_third_axis_corrected_rerun",
  zip_path = "02_final_result_zips/ECM_third_axis_corrected_rerun_results.zip",
  expr152 = "03_result_directories/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE152004_VST_gene_expression.csv",
  meta152 = "03_result_directories/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE152004_metadata.csv",
  expr118 = "03_result_directories/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE118761_VST_gene_expression.csv",
  meta118 = "03_result_directories/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE118761_metadata.csv",
  random_members = "03_result_directories/ECM_third_axis_validation/04_random_gene_set_permutation/random_gene_set_membership.csv.gz",
  previous_stats = "03_result_directories/ECM_third_axis_validation/04_random_gene_set_permutation/permutation_statistics_all_10000.csv.gz",
  repair32 = "01_source_inputs/GSE18965_Codex_full_package/repair_ECM_32_genes.csv",
  modules_gmt = "01_source_inputs/external_ECM_user_inputs/T2_IFN_cluster_metadata_files/article_T2_IFN_Repair_custom_modules.gmt",
  gene_sets_external = "03_result_directories/External_ECM_multicohort_analysis/02_gene_sets/external_ECM_gene_sets_original.csv",
  prompt = "Codex_prompt_ECM_third_axis_corrected_rerun.txt",
  template = "ECM_third_axis_corrected_rerun_template.R"
)

dirs <- c(
  "00_logs", "01_input_audit", "02_corrected_bootstrap",
  "03_external_specificity_permutation/checkpoints",
  "04_fair_mean_z_ssGSEA_models", "05_corrected_projection",
  "06_strict_reclassification", "07_report", "08_code"
)

reset_dir <- function(path) {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  target <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (file.exists(path) && !startsWith(target, paste0(root, "/"))) {
    stop("Refusing to delete outside workspace: ", target)
  }
  if (file.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

dir.create(dirname(CONFIG$zip_path), recursive = TRUE, showWarnings = FALSE)
reset_dir(CONFIG$out_dir)
for (d in dirs) dir.create(file.path(CONFIG$out_dir, d), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(CONFIG$out_dir, "00_logs", "run_log.txt")
perm_log <- file.path(CONFIG$out_dir, "03_external_specificity_permutation", "permutation_run_log.txt")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_perm <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n", file = perm_log, append = TRUE)
}

add_sheet <- function(wb, name, df) {
  nm <- substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", name), 1, 31)
  addWorksheet(wb, nm)
  writeData(wb, nm, df)
  freezePane(wb, nm, firstRow = TRUE)
  if (ncol(df)) setColWidths(wb, nm, cols = 1:min(ncol(df), 25), widths = "auto")
}
write_wb <- function(path, sheets) {
  wb <- createWorkbook()
  for (nm in names(sheets)) add_sheet(wb, nm, sheets[[nm]])
  saveWorkbook(wb, path, overwrite = TRUE)
}

md5_or_na <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

read_expr <- function(path) {
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  gene_col <- names(x)[1]
  genes <- toupper(trimws(as.character(x[[gene_col]])))
  x[[gene_col]] <- NULL
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  rownames(m) <- genes
  m <- m[rownames(m) != "" & !is.na(rownames(m)), , drop = FALSE]
  if (any(duplicated(rownames(m)))) m <- rowsum(m, group = rownames(m), reorder = FALSE)
  m
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  for (ln in lines) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(parts) >= 3) out[[parts[1]]] <- toupper(trimws(parts[-c(1, 2)]))
  }
  out
}

score_mean_z <- function(expr, genes, min_genes = 3L, min_coverage = 0.05) {
  genes <- unique(toupper(trimws(genes)))
  present <- intersect(genes, rownames(expr))
  coverage <- length(present) / max(1L, length(genes))
  if (length(present) < min_genes || coverage < min_coverage) {
    return(list(score = setNames(rep(NA_real_, ncol(expr)), colnames(expr)),
                genes_used = present, n_used = length(present), coverage = coverage, estimable = FALSE))
  }
  m <- expr[present, , drop = FALSE]
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  keep <- is.finite(sdv) & sdv > 0
  m <- m[keep, , drop = FALSE]
  if (!nrow(m)) {
    return(list(score = setNames(rep(NA_real_, ncol(expr)), colnames(expr)),
                genes_used = character(), n_used = 0L, coverage = 0, estimable = FALSE))
  }
  zm <- t(scale(t(m)))
  list(score = setNames(colMeans(zm, na.rm = TRUE), colnames(expr)),
       genes_used = rownames(zm), n_used = nrow(zm), coverage = nrow(zm) / max(1L, length(genes)), estimable = TRUE)
}

kmeans_full <- function(x, seed = CONFIG$seed) {
  set.seed(seed)
  stats::kmeans(as.matrix(x), centers = CONFIG$k, nstart = CONFIG$full_nstart, iter.max = 100)
}

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], all_permutations(x[-i]))))
}

hungarian_by_centroid_distance <- function(reference_centers, candidate_centers) {
  perms <- all_permutations(seq_len(nrow(candidate_centers)))
  loss <- apply(perms, 1, function(p) sum((reference_centers - candidate_centers[p, , drop = FALSE])^2))
  best <- perms[which.min(loss), ]
  list(order = best, centers = candidate_centers[best, , drop = FALSE], loss = min(loss))
}

assign_to_centers <- function(x, centers) {
  x <- as.matrix(x)
  centers <- as.matrix(centers)
  d <- sapply(seq_len(nrow(centers)), function(j) {
    sqrt(rowSums((x - matrix(centers[j, ], nrow(x), ncol(x), byrow = TRUE))^2))
  })
  max.col(-d, ties.method = "first")
}

normalized_mutual_information <- function(a, b) {
  a <- factor(a)
  b <- factor(b)
  tab <- table(a, b)
  pxy <- tab / sum(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)
  mi <- 0
  for (i in seq_len(nrow(pxy))) {
    for (j in seq_len(ncol(pxy))) {
      if (pxy[i, j] > 0) mi <- mi + pxy[i, j] * log(pxy[i, j] / (px[i] * py[j]))
    }
  }
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx <= 0 || hy <= 0) return(NA_real_)
  mi / sqrt(hx * hy)
}

best_match <- function(reference, candidate) {
  reference <- as.character(reference)
  candidate <- as.character(candidate)
  refs <- sort(unique(reference))
  cands <- sort(unique(candidate))
  tab <- table(factor(reference, levels = refs), factor(candidate, levels = cands))
  perms <- all_permutations(seq_along(cands))
  score <- apply(perms, 1, function(p) sum(tab[cbind(seq_along(refs), p)]))
  p <- perms[which.max(score), ]
  map <- setNames(refs, cands[p])
  unname(map[candidate])
}

cluster_metrics <- function(x, cl) {
  x <- as.matrix(x)
  cl <- factor(cl)
  tab <- table(cl)
  if (nlevels(cl) != CONFIG$k) {
    return(data.frame(mean_silhouette = NA_real_, median_silhouette = NA_real_,
                      calinski_harabasz = NA_real_, davies_bouldin = NA_real_,
                      centroid_min_distance = NA_real_, min_cluster_size = min(tab),
                      cluster_collapse = TRUE))
  }
  sil <- cluster::silhouette(as.integer(cl), stats::dist(x))[, "sil_width"]
  lev <- levels(cl)
  centers <- do.call(rbind, lapply(lev, function(g) colMeans(x[cl == g, , drop = FALSE])))
  overall <- colMeans(x)
  within <- sum(vapply(seq_along(lev), function(i) {
    g <- lev[i]
    sum(rowSums((x[cl == g, , drop = FALSE] - matrix(centers[i, ], sum(cl == g), ncol(x), byrow = TRUE))^2))
  }, numeric(1)))
  between <- sum(vapply(seq_along(lev), function(i) {
    g <- lev[i]
    sum(cl == g) * sum((centers[i, ] - overall)^2)
  }, numeric(1)))
  ch <- (between / (CONFIG$k - 1)) / (within / (nrow(x) - CONFIG$k))
  scat <- vapply(seq_along(lev), function(i) {
    g <- lev[i]
    mean(sqrt(rowSums((x[cl == g, , drop = FALSE] - matrix(centers[i, ], sum(cl == g), ncol(x), byrow = TRUE))^2)))
  }, numeric(1))
  cdist <- as.matrix(dist(centers))
  db <- mean(vapply(seq_len(CONFIG$k), function(i) max((scat[i] + scat[-i]) / cdist[i, -i]), numeric(1)))
  data.frame(mean_silhouette = mean(sil), median_silhouette = median(sil),
             calinski_harabasz = ch, davies_bouldin = db,
             centroid_min_distance = min(cdist[upper.tri(cdist)]),
             min_cluster_size = min(tab), cluster_collapse = min(tab) < 5L)
}

make_bootstrap_indices <- function(n, reps, seed = CONFIG$seed + 800000L) {
  set.seed(seed)
  lapply(seq_len(reps), function(i) sample.int(n, size = n, replace = TRUE))
}

corrected_bootstrap_refit <- function(x, reference_cluster, reference_centers, bootstrap_indices, min_cluster_size = 5L) {
  x <- as.matrix(x)
  reference_cluster <- as.integer(reference_cluster)
  out <- vector("list", length(bootstrap_indices))
  for (b in seq_along(bootstrap_indices)) {
    idx <- bootstrap_indices[[b]]
    xb <- x[idx, , drop = FALSE]
    fit <- try(stats::kmeans(xb, centers = reference_centers, iter.max = 100, algorithm = "Lloyd"), silent = TRUE)
    if (inherits(fit, "try-error") || length(unique(fit$cluster)) != CONFIG$k) {
      out[[b]] <- data.frame(bootstrap = b, ARI = NA_real_, NMI = NA_real_, agreement = NA_real_,
                             min_cluster_size = if (inherits(fit, "try-error")) 0L else min(table(fit$cluster)),
                             collapse = TRUE)
      next
    }
    matched <- hungarian_by_centroid_distance(reference_centers, fit$centers)
    predicted_all <- assign_to_centers(x, matched$centers)
    tab <- table(predicted_all)
    collapse <- length(tab) != CONFIG$k || min(tab) < min_cluster_size
    if (collapse) {
      out[[b]] <- data.frame(bootstrap = b, ARI = NA_real_, NMI = NA_real_, agreement = NA_real_,
                             min_cluster_size = if (length(tab)) min(tab) else 0L, collapse = TRUE)
    } else {
      out[[b]] <- data.frame(bootstrap = b,
                             ARI = mclust::adjustedRandIndex(reference_cluster, predicted_all),
                             NMI = normalized_mutual_information(reference_cluster, predicted_all),
                             agreement = mean(reference_cluster == predicted_all),
                             min_cluster_size = min(tab), collapse = FALSE)
    }
  }
  bind_rows(out)
}

summarise_bootstrap <- function(x, model_name) {
  valid <- x[!x$collapse & is.finite(x$ARI), , drop = FALSE]
  data.frame(model = model_name, reps_requested = nrow(x), reps_valid = nrow(valid),
             collapse_n = sum(x$collapse), collapse_rate = mean(x$collapse),
             mean_ARI = mean(valid$ARI), median_ARI = median(valid$ARI),
             sd_ARI = sd(valid$ARI), mean_NMI = mean(valid$NMI),
             mean_agreement = mean(valid$agreement), min_ARI = min(valid$ARI),
             max_ARI = max(valid$ARI))
}

run_bootstrap_unit_tests <- function() {
  set.seed(123)
  x <- rbind(matrix(rnorm(150, -1.2, 0.9), ncol = 3),
             matrix(rnorm(150, 0, 0.9), ncol = 3),
             matrix(rnorm(150, 1.2, 0.9), ncol = 3))
  ref <- kmeans(x, centers = 3, nstart = 100)
  idx <- make_bootstrap_indices(nrow(x), 50, seed = 777)
  boot <- corrected_bootstrap_refit(x, ref$cluster, ref$centers, idx)
  same_ari <- adjustedRandIndex(ref$cluster, ref$cluster)
  shuffled_ari <- adjustedRandIndex(ref$cluster, sample(ref$cluster))
  list(same_cluster_ARI = same_ari,
       shuffled_cluster_ARI = shuffled_ari,
       bootstrap_ARI_variance = var(boot$ARI, na.rm = TRUE),
       bootstrap_all_one = all(boot$ARI == 1, na.rm = TRUE),
       passed = abs(same_ari - 1) < 1e-12 && abs(shuffled_ari) < 0.2 &&
         sd(boot$ARI, na.rm = TRUE) > 0 && !all(boot$ARI == 1, na.rm = TRUE),
       bootstrap = boot)
}

projection_margin <- function(train_x, train_cluster, valid_x) {
  train_x <- as.matrix(train_x)
  valid_x <- as.matrix(valid_x)
  centers <- do.call(rbind, lapply(sort(unique(train_cluster)), function(g) colMeans(train_x[train_cluster == g, , drop = FALSE])))
  d <- sapply(seq_len(nrow(centers)), function(j) {
    sqrt(rowSums((valid_x - matrix(centers[j, ], nrow(valid_x), ncol(valid_x), byrow = TRUE))^2))
  })
  ord <- t(apply(d, 1, order))
  nearest <- d[cbind(seq_len(nrow(d)), ord[, 1])]
  second <- d[cbind(seq_len(nrow(d)), ord[, 2])]
  margin <- second - nearest
  data.frame(projected_cluster = ord[, 1], nearest_distance = nearest,
             second_distance = second, assignment_margin = margin,
             low_confidence = margin <= CONFIG$low_confidence_margin)
}

empirical_p_larger <- function(random_values, true_value) {
  random_values <- random_values[is.finite(random_values)]
  n <- length(random_values)
  p <- (1 + sum(random_values >= true_value)) / (1 + n)
  data.frame(true_value = true_value, random_n = n, random_mean = mean(random_values),
             random_sd = sd(random_values),
             q2.5 = unname(quantile(random_values, 0.025)),
             q50 = unname(quantile(random_values, 0.50)),
             q97.5 = unname(quantile(random_values, 0.975)),
             empirical_p = p, monte_carlo_se = sqrt(p * (1 - p) / (n + 1)),
             true_percentile = mean(random_values <= true_value))
}

eta2_oneway <- function(y, group) {
  y <- as.numeric(y)
  group <- factor(group)
  ok <- is.finite(y) & !is.na(group)
  y <- y[ok]
  group <- droplevels(group[ok])
  if (length(unique(group)) < 2) return(NA_real_)
  grand <- mean(y)
  ss_between <- sum(vapply(levels(group), function(g) sum(group == g) * (mean(y[group == g]) - grand)^2, numeric(1)))
  ss_total <- sum((y - grand)^2)
  ss_between / ss_total
}

fit_fair_logistic_pair <- function(metadata, outcome, mean_z_score, ssgsea_score, covariates, dataset, tissue, ECM_definition) {
  d <- metadata
  d$mean_z_score <- zscore(mean_z_score)
  d$ssgsea_score <- zscore(ssgsea_score)
  if ("sex" %in% names(d)) d$sex <- factor(d$sex)
  needed <- unique(c(outcome, "mean_z_score", "ssgsea_score", covariates))
  d <- d[complete.cases(d[, needed, drop = FALSE]), needed, drop = FALSE]
  if (nrow(d) < 20 || length(unique(d[[outcome]])) != 2) {
    return(data.frame(dataset = dataset, tissue = tissue, ECM_definition = ECM_definition,
                      outcome = outcome, method = c("mean_z", "ssGSEA"),
                      covariates = paste(covariates, collapse = "+"), n = nrow(d),
                      events = if (nrow(d)) sum(d[[outcome]] == 1) else NA_integer_,
                      estimable = FALSE))
  }
  fit_one <- function(score_name, method) {
    f <- reformulate(c(score_name, covariates), response = outcome)
    fit <- glm(f, data = d, family = binomial())
    co <- summary(fit)$coefficients[score_name, ]
    data.frame(dataset = dataset, tissue = tissue, ECM_definition = ECM_definition,
               outcome = outcome, method = method, covariates = paste(covariates, collapse = "+"),
               n = nrow(d), events = sum(d[[outcome]] == 1), estimable = TRUE,
               beta = unname(co["Estimate"]), SE = unname(co["Std. Error"]),
               OR_per_SD = exp(unname(co["Estimate"])),
               CI_low = exp(unname(co["Estimate"]) - 1.96 * unname(co["Std. Error"])),
               CI_high = exp(unname(co["Estimate"]) + 1.96 * unname(co["Std. Error"])),
               p_value = unname(co["Pr(>|z|)"]), AIC = AIC(fit))
  }
  bind_rows(fit_one("mean_z_score", "mean_z"), fit_one("ssgsea_score", "ssGSEA"))
}

compare_fair_models <- function(model_table) {
  model_table |>
    filter(estimable) |>
    dplyr::select(dataset, tissue, ECM_definition, outcome, method, n, events, beta, SE, OR_per_SD, CI_low, CI_high, p_value, AIC) |>
    tidyr::pivot_wider(names_from = method, values_from = c(n, events, beta, SE, OR_per_SD, CI_low, CI_high, p_value, AIC)) |>
    mutate(same_n = n_mean_z == n_ssGSEA,
           same_events = events_mean_z == events_ssGSEA,
           direction_consistent = sign(beta_mean_z) == sign(beta_ssGSEA),
           materially_consistent = pmax(CI_low_mean_z, CI_low_ssGSEA) <= pmin(CI_high_mean_z, CI_high_ssGSEA),
           beta_difference = beta_ssGSEA - beta_mean_z)
}

strict_classification <- function(criteria_table) {
  flags <- setNames(criteria_table$passed, criteria_table$criterion)
  if (all(flags)) {
    classification <- "Main third axis supported"
  } else if (sum(flags) >= 4 && isTRUE(flags["Criterion_2"]) && (isTRUE(flags["Criterion_3"]) || isTRUE(flags["Criterion_5"]))) {
    classification <- "Partial support"
  } else {
    classification <- "Not supported"
  }
  if (!isTRUE(flags["Criterion_3"]) && !isTRUE(flags["Criterion_5"])) {
    classification <- "Not supported for external specificity"
  }
  data.frame(classification = classification, criteria_passed = sum(flags),
             criterion_2_required = flags["Criterion_2"],
             criterion_3_or_5_required = flags["Criterion_3"] || flags["Criterion_5"])
}

log_msg("Started corrected ECM third-axis rerun.")

unit <- run_bootstrap_unit_tests()
unit_lines <- c(
  "# CORRECTED_BOOTSTRAP_UNIT_TEST",
  "",
  "- Previous erroneous expression `adjustedRandIndex(cl[idx], cl[idx])` is not used.",
  paste0("- Identical-cluster ARI: ", signif(unit$same_cluster_ARI, 6)),
  paste0("- Shuffled-cluster ARI: ", signif(unit$shuffled_cluster_ARI, 6)),
  paste0("- Corrected bootstrap ARI variance: ", signif(unit$bootstrap_ARI_variance, 6)),
  paste0("- Bootstrap ARI all equal 1: ", unit$bootstrap_all_one),
  paste0("- Unit test passed: ", unit$passed)
)
writeLines(unit_lines, file.path(CONFIG$out_dir, "CORRECTED_BOOTSTRAP_UNIT_TEST.md"))
writeLines(unit_lines, file.path(CONFIG$out_dir, "02_corrected_bootstrap", "CORRECTED_BOOTSTRAP_UNIT_TEST.md"))
if (!isTRUE(unit$passed)) stop("Corrected bootstrap unit test failed.")

log_msg("Loading fixed inputs.")
expr152 <- read_expr(CONFIG$expr152)
meta152 <- fread(CONFIG$meta152, data.table = FALSE)
expr118 <- read_expr(CONFIG$expr118)
meta118 <- fread(CONFIG$meta118, data.table = FALSE)
if (!setequal(colnames(expr152), meta152$sample_id)) stop("GSE152004 expression/metadata sample IDs mismatch")
expr152 <- expr152[, meta152$sample_id, drop = FALSE]
if (!setequal(colnames(expr118), meta118$sample_id)) stop("GSE118761 expression/metadata sample IDs mismatch")
expr118 <- expr118[, meta118$sample_id, drop = FALSE]

modules <- read_gmt(CONFIG$modules_gmt)
repair32 <- unique(toupper(trimws(fread(CONFIG$repair32, data.table = FALSE)[[1]])))
external_genes <- fread(CONFIG$gene_sets_external, data.table = FALSE)
gene_sets <- list(
  original32 = repair32,
  Reactome_ECM = unique(toupper(trimws(external_genes$gene_symbol[external_genes$gene_set == "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"]))),
  NABA_core = unique(toupper(trimws(external_genes$gene_symbol[external_genes$gene_set == "NABA_CORE_MATRISOME"])))
)
external_sets <- c("Reactome_ECM", "NABA_core")
t2_genes <- modules[["T2_SIGNATURE_3GENE"]]
ifn_genes <- modules[["IFN_MODULE_CURRENT"]]

random_members <- read.csv(gzfile(CONFIG$random_members), check.names = FALSE)
random_members <- random_members[random_members$ECM_definition %in% external_sets, , drop = FALSE]
random_members$iter <- as.integer(random_members$iter)
random_counts <- as.data.frame(table(random_members$ECM_definition))
names(random_counts) <- c("ECM_definition", "n_iterations")
if (!all(setNames(random_counts$n_iterations, random_counts$ECM_definition)[external_sets] == CONFIG$permutation_n)) {
  stop("Random gene-set membership does not contain exactly 10,000 Reactome and 10,000 NABA iterations.")
}

input_files <- data.frame(
  role = names(CONFIG)[vapply(CONFIG, is.character, logical(1))],
  path = unlist(CONFIG[vapply(CONFIG, is.character, logical(1))], use.names = FALSE),
  stringsAsFactors = FALSE
) |>
  filter(file.exists(path)) |>
  mutate(size_bytes = file.info(path)$size, md5 = vapply(path, md5_or_na, character(1)))

coverage_rows <- list()
for (cohort in c("GSE152004", "GSE118761")) {
  expr <- if (cohort == "GSE152004") expr152 else expr118
  for (nm in names(gene_sets)) {
    present <- intersect(gene_sets[[nm]], rownames(expr))
    coverage_rows[[length(coverage_rows) + 1]] <- data.frame(
      cohort = cohort, ECM_definition = nm, fixed_members = length(gene_sets[[nm]]),
      present = length(present), coverage = length(present) / length(gene_sets[[nm]])
    )
  }
}
coverage_df <- bind_rows(coverage_rows)
audit_lines <- c(
  "# CORRECTION_INPUT_AUDIT",
  "",
  paste0("- GSE152004 expression: ", CONFIG$expr152, "; samples=", ncol(expr152), "; genes=", nrow(expr152), "."),
  paste0("- GSE152004 metadata columns: ", paste(names(meta152), collapse = ", ")),
  paste0("- GSE118761 expression: ", CONFIG$expr118, "; samples=", ncol(expr118), "; genes=", nrow(expr118), "."),
  paste0("- GSE118761 metadata columns: ", paste(names(meta118), collapse = ", ")),
  paste0("- Random membership file: ", CONFIG$random_members, "; Reactome/NABA iterations=", paste(random_counts$ECM_definition, random_counts$n_iterations, sep = "=", collapse = "; "), "."),
  paste0("- Fixed gene-set sizes: original32=", length(gene_sets$original32), "; Reactome=", length(gene_sets$Reactome_ECM), "; NABA=", length(gene_sets$NABA_core), "."),
  "- original32 is retained only as historical/descriptive baseline and is excluded from external empirical FDR.",
  "- No manuscript, Figure, title, or response letter is modified."
)
writeLines(audit_lines, file.path(CONFIG$out_dir, "CORRECTION_INPUT_AUDIT.md"))
writeLines(audit_lines, file.path(CONFIG$out_dir, "01_input_audit", "CORRECTION_INPUT_AUDIT.md"))
write_wb(file.path(CONFIG$out_dir, "01_input_audit", "input_audit_tables.xlsx"), list(files = input_files, coverage = coverage_df, random_counts = random_counts))

base_scores <- function(expr, meta) {
  keep <- intersect(names(meta), c("sample_id", "subject_id", "tissue", "asthma", "atopy", "wheeze", "T2_high", "age", "sex", "batch"))
  out <- meta[, keep, drop = FALSE]
  out$T2 <- score_mean_z(expr, t2_genes)$score[out$sample_id]
  out$IFN <- score_mean_z(expr, ifn_genes)$score[out$sample_id]
  for (nm in names(gene_sets)) out[[nm]] <- score_mean_z(expr, gene_sets[[nm]])$score[out$sample_id]
  out$T2_z <- zscore(out$T2)
  out$IFN_z <- zscore(out$IFN)
  for (nm in names(gene_sets)) out[[paste0(nm, "_z")]] <- zscore(out[[nm]])
  out
}
scores152 <- base_scores(expr152, meta152)
scores152$tissue <- "nasal"
scores118 <- base_scores(expr118, meta118)
scores118$wheeze_duplicate_asthma <- if ("wheeze" %in% names(scores118)) scores118$wheeze == scores118$asthma else NA

log_msg("Fitting full-data clusters and corrected real bootstrap.")
model_x <- list(
  two_axis = scores152[, c("T2_z", "IFN_z")],
  original32 = setNames(scores152[, c("T2_z", "IFN_z", "original32_z")], c("T2_z", "IFN_z", "ECM_z")),
  Reactome_ECM = setNames(scores152[, c("T2_z", "IFN_z", "Reactome_ECM_z")], c("T2_z", "IFN_z", "ECM_z")),
  NABA_core = setNames(scores152[, c("T2_z", "IFN_z", "NABA_core_z")], c("T2_z", "IFN_z", "ECM_z"))
)
fits <- lapply(model_x, kmeans_full)
for (nm in names(fits)) scores152[[paste0(nm, "_cluster")]] <- fits[[nm]]$cluster
full_metrics <- bind_rows(lapply(names(model_x), function(nm) {
  cluster_metrics(model_x[[nm]], fits[[nm]]$cluster) |> mutate(model = nm, .before = 1)
}))

boot_idx_100 <- make_bootstrap_indices(nrow(scores152), CONFIG$permutation_boot_reps)
saveRDS(boot_idx_100, file.path(CONFIG$out_dir, "02_corrected_bootstrap", "bootstrap_indices.rds"))
real_boot_100_detail <- list()
real_boot_100_summary <- list()
for (nm in names(model_x)) {
  detail <- corrected_bootstrap_refit(model_x[[nm]], fits[[nm]]$cluster, fits[[nm]]$centers, boot_idx_100)
  detail$model <- nm
  real_boot_100_detail[[nm]] <- detail
  real_boot_100_summary[[nm]] <- summarise_bootstrap(detail, nm)
}
real_boot_100_summary_df <- bind_rows(real_boot_100_summary)
real_boot_100_detail_df <- bind_rows(real_boot_100_detail)

boot_idx_1000 <- make_bootstrap_indices(nrow(scores152), CONFIG$descriptive_boot_reps, seed = CONFIG$seed + 900000L)
real_boot_1000_summary <- list()
for (nm in names(model_x)) {
  detail <- corrected_bootstrap_refit(model_x[[nm]], fits[[nm]]$cluster, fits[[nm]]$centers, boot_idx_1000)
  real_boot_1000_summary[[nm]] <- summarise_bootstrap(detail, nm)
}
real_boot_1000_summary_df <- bind_rows(real_boot_1000_summary)

two_boot_mean <- real_boot_100_summary_df$mean_ARI[real_boot_100_summary_df$model == "two_axis"]
real_delta <- full_metrics |>
  filter(model != "two_axis") |>
  transmute(ECM_definition = model,
            delta_silhouette = mean_silhouette - full_metrics$mean_silhouette[full_metrics$model == "two_axis"]) |>
  left_join(real_boot_100_summary_df |> filter(model != "two_axis") |> transmute(ECM_definition = model, bootstrap_mean_ARI = mean_ARI), by = "ECM_definition") |>
  mutate(delta_bootstrap_ARI = bootstrap_mean_ARI - two_boot_mean)

write_wb(file.path(CONFIG$out_dir, "corrected_real_bootstrap_results.xlsx"),
         list(summary_100 = real_boot_100_summary_df, detail_100 = real_boot_100_detail_df, full_metrics = full_metrics, real_delta = real_delta))
write_wb(file.path(CONFIG$out_dir, "02_corrected_bootstrap", "real_models_bootstrap_100.xlsx"),
         list(summary_100 = real_boot_100_summary_df, detail_100 = real_boot_100_detail_df, full_metrics = full_metrics, real_delta = real_delta))
write_wb(file.path(CONFIG$out_dir, "02_corrected_bootstrap", "real_models_bootstrap_1000_descriptive.xlsx"),
         list(summary_1000 = real_boot_1000_summary_df))

log_msg("Running residual nonredundancy audit.")
residual_rows <- list()
for (nm in external_sets) {
  fit <- lm(as.formula(paste0(nm, "_z ~ T2_z + IFN_z")), data = scores152)
  resid_z <- zscore(resid(fit))
  x_res <- data.frame(T2_z = scores152$T2_z, IFN_z = scores152$IFN_z, ECM_resid_z = resid_z)
  kmr <- kmeans_full(x_res)
  residual_rows[[nm]] <- data.frame(
    ECM_definition = nm,
    R2_T2_IFN = summary(fit)$r.squared,
    residual_cluster_ARI_vs_original = adjustedRandIndex(kmr$cluster, fits[[nm]]$cluster),
    residual_cluster_NMI_vs_original = normalized_mutual_information(kmr$cluster, fits[[nm]]$cluster),
    residual_ECM_eta2 = eta2_oneway(resid_z, kmr$cluster)
  )
}
residual_df <- bind_rows(residual_rows)

pair_ari <- adjustedRandIndex(fits$Reactome_ECM$cluster, fits$NABA_core$cluster)
pair_nmi <- normalized_mutual_information(fits$Reactome_ECM$cluster, fits$NABA_core$cluster)
matched_naba <- best_match(fits$Reactome_ECM$cluster, fits$NABA_core$cluster)
pair_agreement <- mean(as.character(fits$Reactome_ECM$cluster) == matched_naba)
pairwise_external <- data.frame(definition_a = "Reactome_ECM", definition_b = "NABA_core",
                                ARI = pair_ari, NMI = pair_nmi, Hungarian_agreement = pair_agreement,
                                changed_count_after_Hungarian = sum(as.character(fits$Reactome_ECM$cluster) != matched_naba))

log_msg("Official GSVA ssGSEA and fair model comparison.")
unit_ss <- GSVA::gsva(matrix(rnorm(200), nrow = 20, dimnames = list(paste0("G", 1:20), paste0("S", 1:10))),
                     list(test = paste0("G", 1:8)), method = "ssgsea", min.sz = 1, max.sz = Inf, ssgsea.norm = TRUE, verbose = FALSE)
if (!all(dim(unit_ss) == c(1, 10))) stop("Official GSVA ssGSEA unit test failed.")

ss_exprs <- list(GSE152004 = expr152, GSE118761 = expr118)
ss_scores <- list()
for (cohort in names(ss_exprs)) {
  expr <- ss_exprs[[cohort]]
  gs <- lapply(gene_sets, intersect, y = rownames(expr))
  ans <- GSVA::gsva(expr, gs, method = "ssgsea", min.sz = 3, max.sz = Inf, ssgsea.norm = TRUE, verbose = FALSE)
  for (nm in rownames(ans)) {
    ss_scores[[paste(cohort, nm, sep = "::")]] <- setNames(as.numeric(ans[nm, ]), colnames(ans))
  }
}

fair_models <- list()
sample_audit <- list()
for (nm in names(gene_sets)) {
  sample_audit[[length(sample_audit) + 1]] <- data.frame(dataset = "GSE152004", tissue = "nasal", ECM_definition = nm,
                                                         outcome = "T2_high", model = "T2_high ~ ECM + IFN_z",
                                                         n_complete = sum(complete.cases(data.frame(scores152$T2_high, scores152[[paste0(nm, "_z")]], ss_scores[[paste("GSE152004", nm, sep = "::")]][scores152$sample_id], scores152$IFN_z))),
                                                         events = sum(scores152$T2_high == 1, na.rm = TRUE))
  fair_models[[length(fair_models) + 1]] <- fit_fair_logistic_pair(scores152, "T2_high",
                                                                   scores152[[paste0(nm, "_z")]],
                                                                   ss_scores[[paste("GSE152004", nm, sep = "::")]][scores152$sample_id],
                                                                   c("IFN_z"), "GSE152004", "nasal", nm)
  sample_audit[[length(sample_audit) + 1]] <- data.frame(dataset = "GSE152004", tissue = "nasal", ECM_definition = nm,
                                                         outcome = "asthma", model = "asthma ~ ECM + IFN_z",
                                                         n_complete = sum(complete.cases(data.frame(scores152$asthma, scores152[[paste0(nm, "_z")]], ss_scores[[paste("GSE152004", nm, sep = "::")]][scores152$sample_id], scores152$IFN_z))),
                                                         events = sum(scores152$asthma == 1, na.rm = TRUE))
  fair_models[[length(fair_models) + 1]] <- fit_fair_logistic_pair(scores152, "asthma",
                                                                   scores152[[paste0(nm, "_z")]],
                                                                   ss_scores[[paste("GSE152004", nm, sep = "::")]][scores152$sample_id],
                                                                   c("IFN_z"), "GSE152004", "nasal", nm)
  for (tiss in c("nasal", "tracheal")) {
    d <- scores118 |> filter(tissue == tiss)
    for (outcome in c("atopy", "asthma")) {
      sample_audit[[length(sample_audit) + 1]] <- data.frame(dataset = "GSE118761", tissue = tiss, ECM_definition = nm,
                                                             outcome = outcome, model = paste0(outcome, " ~ ECM + IFN_z + age + sex"),
                                                             n_complete = sum(complete.cases(data.frame(d[[outcome]], d[[paste0(nm, "_z")]], ss_scores[[paste("GSE118761", nm, sep = "::")]][d$sample_id], d$IFN_z, d$age, d$sex))),
                                                             events = sum(d[[outcome]] == 1, na.rm = TRUE))
      fair_models[[length(fair_models) + 1]] <- fit_fair_logistic_pair(d, outcome,
                                                                       d[[paste0(nm, "_z")]],
                                                                       ss_scores[[paste("GSE118761", nm, sep = "::")]][d$sample_id],
                                                                       c("IFN_z", "age", "sex"), "GSE118761", tiss, nm)
    }
  }
}
fair_models_df <- bind_rows(fair_models)
fair_compare <- compare_fair_models(fair_models_df)
ss_cor <- bind_rows(lapply(names(gene_sets), function(nm) {
  data.frame(ECM_definition = nm,
             GSE152004_pearson = cor(scores152[[paste0(nm, "_z")]], zscore(ss_scores[[paste("GSE152004", nm, sep = "::")]][scores152$sample_id]), use = "pairwise.complete.obs"),
             GSE152004_spearman = cor(scores152[[paste0(nm, "_z")]], zscore(ss_scores[[paste("GSE152004", nm, sep = "::")]][scores152$sample_id]), method = "spearman", use = "pairwise.complete.obs"),
             GSE118761_pearson = cor(scores118[[paste0(nm, "_z")]], zscore(ss_scores[[paste("GSE118761", nm, sep = "::")]][scores118$sample_id]), use = "pairwise.complete.obs"),
             GSE118761_spearman = cor(scores118[[paste0(nm, "_z")]], zscore(ss_scores[[paste("GSE118761", nm, sep = "::")]][scores118$sample_id]), method = "spearman", use = "pairwise.complete.obs"))
}))
sample_audit_df <- bind_rows(sample_audit)
write_wb(file.path(CONFIG$out_dir, "fair_mean_z_vs_ssGSEA_models.xlsx"),
         list(models = fair_models_df, comparison = fair_compare, score_correlations = ss_cor))
write_wb(file.path(CONFIG$out_dir, "04_fair_mean_z_ssGSEA_models", "fair_mean_z_vs_ssGSEA_models.xlsx"),
         list(models = fair_models_df, comparison = fair_compare, score_correlations = ss_cor))
write_wb(file.path(CONFIG$out_dir, "04_fair_mean_z_ssGSEA_models", "fair_mean_z_vs_ssGSEA_sample_audit.xlsx"),
         list(sample_audit = sample_audit_df))

projection_real_rows <- list()
for (nm in external_sets) {
  train_x <- setNames(scores152[, c("T2_z", "IFN_z", paste0(nm, "_z"))], c("T2_z", "IFN_z", "ECM_z"))
  for (tiss in c("nasal", "tracheal")) {
    valid <- scores118 |> filter(tissue == tiss)
    valid_x <- setNames(valid[, c("T2_z", "IFN_z", paste0(nm, "_z"))], c("T2_z", "IFN_z", "ECM_z"))
    pr <- projection_margin(train_x, fits[[nm]]$cluster, valid_x)
    projection_real_rows[[length(projection_real_rows) + 1]] <- data.frame(
      ECM_definition = nm, tissue = tiss, sample_id = valid$sample_id, pr
    )
  }
}
projection_real <- bind_rows(projection_real_rows)
projection_real_summary <- projection_real |>
  group_by(ECM_definition, tissue) |>
  summarise(mean_assignment_margin = mean(assignment_margin), low_confidence_rate = mean(low_confidence), .groups = "drop")

log_msg("Running corrected random external specificity permutation.")
random_stats_file <- file.path(CONFIG$out_dir, "corrected_random_bootstrap_statistics.csv.gz")
if (file.exists(random_stats_file)) unlink(random_stats_file)
collapse_audit <- list()

process_random_row <- function(set_name, iter, genes_string) {
  genes <- strsplit(genes_string, ";", fixed = TRUE)[[1]]
  sc152 <- score_mean_z(expr152, genes)
  sc118 <- score_mean_z(expr118, genes)
  if (!sc152$estimable || !sc118$estimable) {
    return(list(stat = data.frame(ECM_definition = set_name, iter = iter, valid = FALSE,
                                  reason = "score_not_estimable", n_genes_152 = sc152$n_used, n_genes_118 = sc118$n_used),
                boot = NULL))
  }
  rnd_z <- zscore(sc152$score[scores152$sample_id])
  rnd118_z <- zscore(sc118$score[scores118$sample_id])
  names(rnd_z) <- scores152$sample_id
  names(rnd118_z) <- scores118$sample_id
  x <- data.frame(T2_z = scores152$T2_z, IFN_z = scores152$IFN_z, ECM_z = rnd_z)
  if (!all(complete.cases(x))) {
    return(list(stat = data.frame(ECM_definition = set_name, iter = iter, valid = FALSE,
                                  reason = "missing_train_feature", n_genes_152 = sc152$n_used, n_genes_118 = sc118$n_used),
                boot = NULL))
  }
  km <- kmeans_full(x)
  met <- cluster_metrics(x, km$cluster)
  boot <- corrected_bootstrap_refit(x, km$cluster, km$centers, boot_idx_100)
  boot_sum <- summarise_bootstrap(boot, paste0(set_name, "_random_", iter))
  fixed_other <- if (set_name == "Reactome_ECM") fits$NABA_core$cluster else fits$Reactome_ECM$cluster
  valid_nasal <- scores118 |> filter(tissue == "nasal")
  valid_x <- data.frame(T2_z = valid_nasal$T2_z, IFN_z = valid_nasal$IFN_z,
                        ECM_z = rnd118_z[valid_nasal$sample_id])
  pr <- projection_margin(x, km$cluster, valid_x)
  stat <- data.frame(
    ECM_definition = set_name, iter = iter, valid = TRUE, reason = NA_character_,
    n_genes_152 = sc152$n_used, n_genes_118 = sc118$n_used,
    coverage_152 = sc152$coverage, coverage_118 = sc118$coverage,
    mean_silhouette = met$mean_silhouette,
    delta_silhouette = met$mean_silhouette - full_metrics$mean_silhouette[full_metrics$model == "two_axis"],
    bootstrap_mean_ARI = boot_sum$mean_ARI,
    bootstrap_sd_ARI = boot_sum$sd_ARI,
    bootstrap_reps_valid = boot_sum$reps_valid,
    bootstrap_collapse_rate = boot_sum$collapse_rate,
    delta_bootstrap_ARI = boot_sum$mean_ARI - two_boot_mean,
    cross_definition_cluster_ARI = adjustedRandIndex(km$cluster, fixed_other),
    cross_definition_cluster_NMI = normalized_mutual_information(km$cluster, fixed_other),
    nasal_projection_assignment_margin = mean(pr$assignment_margin),
    nasal_projection_low_confidence_rate = mean(pr$low_confidence),
    calinski_harabasz = met$calinski_harabasz,
    davies_bouldin = met$davies_bouldin,
    centroid_separation = met$centroid_min_distance,
    min_cluster_size = met$min_cluster_size
  )
  list(stat = stat, boot = data.frame(ECM_definition = set_name, iter = iter, boot))
}

all_random_stats <- list()
for (set_name in external_sets) {
  log_msg("Corrected permutation start for ", set_name)
  rows <- random_members |> filter(ECM_definition == set_name) |> arrange(iter)
  stats_list <- list()
  boot_collapse_list <- list()
  for (ii in seq_len(nrow(rows))) {
    res <- process_random_row(rows$ECM_definition[ii], rows$iter[ii], rows$genes[ii])
    stats_list[[length(stats_list) + 1]] <- res$stat
    if (!is.null(res$boot)) boot_collapse_list[[length(boot_collapse_list) + 1]] <- res$boot
    if (ii %% CONFIG$checkpoint_every == 0) {
      chk_stats <- bind_rows(stats_list)
      chk_path <- file.path(CONFIG$out_dir, "03_external_specificity_permutation", "checkpoints", paste0(set_name, "_corrected_stats_", ii, ".csv.gz"))
      fwrite(chk_stats, chk_path)
      if (length(boot_collapse_list)) {
        chk_boot <- bind_rows(boot_collapse_list)
        fwrite(chk_boot, file.path(CONFIG$out_dir, "03_external_specificity_permutation", "checkpoints", paste0(set_name, "_bootstrap_detail_", ii, ".csv.gz")))
      }
      log_perm("Checkpoint ", set_name, " valid=", sum(chk_stats$valid), " iter=", ii)
    }
  }
  set_stats <- bind_rows(stats_list)
  if (sum(set_stats$valid) != CONFIG$permutation_n) stop(set_name, " does not have exactly 10,000 valid iterations.")
  fwrite(set_stats, random_stats_file, append = file.exists(random_stats_file))
  all_random_stats[[set_name]] <- set_stats
  if (length(boot_collapse_list)) collapse_audit[[set_name]] <- bind_rows(boot_collapse_list)
  log_msg("Corrected permutation completed for ", set_name)
}
random_stats <- bind_rows(all_random_stats)
collapse_detail <- bind_rows(collapse_audit)
write_wb(file.path(CONFIG$out_dir, "02_corrected_bootstrap", "bootstrap_collapse_audit.xlsx"),
         list(random_bootstrap = collapse_detail |> group_by(ECM_definition, iter) |> summarise(collapse_rate = mean(collapse), valid_reps = sum(!collapse), .groups = "drop"),
              real_bootstrap = real_boot_100_detail_df |> group_by(model) |> summarise(collapse_rate = mean(collapse), valid_reps = sum(!collapse), .groups = "drop")))

true_metrics <- bind_rows(lapply(external_sets, function(nm) {
  other <- if (nm == "Reactome_ECM") "NABA_core" else "Reactome_ECM"
  data.frame(
    ECM_definition = nm,
    delta_bootstrap_ARI = real_delta$delta_bootstrap_ARI[real_delta$ECM_definition == nm],
    delta_silhouette = real_delta$delta_silhouette[real_delta$ECM_definition == nm],
    cross_definition_cluster_ARI = adjustedRandIndex(fits[[nm]]$cluster, fits[[other]]$cluster),
    nasal_projection_assignment_margin = projection_real_summary$mean_assignment_margin[projection_real_summary$ECM_definition == nm & projection_real_summary$tissue == "nasal"],
    nasal_projection_low_confidence_rate = projection_real_summary$low_confidence_rate[projection_real_summary$ECM_definition == nm & projection_real_summary$tissue == "nasal"]
  )
}))
metric_names <- c("delta_bootstrap_ARI", "delta_silhouette", "cross_definition_cluster_ARI", "nasal_projection_assignment_margin")
emp_rows <- list()
for (nm in external_sets) {
  rnd <- random_stats |> filter(ECM_definition == nm)
  for (metric in metric_names) {
    ep <- empirical_p_larger(rnd[[metric]], true_metrics[[metric]][true_metrics$ECM_definition == nm])
    emp_rows[[length(emp_rows) + 1]] <- cbind(ECM_definition = nm, metric = metric, ep)
  }
}
empirical_df <- bind_rows(emp_rows)
empirical_df$BH_FDR_external8 <- p.adjust(empirical_df$empirical_p, "BH")
write_wb(file.path(CONFIG$out_dir, "corrected_permutation_empirical_pvalues.xlsx"), list(empirical_pvalues = empirical_df, true_metrics = true_metrics))
write_wb(file.path(CONFIG$out_dir, "corrected_external8_FDR.xlsx"), list(external8 = empirical_df))
write_wb(file.path(CONFIG$out_dir, "03_external_specificity_permutation", "corrected_permutation_empirical_pvalues.xlsx"), list(empirical_pvalues = empirical_df, true_metrics = true_metrics))

projection_specificity <- empirical_df |>
  filter(metric == "nasal_projection_assignment_margin") |>
  left_join(projection_real_summary |> filter(tissue == "nasal"), by = "ECM_definition") |>
  left_join(random_stats |> group_by(ECM_definition) |> summarise(random_median_low_confidence = median(nasal_projection_low_confidence_rate), .groups = "drop"), by = "ECM_definition") |>
  mutate(projection_specificity_supported = empirical_p <= 0.05 & BH_FDR_external8 < 0.05 & low_confidence_rate <= random_median_low_confidence)
write_wb(file.path(CONFIG$out_dir, "corrected_projection_specificity.xlsx"),
         list(projection_specificity = projection_specificity, real_projection_assignments = projection_real, real_projection_summary = projection_real_summary))
write_wb(file.path(CONFIG$out_dir, "05_corrected_projection", "corrected_projection_specificity.xlsx"),
         list(projection_specificity = projection_specificity, real_projection_assignments = projection_real, real_projection_summary = projection_real_summary))

pipeline_lines <- c(
  "# PIPELINE_EQUALITY_AUDIT",
  "",
  "- Real Reactome, real NABA, and matched random sets use cohort-internal mean-z scoring.",
  "- Real and random models use fixed k=3 and full-data k-means nstart=500.",
  "- Real and random bootstrap use the same 100 bootstrap index vectors saved in bootstrap_indices.rds.",
  "- Bootstrap refit uses bootstrap samples with replacement, retaining duplicate samples.",
  "- Bootstrap refit initializes from full-data centroids and uses the same Hungarian centroid matching.",
  "- Collapse rule is identical: any missing cluster or minimum predicted cluster size <5 yields collapse=TRUE and ARI/NMI=NA.",
  "- Reactome random cross-definition ARI is compared against fixed NABA cluster.",
  "- NABA random cross-definition ARI is compared against fixed Reactome cluster.",
  "- Projection formula, nearest-centroid margin, and low-confidence cutoff are identical for real and random sets.",
  "- original32 is excluded from external specificity empirical P values and BH_FDR_external8."
)
writeLines(pipeline_lines, file.path(CONFIG$out_dir, "PIPELINE_EQUALITY_AUDIT.md"))
writeLines(pipeline_lines, file.path(CONFIG$out_dir, "03_external_specificity_permutation", "PIPELINE_EQUALITY_AUDIT.md"))

log_msg("Strict reclassification.")
criterion1_by_set <- real_delta |>
  mutate(internal_no_obvious_joint_deterioration = !(delta_bootstrap_ARI < -0.05 & delta_silhouette < -0.05),
         set_pass = (delta_bootstrap_ARI > 0 | delta_silhouette > 0) & internal_no_obvious_joint_deterioration)
criterion1_pass <- sum(criterion1_by_set$set_pass) >= 2
criterion2_by_set <- residual_df |>
  mutate(set_pass = R2_T2_IFN < 0.25 & residual_cluster_ARI_vs_original >= 0.50 & residual_ECM_eta2 > 0.05)
criterion2_pass <- any(criterion2_by_set$set_pass)
criterion3_rows <- empirical_df |>
  filter(metric %in% c("delta_bootstrap_ARI", "delta_silhouette", "cross_definition_cluster_ARI")) |>
  mutate(set_pass = BH_FDR_external8 < 0.05)
criterion3_pass <- any(criterion3_rows$set_pass)
criterion4_pass <- pair_ari >= 0.50 && pair_nmi >= 0.50 && pair_agreement >= 0.80
criterion5_pass <- any(projection_specificity$projection_specificity_supported)
criterion6_base <- fair_compare |>
  filter(ECM_definition %in% external_sets)
criterion6_pass <- mean(criterion6_base$direction_consistent, na.rm = TRUE) >= 0.75

criteria_table <- data.frame(
  criterion = paste0("Criterion_", 1:6),
  description = c("Two-axis to three-axis increment",
                  "Residual nonredundancy",
                  "External specificity from non-original32 empirical tests",
                  "Reactome/NABA external-definition consistency",
                  "Cross-cohort nasal projection specificity",
                  "Fair mean-z versus official ssGSEA model direction consistency"),
  threshold = c("At least two of original32/Reactome/NABA have corrected delta_bootstrap_ARI>0 or delta_silhouette>0 without joint obvious deterioration",
                "Reactome or NABA: R2<0.25, residual ARI>=0.50, residual ECM eta2>0.05",
                "Reactome or NABA: at least one clustering stability/separation metric BH_FDR_external8<0.05",
                "Reactome vs NABA ARI>=0.50, NMI>=0.50, Hungarian agreement>=0.80",
                "Reactome or NABA projection margin empirical P<=0.05, BH_FDR_external8<0.05, low-confidence<=random median",
                ">=75% direction consistency in prespecified Reactome/NABA fair models"),
  passed = c(criterion1_pass, criterion2_pass, criterion3_pass, criterion4_pass, criterion5_pass, criterion6_pass),
  affected_by_original32 = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  source_file = c("corrected_real_bootstrap_results.xlsx",
                  "ECM_third_axis_corrected_rerun_full_results.xlsx",
                  "corrected_external8_FDR.xlsx",
                  "ECM_third_axis_corrected_rerun_full_results.xlsx",
                  "corrected_projection_specificity.xlsx",
                  "fair_mean_z_vs_ssGSEA_models.xlsx")
)
classification <- strict_classification(criteria_table)
write_wb(file.path(CONFIG$out_dir, "STRICT_CLASSIFICATION_AUDIT.xlsx"),
         list(criteria = criteria_table, classification = classification,
              criterion1 = criterion1_by_set, criterion2 = criterion2_by_set,
              criterion3 = criterion3_rows, criterion4 = pairwise_external,
              criterion5 = projection_specificity, criterion6 = criterion6_base))
write_wb(file.path(CONFIG$out_dir, "06_strict_reclassification", "STRICT_CLASSIFICATION_AUDIT.xlsx"),
         list(criteria = criteria_table, classification = classification,
              criterion1 = criterion1_by_set, criterion2 = criterion2_by_set,
              criterion3 = criterion3_rows, criterion4 = pairwise_external,
              criterion5 = projection_specificity, criterion6 = criterion6_base))

full_results_path <- file.path(CONFIG$out_dir, "ECM_third_axis_corrected_rerun_full_results.xlsx")
write_wb(full_results_path, list(
  input_files = input_files,
  gene_coverage = coverage_df,
  full_cluster_metrics = full_metrics,
  real_bootstrap_summary = real_boot_100_summary_df,
  real_delta = real_delta,
  residual_nonredundancy = residual_df,
  external_pairwise = pairwise_external,
  empirical_external8 = empirical_df,
  projection_specificity = projection_specificity,
  fair_models = fair_models_df,
  fair_model_comparison = fair_compare,
  strict_criteria = criteria_table,
  classification = classification
))

manuscript_implication <- if (classification$classification == "Main third axis supported") {
  "retain repair-ECM as equal principal third axis"
} else if (classification$classification == "Partial support") {
  "retain as nonredundant continuous axis but not proven cluster-improving axis"
} else {
  "remove from primary framework or demote; corrected evidence is not sufficient for external specificity"
}
summary_lines <- c(
  "# ECM_third_axis_corrected_rerun_report",
  "",
  "## Fixed Opening Summary",
  "- Previous bootstrap bug removed: Yes",
  paste0("- Random bootstrap ARI variance nonzero: ", all(random_stats$bootstrap_sd_ARI > 0, na.rm = TRUE)),
  "- External sets tested for specificity: Reactome/NABA",
  "- Original32 included in external empirical FDR: No",
  paste0("- Corrected two-axis bootstrap mean ARI: ", signif(two_boot_mean, 5)),
  paste0("- Corrected Reactome bootstrap mean ARI and delta: ", signif(real_delta$bootstrap_mean_ARI[real_delta$ECM_definition == "Reactome_ECM"], 5), " / ", signif(real_delta$delta_bootstrap_ARI[real_delta$ECM_definition == "Reactome_ECM"], 5)),
  paste0("- Corrected NABA bootstrap mean ARI and delta: ", signif(real_delta$bootstrap_mean_ARI[real_delta$ECM_definition == "NABA_core"], 5), " / ", signif(real_delta$delta_bootstrap_ARI[real_delta$ECM_definition == "NABA_core"], 5)),
  paste0("- Reactome/NABA delta silhouette: ", signif(real_delta$delta_silhouette[real_delta$ECM_definition == "Reactome_ECM"], 5), " / ", signif(real_delta$delta_silhouette[real_delta$ECM_definition == "NABA_core"], 5)),
  paste0("- Reactome/NABA cross-definition cluster ARI: ", signif(pair_ari, 5)),
  paste0("- Reactome four empirical P and BH_FDR_external8: ", paste(empirical_df |> filter(ECM_definition == "Reactome_ECM") |> transmute(x = paste0(metric, "=", signif(empirical_p, 4), "/", signif(BH_FDR_external8, 4))) |> pull(x), collapse = "; ")),
  paste0("- NABA four empirical P and BH_FDR_external8: ", paste(empirical_df |> filter(ECM_definition == "NABA_core") |> transmute(x = paste0(metric, "=", signif(empirical_p, 4), "/", signif(BH_FDR_external8, 4))) |> pull(x), collapse = "; ")),
  paste0("- Reactome/NABA nasal projection percentile and FDR: ",
         paste(projection_specificity |> transmute(x = paste0(ECM_definition, "=", signif(true_percentile, 4), "/", signif(BH_FDR_external8, 4))) |> pull(x), collapse = "; ")),
  paste0("- Fair mean-z versus ssGSEA direction consistency: ", signif(mean(criterion6_base$direction_consistent, na.rm = TRUE), 4)),
  paste0("- Criterion 1-6 results: ", paste(criteria_table$criterion, criteria_table$passed, sep = "=", collapse = "; ")),
  paste0("- Strict final classification: ", classification$classification),
  paste0("- Manuscript implication: ", manuscript_implication),
  "",
  "## Notes",
  "- This corrected rerun writes to a new directory and does not modify prior result files.",
  "- original32 is descriptive only and does not enter the external8 FDR family.",
  "- Official GSVA was used for ssGSEA; no custom fallback was used.",
  "- Reactome and NABA each have exactly 10,000 valid corrected random iterations.",
  "- No manuscript, Figure, title, or response letter was modified."
)
writeLines(summary_lines, file.path(CONFIG$out_dir, "ECM_third_axis_corrected_rerun_report.md"))
writeLines(summary_lines, file.path(CONFIG$out_dir, "07_report", "ECM_third_axis_corrected_rerun_report.md"))
html <- paste0("<!doctype html><html><head><meta charset='utf-8'><title>ECM corrected rerun</title>",
               "<style>body{font-family:Arial,sans-serif;max-width:1100px;margin:32px auto;line-height:1.45} pre{white-space:pre-wrap;background:#f7f7f7;border:1px solid #ddd;padding:12px}</style>",
               "</head><body><h1>ECM third-axis corrected rerun</h1><pre>",
               htmltools::htmlEscape(paste(summary_lines, collapse = "\n")), "</pre></body></html>")
writeLines(html, file.path(CONFIG$out_dir, "ECM_third_axis_corrected_rerun_report.html"))
writeLines(html, file.path(CONFIG$out_dir, "07_report", "ECM_third_axis_corrected_rerun_report.html"))

changelog <- c(
  "# CHANGELOG_corrected_rerun",
  "",
  "- Removed erroneous assignment-bootstrap logic; no `adjustedRandIndex(cl[idx], cl[idx])` is used.",
  "- Implemented true bootstrap refit with replacement, duplicate samples retained, full-data centroid initialization, Hungarian centroid matching, whole-sample reassignment, and ARI/NMI versus reference cluster.",
  "- Used one shared set of 100 bootstrap indices for real Reactome, real NABA, and all matched random Reactome/NABA sets.",
  "- Reused prior locked random_gene_set_membership.csv.gz; no random gene sets were redrawn.",
  "- Excluded original32 from external empirical P values and BH_FDR_external8.",
  "- Replaced original32 self-reference with Reactome/NABA cross-definition cluster ARI.",
  "- Applied BH_FDR_external8 only to Reactome/NABA x four primary metrics.",
  "- Refit fair mean-z and official GSVA ssGSEA models using identical samples and covariates.",
  "- Kept k=3, seed=20260626, and full-data k-means nstart=500.",
  "- No manuscript, Figure, title, or response letter was modified."
)
writeLines(changelog, file.path(CONFIG$out_dir, "CHANGELOG_corrected_rerun.md"))

writeLines(capture.output(sessionInfo()), file.path(CONFIG$out_dir, "sessionInfo.txt"))
file.copy("04_scripts/ECM_third_axis_corrected_rerun.R", file.path(CONFIG$out_dir, "08_code", "ECM_third_axis_corrected_rerun.R"), overwrite = TRUE)

log_msg("Completed corrected ECM third-axis rerun.")

