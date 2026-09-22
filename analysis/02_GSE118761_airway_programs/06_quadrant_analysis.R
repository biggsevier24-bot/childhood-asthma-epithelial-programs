options(stringsAsFactors = FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "../.."), winslash = "/", mustWork = TRUE)
d <- read.csv(file.path(repo, "outputs/02_GSE118761/GSE118761_recomputed_module_scores.csv"), check.names = FALSE)
d <- d[tolower(d$tissue) == "nasal", ]
d$IFN_high <- d$IFN_z >= median(d$IFN_z, na.rm = TRUE)
d$ECM_high <- d$repair_ECM_z >= median(d$repair_ECM_z, na.rm = TRUE)
d$quadrant <- paste0(ifelse(d$IFN_high, "IFN_high", "IFN_low"), "_", ifelse(d$ECM_high, "ECM_high", "ECM_low"))
write.csv(d, file.path(repo, "outputs/02_GSE118761/GSE118761_nasal_quadrants_recomputed.csv"), row.names = FALSE)
write.csv(d, file.path(repo, "outputs/02_GSE118761/GSE118761_nasal_module_scores_quadrants.csv"), row.names = FALSE)

summary_rows <- do.call(rbind, lapply(sort(unique(d$quadrant)), function(q) {
  x <- d[d$quadrant == q, ]
  data.frame(quadrant = q, n = nrow(x), atopy_rate = mean(x$atopy == 1), asthma_rate = mean(x$asthma == 1),
             wheeze_rate = mean(x$wheeze == 1), mean_repair = mean(x$repair_ECM_z), mean_ifn = mean(x$IFN_z))
}))
write.csv(summary_rows, file.path(repo, "outputs/02_GSE118761/GSE118761_quadrant_summary.csv"), row.names = FALSE)

overlap <- data.frame(nasal_n = nrow(d), identical_asthma_wheeze = identical(as.integer(d$asthma), as.integer(d$wheeze)),
                      discordant_count = sum(d$asthma != d$wheeze, na.rm = TRUE))
write.csv(overlap, file.path(repo, "outputs/02_GSE118761/GSE118761_asthma_wheeze_overlap_audit.csv"), row.names = FALSE)

calc_auc <- function(y, pred) {
  ok <- !is.na(y) & is.finite(pred)
  if (length(unique(y[ok])) < 2L) return(NA_real_)
  if (!requireNamespace("pROC", quietly = TRUE)) stop("Package pROC is required")
  as.numeric(pROC::auc(pROC::roc(y[ok], pred[ok], quiet = TRUE)))
}
stratified_folds <- function(y, k = 5L, seed = 42L) {
  set.seed(seed)
  folds <- integer(length(y))
  for (cls in unique(y)) {
    idx <- which(y == cls)
    folds[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  folds
}

roc_data <- d[d$wheeze == 1 & !is.na(d$atopy) & is.finite(d$IFN_z) & is.finite(d$repair_ECM_z), ]
roc_data$fold <- stratified_folds(roc_data$atopy, 5L, 42L)
roc_data$pred_cv <- NA_real_
fold_rows <- list()
for (fold in sort(unique(roc_data$fold))) {
  train <- roc_data[roc_data$fold != fold, ]
  test <- roc_data[roc_data$fold == fold, ]
  fit <- glm(atopy ~ scale(IFN_z) * scale(repair_ECM_z), data = train, family = binomial())
  pred <- predict(fit, newdata = test, type = "response")
  roc_data$pred_cv[roc_data$fold == fold] <- pred
  fold_rows[[length(fold_rows) + 1L]] <- data.frame(fold = fold, n = nrow(test), events = sum(test$atopy == 1), AUC = calc_auc(test$atopy, pred))
}
fold_auc <- do.call(rbind, fold_rows)
finite_n <- sum(is.finite(fold_auc$AUC))
se <- sd(fold_auc$AUC, na.rm = TRUE) / sqrt(finite_n)
roc_summary <- data.frame(mean_auc = mean(fold_auc$AUC, na.rm = TRUE), sd_auc = sd(fold_auc$AUC, na.rm = TRUE),
                          ci_low = mean(fold_auc$AUC, na.rm = TRUE) - qt(0.975, finite_n - 1) * se,
                          ci_high = mean(fold_auc$AUC, na.rm = TRUE) + qt(0.975, finite_n - 1) * se,
                          pooled_cv_auc = calc_auc(roc_data$atopy, roc_data$pred_cv))
write.csv(fold_auc, file.path(repo, "outputs/02_GSE118761/GSE118761_wheezer_atopy_CV_fold_AUC.csv"), row.names = FALSE)
write.csv(roc_summary, file.path(repo, "outputs/02_GSE118761/GSE118761_wheezer_atopy_CV_ROC_summary.csv"), row.names = FALSE)
