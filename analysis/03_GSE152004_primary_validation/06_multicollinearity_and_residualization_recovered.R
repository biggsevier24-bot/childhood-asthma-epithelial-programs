options(stringsAsFactors = FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "../.."), winslash = "/", mustWork = TRUE)

vif <- function(x) {
  setNames(vapply(colnames(x), function(v) {
    others <- setdiff(colnames(x), v)
    1 / (1 - summary(lm(reformulate(others, response = v), data = as.data.frame(x)))$r.squared)
  }, numeric(1)), colnames(x))
}
diagnose <- function(d, cohort, tissue, t2, ifn, repair) {
  x <- data.frame(T2 = as.numeric(scale(d[[t2]])), IFN = as.numeric(scale(d[[ifn]])), repair_ECM = as.numeric(scale(d[[repair]])))
  fit <- lm(repair_ECM ~ T2 + IFN, data = x)
  resid <- residuals(fit)
  pearson <- cor(x, method = "pearson", use = "complete.obs")
  spearman <- cor(x, method = "spearman", use = "complete.obs")
  vv <- vif(x)
  data.frame(cohort = cohort, tissue = tissue, n = nrow(x), score_branch = "VST_gene_z_mean",
             Pearson_T2_IFN = pearson["T2", "IFN"], Pearson_T2_repair = pearson["T2", "repair_ECM"],
             Pearson_IFN_repair = pearson["IFN", "repair_ECM"], Spearman_T2_IFN = spearman["T2", "IFN"],
             Spearman_T2_repair = spearman["T2", "repair_ECM"], Spearman_IFN_repair = spearman["IFN", "repair_ECM"],
             VIF_T2 = vv["T2"], VIF_IFN = vv["IFN"], VIF_repair_ECM = vv["repair_ECM"],
             condition_number = kappa(as.matrix(x), exact = TRUE), R2_repair_ECM_on_T2_IFN = summary(fit)$r.squared,
             residual_variance_retained = 1 - summary(fit)$r.squared,
             residual_cor_T2 = cor(resid, x$T2), residual_cor_IFN = cor(resid, x$IFN), stringsAsFactors = FALSE)
}

d152 <- read.csv(file.path(repo, "outputs/03_GSE152004/GSE152004_module_scores_recomputed.csv"), check.names = FALSE)
rows <- list(diagnose(d152, "GSE152004", "nasal", "T2_z", "IFN_z", "repair_ECM_z"))
d118 <- read.csv(file.path(repo, "outputs/02_GSE118761/GSE118761_recomputed_module_scores.csv"), check.names = FALSE)
for (tissue in c("nasal", "tracheal")) {
  z <- d118[tolower(d118$tissue) == tissue, ]
  rows[[length(rows) + 1L]] <- diagnose(z, "GSE118761", tissue, "T2_z", "IFN_z", "repair_ECM_z")
}
out <- do.call(rbind, rows)

fit_r <- lm(repair_ECM_z ~ T2_z + IFN_z, data = d152)
d152$repair_ECM_resid <- as.numeric(scale(residuals(fit_r)))
fit_y <- glm(T2_high ~ repair_ECM_resid + IFN_z, data = d152, family = binomial())
co <- coef(summary(fit_y))["repair_ECM_resid", ]
ci <- confint.default(fit_y, "repair_ECM_resid")
out$residual_repair_ECM_association_outcome <- NA_character_
out$OR <- out$CI_low <- out$CI_high <- out$P <- NA_real_
out$residual_repair_ECM_association_outcome[1] <- "molecular_T2_high"
out$OR[1] <- exp(co["Estimate"]); out$CI_low[1] <- exp(ci[1]); out$CI_high[1] <- exp(ci[2]); out$P[1] <- co["Pr(>|z|)"]
write.csv(out, file.path(repo, "outputs/03_GSE152004/Supplementary_Table_S13_multicollinearity_residualization.csv"), row.names = FALSE)
write.csv(d152[, c("sample_id", "repair_ECM_resid")], file.path(repo, "outputs/03_GSE152004/GSE152004_repair_ECM_residual_scores.csv"), row.names = FALSE)
