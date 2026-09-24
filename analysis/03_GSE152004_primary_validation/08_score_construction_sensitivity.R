options(stringsAsFactors = FALSE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
source(file.path(repo, "R/common/module_scoring.R"))
raw <- file.path(repo, "data/public/GSE152004/raw/GSE152004_695_raw_counts.txt.gz")
if (!file.exists(raw)) stop("Official GSE152004 raw count file is required")
x <- read.delim(gzfile(raw), check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
m <- as.matrix(x); storage.mode(m) <- "numeric"; rownames(m) <- toupper(rownames(m))
if (anyDuplicated(rownames(m))) m <- rowsum(m, rownames(m), reorder = FALSE)
cpm <- t(t(m) / colSums(m)) * 1e6
log_cpm <- log2(cpm + 1)
d <- read.csv(file.path(repo, "outputs/03_GSE152004/GSE152004_module_scores_recomputed.csv"), check.names = FALSE)
if (anyDuplicated(d$sample_id)) stop("Duplicate sample IDs in primary score table")
log_cpm <- log_cpm[, d$sample_id, drop = FALSE]
genes <- read_gene_set(file.path(repo, "data/reference_gene_sets/repair_ECM_fixed32.txt"))
present <- intersect(unique(toupper(genes)), rownames(log_cpm))
if (length(genes) != 32L || length(present) != 32L) stop("Fixed32 coverage is incomplete")
repair_raw <- colMeans(log_cpm[present, , drop = FALSE], na.rm = TRUE)
repair <- as.numeric(scale(repair_raw))
dat <- data.frame(T2_high = d$T2_high, repair_alt_z = repair, IFN_z = d$IFN_z)
if (anyNA(dat)) stop("Missing outcome or predictor values in canonical model")
fit <- glm(T2_high ~ repair_alt_z + IFN_z, data = dat, family = binomial())
sm <- coef(summary(fit)); ci <- confint.default(fit)
out <- data.frame(
  term = "alternative_repair_per_SD",
  branch = "canonical_public_rebuild",
  N = nobs(fit),
  events = sum(dat$T2_high),
  repair_genes_expected = length(genes),
  repair_genes_used = length(present),
  repair_score_SD = sd(repair_raw),
  beta = unname(sm["repair_alt_z", "Estimate"]),
  SE = unname(sm["repair_alt_z", "Std. Error"]),
  OR = exp(unname(sm["repair_alt_z", "Estimate"])),
  CI_low = exp(ci["repair_alt_z", 1]),
  CI_high = exp(ci["repair_alt_z", 2]),
  p_value = unname(sm["repair_alt_z", "Pr(>|z|)"]),
  IFN_beta = unname(sm["IFN_z", "Estimate"]),
  IFN_SE = unname(sm["IFN_z", "Std. Error"]),
  IFN_OR = exp(unname(sm["IFN_z", "Estimate"])),
  IFN_CI_low = exp(ci["IFN_z", 1]),
  IFN_CI_high = exp(ci["IFN_z", 2]),
  IFN_p_value = unname(sm["IFN_z", "Pr(>|z|)"]),
  score_definition = "mean log2(CPM+1) over fixed32; no gene-wise standardization; predictor standardized per 1 SD",
  adjustment = "primary VST-derived IFN_z",
  stringsAsFactors = FALSE
)
write.csv(out, file.path(repo, "outputs/03_GSE152004/GSE152004_alternative_repair_recomputed.csv"), row.names = FALSE)
