# Recovered source note:
# This file and its paired RNA-seq/gene-set workflow script both derive from
# the same recovered comprehensive source, 04_scripts/External_ECM_multicohort_analysis.R.
# They are intentionally duplicated for public workflow mapping; they should
# not be interpreted as two independent original source scripts.
options(stringsAsFactors = FALSE)
if (file.exists("config/use_local_r_libs.R")) source("config/use_local_r_libs.R")
set.seed(20260626)

required <- c(
  "data.table", "dplyr", "tidyr", "readr", "stringr", "ggplot2",
  "openxlsx", "metafor", "mclust", "GSVA", "GSEABase", "DESeq2",
  "AnnotationDbi", "org.Hs.eg.db", "limma", "cluster"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(metafor)
  library(mclust)
  library(GSVA)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

CONFIG <- list(
  out_dir = "outputs/user_rerun/External_ECM_multicohort_analysis",
  instruction_dir = "data_manifest",
  user_input_dir = "data/public/metadata",
  gse18965_series_expr = "data/public/GSE18965/series_gene_level_expression.csv",
  gse18965_series_meta = "data/public/GSE18965/series_matrix_sample_metadata.csv",
  gse18965_cel_expr = "data/public/GSE18965/GSE18965_GCRMA_gene_expression.csv",
  gse18965_cel_meta = "data/public/GSE18965/sample_metadata.csv",
  gse118_counts = "data/public/GSE118761/GSE118761_genecounts.csv.gz",
  gse118_series = "data/public/GSE118761/GSE118761_series_matrix.txt.gz",
  gse118_scores = "data/public/GSE118761/GSE118761_IFN_T2_original32_PCA_cluster_scores.csv",
  gse152_counts = "data/public/GSE152004/GSE152004_695_raw_counts.txt.gz",
  gse152_series = "data/public/GSE152004/GSE152004_series_matrix.txt.gz",
  gse152_scores = "data/public/GSE152004/GSE152004_metadata_with_T2_IFN_cluster.csv",
  reactome_gmt = "gene_sets/Reactome_ECM_manifest.csv",
  hallmark_gmt = "gene_sets/Hallmark_EMT_manifest.csv",
  kegg_gmt = "gene_sets/KEGG_ECM_receptor_manifest.csv",
  c2_rds = "gene_sets/msigdb_C2_optional.rds",
  manifest = "data_manifest/public_data_manifest.csv",
  original32_file = "gene_sets/repair_ECM_fixed32.txt",
  custom_modules_gmt = "gene_sets/article_T2_IFN_Repair_custom_modules.gmt"
)

subdirs <- c(
  "00_logs", "01_input_mapping", "02_gene_sets", "02_gene_sets/standardized_inputs",
  "03_GSE18965_series", "04_GSE18965_CEL_sensitivity", "05_GSE118761_nasal",
  "06_GSE118761_tracheal", "07_GSE152004", "08_original32_comparison",
  "09_cross_cohort_meta", "10_downstream_sensitivity", "11_figures_diagnostic",
  "12_report"
)

safe_reset_dir <- function(path) {
  target <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(path) && !startsWith(target, paste0(root, "/"))) {
    stop("Refusing to delete output outside workspace: ", target)
  }
  if (file.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
safe_reset_dir(CONFIG$out_dir)
for (d in subdirs) dir.create(file.path(CONFIG$out_dir, d), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(CONFIG$out_dir, "00_logs", "run_log.txt")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_msg("Started External ECM multicohort analysis.")

file_manifest <- data.frame(
  role = names(CONFIG),
  path = unlist(CONFIG, use.names = FALSE),
  stringsAsFactors = FALSE
) |>
  subset(grepl("\\.(csv|txt|gz|gmt|rds|xlsx|zip)$|_dir$|file$|manifest$", path, ignore.case = TRUE))
file_manifest$exists <- file.exists(file_manifest$path)
file_manifest$size_bytes <- ifelse(file_manifest$exists, file.info(file_manifest$path)$size, NA_real_)
file_manifest$md5 <- ifelse(file_manifest$exists, as.character(tools::md5sum(file_manifest$path)), NA_character_)
write.csv(file_manifest, file.path(CONFIG$out_dir, "00_logs", "input_file_manifest.csv"), row.names = FALSE)

must_exist <- c(
  "gse18965_series_expr", "gse18965_series_meta", "gse18965_cel_expr", "gse18965_cel_meta",
  "gse118_counts", "gse118_series", "gse118_scores", "gse152_counts", "gse152_series",
  "gse152_scores", "reactome_gmt", "hallmark_gmt", "kegg_gmt", "c2_rds",
  "manifest", "original32_file", "custom_modules_gmt"
)
if (any(!file.exists(unlist(CONFIG[must_exist])))) {
  miss <- must_exist[!file.exists(unlist(CONFIG[must_exist]))]
  stop("Missing required input files: ", paste(miss, collapse = ", "))
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  out <- lapply(parts, function(z) unique(toupper(trimws(z[-c(1, 2)]))))
  names(out) <- vapply(parts, function(z) z[1], character(1))
  out
}

read_expr_csv <- function(path) {
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  gene_col <- names(x)[1]
  genes <- toupper(trimws(as.character(x[[gene_col]])))
  x[[gene_col]] <- NULL
  mat <- as.matrix(x)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes
  if (anyDuplicated(rownames(mat))) {
    means <- rowMeans(mat, na.rm = TRUE)
    keep <- unlist(tapply(seq_along(genes), genes, function(idx) idx[which.max(means[idx])]), use.names = FALSE)
    mat <- mat[sort(keep), , drop = FALSE]
  }
  if (anyNA(mat)) stop("NA in expression file: ", path)
  mat
}

normalize_binary <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.integer(x))
  z <- tolower(trimws(as.character(x)))
  out <- ifelse(z %in% c("1", "yes", "true", "atopic", "asthmatic", "wheezer", "high"), 1L,
                ifelse(z %in% c("0", "no", "false", "non-atopic", "non-asthmatic", "non-wheezer", "low"), 0L, NA_integer_))
  out
}

sample_id_from <- function(df, label) {
  if ("sample_id" %in% names(df)) return(as.character(df$sample_id))
  if ("gsm" %in% names(df)) return(as.character(df$gsm))
  if ("GSM" %in% names(df)) return(as.character(df$GSM))
  stop("No sample_id/gsm column in metadata: ", label)
}

parse_gse118_series <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[grepl("^!Sample_", lines)]
  keys <- sub("\t.*$", "", lines)
  vals <- strsplit(sub("^[^\t]+\t", "", lines), "\t")
  vals <- lapply(vals, function(v) gsub("^\"|\"$", "", v))
  n <- max(lengths(vals))
  df <- as.data.frame(matrix(NA_character_, n, length(vals)), stringsAsFactors = FALSE)
  names(df) <- make.unique(sub("^!Sample_", "", keys))
  for (i in seq_along(vals)) df[seq_along(vals[[i]]), i] <- vals[[i]]
  df$sample_id <- sub("^Sample name: ", "", df$description)
  df$diagnosis <- sub("^diagnosis: ", "", df$characteristics_ch1)
  df$tissue <- sub("^tissue: ", "", df$characteristics_ch1.1)
  df$subject_id <- sub("^subject_id: ", "", df$characteristics_ch1.2)
  df$sex <- sub("^gender: ", "", df$characteristics_ch1.3)
  df$age <- as.numeric(sub("^age \\(yrs\\): ", "", df$characteristics_ch1.4))
  df$atopy <- as.integer(grepl("^Atopic", df$diagnosis))
  df$asthma <- as.integer(grepl("asthmatic$", df$diagnosis))
  df$wheeze <- as.integer(grepl("wheezer", df$title) & !grepl("non-wheezer", df$title))
  df
}

read_count_cols <- function(path, sep) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  line <- readLines(con, n = 1, warn = FALSE)
  strsplit(line, sep, fixed = TRUE)[[1]] |> gsub("^\"|\"$", "", x = _)
}

prepare_count_vst <- function(counts, symbols, label, min_nonzero_fraction = 0.20) {
  symbols <- toupper(trimws(as.character(symbols)))
  keep_symbol <- !is.na(symbols) & nzchar(symbols) & symbols != "NA" & !grepl("///", symbols, fixed = TRUE)
  counts <- counts[keep_symbol, , drop = FALSE]
  symbols <- symbols[keep_symbol]
  keep_expr <- rowSums(counts > 0, na.rm = TRUE) >= ceiling(min_nonzero_fraction * ncol(counts))
  counts <- counts[keep_expr, , drop = FALSE]
  symbols <- symbols[keep_expr]
  counts <- rowsum(counts, group = symbols, reorder = FALSE)
  counts <- round(as.matrix(counts))
  storage.mode(counts) <- "integer"
  coldata <- data.frame(row.names = colnames(counts), sample_id = colnames(counts))
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  vst_mat <- assay(varianceStabilizingTransformation(dds, blind = TRUE))
  if (anyNA(vst_mat)) stop("NA in VST matrix: ", label)
  list(counts = counts, vst = vst_mat)
}

read_gse118_counts_vst <- function(path) {
  log_msg("Reading GSE118761 raw counts.")
  x <- read.csv(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  samples <- setdiff(names(x), c("EnsemblID", "Symbol", "Description"))
  counts <- as.matrix(x[, samples, drop = FALSE])
  storage.mode(counts) <- "numeric"
  list(raw = counts, symbols = x$Symbol, vst = prepare_count_vst(counts, x$Symbol, "GSE118761")$vst)
}

read_gse152_counts_vst <- function(path) {
  log_msg("Reading GSE152004 raw counts.")
  x <- read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
  counts <- as.matrix(x)
  storage.mode(counts) <- "numeric"
  list(raw = counts, symbols = rownames(counts), vst = prepare_count_vst(counts, rownames(counts), "GSE152004")$vst)
}

score_mean_z <- function(expr, genes, min_genes = 15, min_coverage = 0.20) {
  genes <- unique(toupper(genes))
  present <- intersect(genes, rownames(expr))
  coverage <- length(present) / length(genes)
  if (length(present) < min_genes || coverage < min_coverage) {
    return(list(score = rep(NA_real_, ncol(expr)), genes_present = present, n_present = length(present),
                n_total = length(genes), coverage = coverage, estimable = FALSE))
  }
  m <- expr[present, , drop = FALSE]
  mu <- rowMeans(m)
  sdv <- apply(m, 1, sd)
  keep <- is.finite(sdv) & sdv > 0
  z <- sweep(m[keep, , drop = FALSE], 1, mu[keep], "-")
  z <- sweep(z, 1, sdv[keep], "/")
  list(score = colMeans(z, na.rm = TRUE), genes_present = rownames(z), n_present = nrow(z),
       n_total = length(genes), coverage = nrow(z) / length(genes), estimable = TRUE)
}

score_ssgsea <- function(expr, genes, min_genes = 15, min_coverage = 0.20) {
  genes <- unique(toupper(genes))
  present <- intersect(genes, rownames(expr))
  coverage <- length(present) / length(genes)
  if (length(present) < min_genes || coverage < min_coverage) return(rep(NA_real_, ncol(expr)))
  gs <- list(external_set = present)
  ans <- try(GSVA::gsva(expr, gs, method = "ssgsea", min.sz = min_genes, max.sz = Inf,
                        ssgsea.norm = TRUE, verbose = FALSE), silent = TRUE)
  if (!inherits(ans, "try-error")) {
    return(as.numeric(ans[1, ]))
  }
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    par <- GSVA::ssgseaParam(exprData = expr, geneSets = gs, minSize = min_genes,
                             maxSize = Inf, normalize = TRUE)
    ans <- GSVA::gsva(par, verbose = FALSE)
    return(as.numeric(ans[1, ]))
  }
  log_msg("GSVA ssGSEA call failed; using deterministic in-script ssGSEA fallback for compatibility.")
  raw <- apply(expr, 2, function(x) {
    ord <- order(x, decreasing = TRUE, na.last = NA)
    ranked_genes <- rownames(expr)[ord]
    hits <- ranked_genes %in% present
    nh <- sum(hits)
    nm <- length(hits) - nh
    if (nh == 0 || nm == 0) return(NA_real_)
    ranks <- rank(x, ties.method = "average", na.last = "keep")[ord]
    hit_weight <- ifelse(hits, abs(ranks)^0.25, 0)
    phit <- cumsum(hit_weight) / sum(hit_weight)
    pmiss <- cumsum(ifelse(!hits, 1 / nm, 0))
    sum(phit - pmiss)
  })
  rng <- range(raw, na.rm = TRUE)
  if (is.finite(diff(rng)) && diff(rng) > 0) raw <- (raw - rng[1]) / diff(rng)
  as.numeric(raw)
}

group_test <- function(score, group, positive_level, reference_level, dataset, gene_set, method) {
  ok <- is.finite(score) & !is.na(group)
  score <- score[ok]
  group <- as.character(group[ok])
  if (!all(c(positive_level, reference_level) %in% group)) {
    return(data.frame(dataset, gene_set, method, positive_level, reference_level, n_positive = sum(group == positive_level),
                      n_reference = sum(group == reference_level), estimable = FALSE))
  }
  y1 <- score[group == positive_level]
  y0 <- score[group == reference_level]
  tt <- t.test(y1, y0, var.equal = FALSE)
  sp <- sqrt(((length(y1) - 1) * var(y1) + (length(y0) - 1) * var(y0)) / (length(y1) + length(y0) - 2))
  d <- (mean(y1) - mean(y0)) / sp
  j <- 1 - 3 / (4 * (length(y1) + length(y0)) - 9)
  data.frame(dataset, gene_set, method, positive_level, reference_level,
             n_positive = length(y1), n_reference = length(y0), estimable = TRUE,
             mean_positive = mean(y1), sd_positive = sd(y1),
             mean_reference = mean(y0), sd_reference = sd(y0),
             mean_difference = mean(y1) - mean(y0), ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
             welch_p = tt$p.value, mann_whitney_p = wilcox.test(y1, y0, exact = FALSE)$p.value,
             cohens_d = d, hedges_g = d * j)
}

logistic_test <- function(meta, score, outcome, covariates, dataset, gene_set, method, model_label) {
  d <- meta
  d$score_z <- as.numeric(scale(score))
  vars <- unique(c(outcome, "score_z", covariates))
  missing_vars <- setdiff(vars, names(d))
  if (length(missing_vars)) {
    return(data.frame(dataset, gene_set, method, outcome, model_label, covariates = paste(covariates, collapse = "+"),
                      n = 0, events = NA_integer_, estimable = FALSE, reason = paste("missing variables", paste(missing_vars, collapse = ","))))
  }
  d <- d[complete.cases(d[, vars, drop = FALSE]), vars, drop = FALSE]
  if (nrow(d) < 20 || length(unique(d[[outcome]])) != 2) {
    return(data.frame(dataset, gene_set, method, outcome, model_label, covariates = paste(covariates, collapse = "+"),
                      n = nrow(d), events = if (nrow(d)) sum(d[[outcome]] == 1) else NA_integer_,
                      estimable = FALSE, reason = "n<20 or outcome lacks two levels"))
  }
  fit <- glm(reformulate(c("score_z", covariates), response = outcome), data = d, family = binomial())
  co <- summary(fit)$coefficients
  beta <- co["score_z", "Estimate"]
  se <- co["score_z", "Std. Error"]
  p <- co["score_z", "Pr(>|z|)"]
  data.frame(dataset, gene_set, method, outcome, model_label, covariates = paste(covariates, collapse = "+"),
             n = nrow(d), events = sum(d[[outcome]] == 1), estimable = TRUE, reason = "",
             beta_per_sd = beta, se = se, OR_per_sd = exp(beta), CI_low = exp(beta - 1.96 * se),
             CI_high = exp(beta + 1.96 * se), p_value = p, AIC = AIC(fit))
}

safe_cor <- function(x, y, dataset, gene_set, comparison, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(data.frame(dataset, gene_set, comparison, method, n = sum(ok), estimate = NA_real_, p_value = NA_real_))
  ans <- cor.test(x[ok], y[ok], method = method, exact = FALSE)
  data.frame(dataset, gene_set, comparison, method, n = sum(ok), estimate = unname(ans$estimate), p_value = ans$p.value)
}

kappa2 <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  if (sum(ok) < 5) return(c(agreement = NA_real_, kappa = NA_real_, n = sum(ok)))
  a <- as.character(a[ok]); b <- as.character(b[ok])
  tab <- table(a, b)
  po <- sum(diag(tab)) / sum(tab)
  pe <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
  c(agreement = po, kappa = (po - pe) / (1 - pe), n = sum(tab))
}

nmi <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  a <- as.factor(a[ok]); b <- as.factor(b[ok])
  if (length(a) < 2) return(NA_real_)
  tab <- table(a, b)
  pxy <- tab / sum(tab); px <- rowSums(pxy); py <- colSums(pxy)
  mi <- sum(pxy[pxy > 0] * log(pxy[pxy > 0] / outer(px, py)[pxy > 0]))
  hx <- -sum(px[px > 0] * log(px[px > 0])); hy <- -sum(py[py > 0] * log(py[py > 0]))
  mi / sqrt(hx * hy)
}

add_sheet <- function(wb, name, df) {
  nm <- substr(name, 1, 31)
  addWorksheet(wb, nm)
  writeData(wb, nm, df)
  freezePane(wb, nm, firstRow = TRUE)
  setColWidths(wb, nm, cols = 1:min(12, max(1, ncol(df))), widths = "auto")
}

log_msg("Loading fixed external gene sets.")
manifest <- fread(CONFIG$manifest, data.table = FALSE)
gmt_reactome <- read_gmt(CONFIG$reactome_gmt)
gmt_hallmark <- read_gmt(CONFIG$hallmark_gmt)
gmt_kegg <- read_gmt(CONFIG$kegg_gmt)
c2 <- as.data.table(readRDS(CONFIG$c2_rds))
get_set <- function(name) {
  if (name %in% names(gmt_reactome)) return(gmt_reactome[[name]])
  if (name %in% names(gmt_hallmark)) return(gmt_hallmark[[name]])
  if (name %in% names(gmt_kegg)) return(gmt_kegg[[name]])
  if (name %in% c2$gs_name) return(unique(toupper(c2[gs_name == name, db_gene_symbol])))
  stop("Prespecified gene set not found: ", name)
}
gene_sets <- setNames(lapply(manifest$gene_set_name, get_set), manifest$gene_set_name)
original32 <- unique(toupper(trimws(fread(CONFIG$original32_file, data.table = FALSE)[[1]])))

gene_set_members <- do.call(rbind, lapply(names(gene_sets), function(nm) {
  data.frame(gene_set = nm, gene_symbol = sort(gene_sets[[nm]]), stringsAsFactors = FALSE)
}))
write.csv(gene_set_members, file.path(CONFIG$out_dir, "02_gene_sets", "external_ECM_gene_sets_original.csv"), row.names = FALSE)
wb_gs <- createWorkbook()
add_sheet(wb_gs, "gene_sets_original", gene_set_members)
add_sheet(wb_gs, "manifest", manifest)
saveWorkbook(wb_gs, file.path(CONFIG$out_dir, "02_gene_sets", "external_ECM_gene_sets_original.xlsx"), overwrite = TRUE)

log_msg("Loading and standardizing GSE18965 Series and CEL matrices.")
expr18965_series <- read_expr_csv(CONFIG$gse18965_series_expr)
meta18965_series_raw <- fread(CONFIG$gse18965_series_meta, data.table = FALSE, check.names = FALSE)
meta18965_series <- data.frame(sample_id = sample_id_from(meta18965_series_raw, "GSE18965 Series"),
                               group = meta18965_series_raw$group,
                               tissue = "bronchial epithelial brushing",
                               stringsAsFactors = FALSE)
expr18965_series <- expr18965_series[, meta18965_series$sample_id, drop = FALSE]
expr18965_cel <- read_expr_csv(CONFIG$gse18965_cel_expr)
meta18965_cel_raw <- fread(CONFIG$gse18965_cel_meta, data.table = FALSE, check.names = FALSE)
meta18965_cel <- data.frame(sample_id = sample_id_from(meta18965_cel_raw, "GSE18965 CEL"),
                            group = meta18965_cel_raw$group,
                            tissue = "bronchial epithelial brushing",
                            stringsAsFactors = FALSE)
expr18965_cel <- expr18965_cel[, meta18965_cel$sample_id, drop = FALSE]

log_msg("Loading GSE118761 metadata and raw counts -> VST.")
meta118_scores <- fread(CONFIG$gse118_scores, data.table = FALSE, check.names = FALSE)
meta118 <- data.frame(
  sample_id = meta118_scores$sample_name,
  subject_id = meta118_scores$subject_id,
  tissue = ifelse(grepl("Nasal", meta118_scores$tissue), "nasal", "tracheal"),
  baseline = TRUE,
  asthma = as.integer(meta118_scores$asthma_binary),
  atopy = as.integer(meta118_scores$atopy_binary),
  wheeze = as.integer(meta118_scores$wheeze_binary_exact),
  age = as.numeric(meta118_scores$age_years),
  sex = factor(meta118_scores$gender),
  batch = NA_character_,
  T2_score = as.numeric(meta118_scores$t2_score_historical_log2Count),
  IFN_score = as.numeric(meta118_scores$ifn_score_existing),
  original32_score = as.numeric(meta118_scores$original32_repair_score_log2CPM),
  repair_factor = as.numeric(meta118_scores$Factor1),
  original_projected_cluster = meta118_scores$historical_projected_endotype_nasal_only,
  stringsAsFactors = FALSE
)
for (j in paste0("Factor", 1:10)) meta118[[j]] <- meta118_scores[[j]]
gse118 <- read_gse118_counts_vst(CONFIG$gse118_counts)
if (!setequal(colnames(gse118$vst), meta118$sample_id)) stop("GSE118761 sample IDs do not match metadata.")
gse118$vst <- gse118$vst[, meta118$sample_id, drop = FALSE]

log_msg("Loading GSE152004 metadata and raw counts -> VST.")
meta152_scores <- fread(CONFIG$gse152_scores, data.table = FALSE, check.names = FALSE)
meta152 <- data.frame(
  sample_id = meta152_scores$sample_id,
  subject_id = meta152_scores$sample_id,
  tissue = "nasal",
  asthma = as.integer(meta152_scores$asthma_binary),
  atopy = NA_integer_,
  T2_high = as.integer(meta152_scores$t2_high_historical_median),
  T2_score = as.numeric(meta152_scores$t2_score_historical_log2Count),
  IFN_score = as.numeric(meta152_scores$ifn_score_existing),
  original32_score = as.numeric(meta152_scores$repair_score_existing),
  age = NA_real_,
  sex = NA_character_,
  batch = NA_character_,
  original_cluster = meta152_scores$historical_endotype,
  stringsAsFactors = FALSE
)
gse152 <- read_gse152_counts_vst(CONFIG$gse152_counts)
if (!setequal(colnames(gse152$vst), meta152$sample_id)) stop("GSE152004 sample IDs do not match metadata.")
gse152$vst <- gse152$vst[, meta152$sample_id, drop = FALSE]

write.csv(data.frame(gene_symbol = rownames(expr18965_series), expr18965_series, check.names = FALSE),
          file.path(CONFIG$out_dir, "02_gene_sets/standardized_inputs", "GSE18965_series_gene_expression.csv"), row.names = FALSE)
write.csv(meta18965_series, file.path(CONFIG$out_dir, "02_gene_sets/standardized_inputs", "GSE18965_series_metadata.csv"), row.names = FALSE)
write.csv(data.frame(gene_symbol = rownames(expr18965_cel), expr18965_cel, check.names = FALSE),
          file.path(CONFIG$out_dir, "02_gene_sets/standardized_inputs", "GSE18965_CEL_gene_expression.csv"), row.names = FALSE)
write.csv(meta18965_cel, file.path(CONFIG$out_dir, "02_gene_sets/standardized_inputs", "GSE18965_CEL_metadata.csv"), row.names = FALSE)
write.csv(data.frame(gene_symbol = rownames(gse118$vst), gse118$vst, check.names = FALSE),
          file.path(CONFIG$out_dir, "02_gene_sets/standardized_inputs", "GSE118761_VST_gene_expression.csv"), row.names = FALSE)
write.csv(meta118, file.path(CONFIG$out_dir, "02_gene_sets/standardized_inputs", "GSE118761_metadata.csv"), row.names = FALSE)
write.csv(data.frame(gene_symbol = rownames(gse152$vst), gse152$vst, check.names = FALSE),
          file.path(CONFIG$out_dir, "02_gene_sets/standardized_inputs", "GSE152004_VST_gene_expression.csv"), row.names = FALSE)
write.csv(meta152, file.path(CONFIG$out_dir, "02_gene_sets/standardized_inputs", "GSE152004_metadata.csv"), row.names = FALSE)

cohorts <- list(
  GSE18965_series = list(expr = expr18965_series, meta = meta18965_series, primary = TRUE),
  GSE18965_CEL_sensitivity = list(expr = expr18965_cel, meta = meta18965_cel, primary = FALSE),
  GSE118761_nasal = list(expr = gse118$vst[, meta118$sample_id[meta118$tissue == "nasal"], drop = FALSE],
                         meta = meta118[meta118$tissue == "nasal", ], primary = TRUE),
  GSE118761_tracheal = list(expr = gse118$vst[, meta118$sample_id[meta118$tissue == "tracheal"], drop = FALSE],
                            meta = meta118[meta118$tissue == "tracheal", ], primary = FALSE),
  GSE152004 = list(expr = gse152$vst, meta = meta152, primary = TRUE)
)

input_lines <- c(
  "# INPUT_MAPPING",
  "",
  "## Global Rules",
  "- R runtime: Rscript.",
  "- R library bootstrap: source(\"config/use_local_r_libs.R\").",
  "- R library used: user-managed R library.",
  "- No cross-cohort expression matrix merging and no cross-cohort ComBat were performed.",
  "- External gene sets were locked before outcome modeling.",
  "",
  "## GSE18965 Series Matrix Main Analysis",
  paste0("- Expression matrix: ", CONFIG$gse18965_series_expr),
  "- Matrix type: GEO deposited Series Matrix reproduced gene-level log2 expression; rows are HGNC symbols; columns are GSM sample IDs.",
  "- Duplicate-gene handling: one duplicated IGK row in the prior gene-level export was resolved by retaining the row with the highest mean expression; this does not affect any external ECM gene set member.",
  paste0("- Metadata: ", CONFIG$gse18965_series_meta),
  "- Sample matching: 16/16 expression columns matched metadata sample_id.",
  "- Primary phenotype: group, standardized as AA versus HN; contrast label is atopic asthma versus non-atopic healthy.",
  "- Covariates: age/sex/batch not available in this dataset.",
  "",
  "## GSE18965 CEL Sensitivity",
  paste0("- Expression matrix: ", CONFIG$gse18965_cel_expr),
  "- Matrix type: CEL -> GC-RMA -> GPL96 annotation -> highest mean probe per gene; rows are HGNC symbols; columns are GSM sample IDs.",
  "- Duplicate-gene handling: one duplicated IGK row in the prior gene-level export was resolved by retaining the row with the highest mean expression; this does not affect any external ECM gene set member.",
  paste0("- Metadata: ", CONFIG$gse18965_cel_meta),
  "- Sample matching: 16/16 expression columns matched metadata sample_id.",
  "",
  "## GSE118761",
  paste0("- Raw counts: ", CONFIG$gse118_counts),
  paste0("- GEO Series Matrix metadata: ", CONFIG$gse118_series),
  paste0("- Supplementary metadata/scores: ", CONFIG$gse118_scores),
  "- Matrix type: raw gene-level counts with EnsemblID, Symbol, Description columns; filtered at nonzero in >=20% samples, duplicate symbols summed, then DESeq2 VST(blind=TRUE).",
  "- Sample ID source: count columns AFxxx match supplementary sample_name and GEO Series Matrix description 'Sample name: AFxxx'.",
  paste0("- Sample matching: ", ncol(gse118$vst), "/", nrow(meta118), " matched."),
  paste0("- Nasal samples: ", sum(meta118$tissue == "nasal"), "; tracheal samples: ", sum(meta118$tissue == "tracheal"), "."),
  "- Primary phenotype: atopy_binary -> atopy. Secondary phenotypes: asthma_binary -> asthma; wheeze_binary_exact -> wheeze.",
  "- Covariates present: age_years -> age; gender -> sex. Technical batch not available. Treatment/timepoint not available; all samples treated as baseline cross-sectional GEO samples.",
  "- Existing fixed scores present: ifn_score_existing -> IFN_score; t2_score_historical_log2Count -> T2_score; original32_repair_score_log2CPM -> original32_score.",
  "- PCA factors present: Factor1-Factor10; repair_factor set to Factor1. Nasal projected cluster present in historical_projected_endotype_nasal_only.",
  "- Duplicate subjects: one nasal and one tracheal sample for most subjects; nasal and tracheal are analyzed separately.",
  "",
  "## GSE152004",
  paste0("- Raw counts: ", CONFIG$gse152_counts),
  paste0("- GEO Series Matrix metadata: ", CONFIG$gse152_series),
  paste0("- Supplementary metadata/scores: ", CONFIG$gse152_scores),
  "- Matrix type: raw gene-level counts with HGNC symbol row identifiers; filtered at nonzero in >=20% samples, duplicate symbols summed, then DESeq2 VST(blind=TRUE).",
  "- Sample ID source: count column names HR/SJ IDs match GEO Series Matrix title and supplementary sample_id.",
  paste0("- Sample matching: ", ncol(gse152$vst), "/", nrow(meta152), " matched."),
  "- Primary phenotype: t2_high_historical_median -> T2_high. Secondary phenotype: asthma_binary -> asthma.",
  "- Atopy is not available and is reported as not estimable.",
  "- Age, sex and technical batch are not available; adjusted models are therefore not estimated for this cohort.",
  "- Existing fixed scores present: ifn_score_existing -> IFN_score; t2_score_historical_log2Count -> T2_score; repair_score_existing -> original32_score.",
  "- Historical k=3 cluster label present in historical_endotype -> original_cluster.",
  "",
  "## External Gene Sets",
  paste0("- Manifest: ", CONFIG$manifest),
  paste0("- Reactome GMT: ", CONFIG$reactome_gmt),
  paste0("- Hallmark GMT: ", CONFIG$hallmark_gmt),
  paste0("- KEGG legacy GMT: ", CONFIG$kegg_gmt),
  paste0("- NABA sets: ", CONFIG$c2_rds, " (MSigDB 2026.1.Hs C2 RDS; c2.cgp/c2.all GMT was not present locally)."),
  "- Custom project GMT is used only for fixed IFN/T2/original repair provenance, not as an external ECM set."
)
writeLines(input_lines, file.path(CONFIG$out_dir, "01_input_mapping", "INPUT_MAPPING.md"))
file.copy(file.path(CONFIG$out_dir, "01_input_mapping", "INPUT_MAPPING.md"),
          file.path(CONFIG$out_dir, "INPUT_MAPPING.md"), overwrite = TRUE)

coverage_rows <- list()
score_rows <- list()
score_lookup <- list()
original32_lookup <- list()

log_msg("Scoring external gene sets.")
for (cohort_name in names(cohorts)) {
  expr <- cohorts[[cohort_name]]$expr
  meta <- cohorts[[cohort_name]]$meta
  expr <- expr[, meta$sample_id, drop = FALSE]
  bg <- rownames(expr)
  original32_score <- score_mean_z(expr, original32)
  original32_lookup[[cohort_name]] <- original32_score$score
  for (set_name in names(gene_sets)) {
    genes <- gene_sets[[set_name]]
    overlap32 <- intersect(genes, original32)
    variants <- list(original = genes, excluding_original32 = setdiff(genes, original32))
    for (variant in names(variants)) {
      genes_use <- variants[[variant]]
      sc <- score_mean_z(expr, genes_use)
      miss <- setdiff(genes_use, rownames(expr))
      hyper_p <- phyper(length(intersect(genes_use, bg)) - 1, length(bg), length(union(bg, genes_use)) - length(bg),
                        length(genes_use), lower.tail = FALSE)
      coverage_rows[[length(coverage_rows) + 1]] <- data.frame(
        cohort = cohort_name, gene_set = set_name, variant = variant,
        original_member_count = length(genes), scoring_member_count = length(genes_use),
        detected_member_count = sc$n_present, coverage = sc$coverage,
        estimable = sc$estimable, missing_genes = paste(sort(miss), collapse = ";"),
        overlap_with_original32 = paste(sort(overlap32), collapse = ";"),
        n_overlap_original32 = length(overlap32),
        jaccard_with_original32 = length(overlap32) / length(union(genes, original32)),
        hypergeometric_overlap_p = hyper_p,
        genes_used = paste(sort(sc$genes_present), collapse = ";"),
        stringsAsFactors = FALSE
      )
      score_lookup[[paste(cohort_name, set_name, variant, "mean_z", sep = "|")]] <- sc$score
      score_rows[[length(score_rows) + 1]] <- data.frame(
        cohort = cohort_name, sample_id = meta$sample_id, gene_set = set_name,
        variant = variant, scoring_method = "mean_z", score = sc$score,
        n_genes_used = sc$n_present, coverage = sc$coverage, estimable = sc$estimable,
        genes_used = paste(sort(sc$genes_present), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
    ss <- score_ssgsea(expr, genes)
    score_lookup[[paste(cohort_name, set_name, "original", "ssGSEA", sep = "|")]] <- ss
    score_rows[[length(score_rows) + 1]] <- data.frame(
      cohort = cohort_name, sample_id = meta$sample_id, gene_set = set_name,
      variant = "original", scoring_method = "ssGSEA", score = ss,
      n_genes_used = length(intersect(genes, rownames(expr))),
      coverage = length(intersect(genes, rownames(expr))) / length(genes),
      estimable = all(is.finite(ss)), genes_used = paste(sort(intersect(genes, rownames(expr))), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}
coverage <- bind_rows(coverage_rows)
sample_scores <- bind_rows(score_rows)
write.csv(coverage, file.path(CONFIG$out_dir, "02_gene_sets", "external_ECM_gene_coverage_by_cohort.csv"), row.names = FALSE)
write.csv(sample_scores, file.path(CONFIG$out_dir, "02_gene_sets", "sample_level_external_ECM_scores.csv"), row.names = FALSE)

log_msg("Running prespecified cohort models.")
gse18965_results <- list()
for (cohort_name in c("GSE18965_series", "GSE18965_CEL_sensitivity")) {
  meta <- cohorts[[cohort_name]]$meta
  for (set_name in names(gene_sets)) for (variant in c("original", "excluding_original32")) {
    score <- score_lookup[[paste(cohort_name, set_name, variant, "mean_z", sep = "|")]]
    gse18965_results[[length(gse18965_results) + 1]] <- group_test(score, meta$group, "AA", "HN", cohort_name, set_name, paste0("mean_z_", variant))
    ss <- score_lookup[[paste(cohort_name, set_name, "original", "ssGSEA", sep = "|")]]
    if (variant == "original") gse18965_results[[length(gse18965_results) + 1]] <- group_test(ss, meta$group, "AA", "HN", cohort_name, set_name, "ssGSEA_original")
  }
}
gse18965_results <- bind_rows(gse18965_results)
gse18965_results$BH_FDR_within_dataset <- ave(gse18965_results$welch_p, gse18965_results$dataset, FUN = function(p) p.adjust(p, "BH"))

logistic_results <- list()
for (cohort_name in c("GSE118761_nasal", "GSE118761_tracheal")) {
  meta <- cohorts[[cohort_name]]$meta
  for (set_name in names(gene_sets)) for (variant in c("original", "excluding_original32")) {
    score <- score_lookup[[paste(cohort_name, set_name, variant, "mean_z", sep = "|")]]
    for (outcome in c("atopy", "asthma", "wheeze")) {
      logistic_results[[length(logistic_results) + 1]] <- logistic_test(meta, score, outcome, character(), cohort_name, set_name, paste0("mean_z_", variant), "unadjusted")
      logistic_results[[length(logistic_results) + 1]] <- logistic_test(meta, score, outcome, c("age", "sex"), cohort_name, set_name, paste0("mean_z_", variant), "age_sex_adjusted")
    }
    ss <- score_lookup[[paste(cohort_name, set_name, "original", "ssGSEA", sep = "|")]]
    for (outcome in c("atopy", "asthma", "wheeze")) {
      logistic_results[[length(logistic_results) + 1]] <- logistic_test(meta, ss, outcome, c("age", "sex"), cohort_name, set_name, "ssGSEA_original", "age_sex_adjusted")
    }
  }
}
meta <- cohorts$GSE152004$meta
for (set_name in names(gene_sets)) for (variant in c("original", "excluding_original32")) {
  score <- score_lookup[[paste("GSE152004", set_name, variant, "mean_z", sep = "|")]]
  for (outcome in c("T2_high", "asthma", "atopy")) {
    logistic_results[[length(logistic_results) + 1]] <- logistic_test(meta, score, outcome, character(), "GSE152004", set_name, paste0("mean_z_", variant), "unadjusted")
  }
  ss <- score_lookup[[paste("GSE152004", set_name, "original", "ssGSEA", sep = "|")]]
  for (outcome in c("T2_high", "asthma", "atopy")) {
    logistic_results[[length(logistic_results) + 1]] <- logistic_test(meta, ss, outcome, character(), "GSE152004", set_name, "ssGSEA_original", "unadjusted")
  }
}
logistic_results <- bind_rows(logistic_results)
logistic_results$BH_FDR_within_dataset <- ave(ifelse(logistic_results$estimable, logistic_results$p_value, NA_real_),
                                              logistic_results$dataset,
                                              FUN = function(p) p.adjust(p, "BH"))

primary_fam <- bind_rows(
  gse18965_results |> filter(dataset == "GSE18965_series", gene_set == "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", method == "mean_z_original") |>
    transmute(family = "Family1_primary_REACTOME", dataset, gene_set, endpoint = "AA_vs_HN", effect = mean_difference, p_value = welch_p),
  logistic_results |> filter(dataset == "GSE118761_nasal", gene_set == "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
                             method == "mean_z_original", outcome == "atopy", model_label == "age_sex_adjusted", estimable) |>
    transmute(family = "Family1_primary_REACTOME", dataset, gene_set, endpoint = "atopy", effect = beta_per_sd, p_value),
  logistic_results |> filter(dataset == "GSE152004", gene_set == "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
                             method == "mean_z_original", outcome == "T2_high", model_label == "unadjusted", estimable) |>
    transmute(family = "Family1_primary_REACTOME", dataset, gene_set, endpoint = "T2_high", effect = beta_per_sd, p_value)
)
primary_fam$BH_FDR <- p.adjust(primary_fam$p_value, "BH")

core_fam <- bind_rows(
  gse18965_results |> filter(dataset == "GSE18965_series", gene_set == "NABA_CORE_MATRISOME", method == "mean_z_original") |>
    transmute(family = "Family2_confirmatory_NABA_CORE", dataset, gene_set, endpoint = "AA_vs_HN", effect = mean_difference, p_value = welch_p),
  logistic_results |> filter(dataset == "GSE118761_nasal", gene_set == "NABA_CORE_MATRISOME",
                             method == "mean_z_original", outcome == "atopy", model_label == "age_sex_adjusted", estimable) |>
    transmute(family = "Family2_confirmatory_NABA_CORE", dataset, gene_set, endpoint = "atopy", effect = beta_per_sd, p_value),
  logistic_results |> filter(dataset == "GSE152004", gene_set == "NABA_CORE_MATRISOME",
                             method == "mean_z_original", outcome == "T2_high", model_label == "unadjusted", estimable) |>
    transmute(family = "Family2_confirmatory_NABA_CORE", dataset, gene_set, endpoint = "T2_high", effect = beta_per_sd, p_value)
)
core_fam$BH_FDR <- p.adjust(core_fam$p_value, "BH")
fdr_families <- bind_rows(primary_fam, core_fam)

write.csv(gse18965_results, file.path(CONFIG$out_dir, "03_GSE18965_series", "GSE18965_external_ECM_results_combined.csv"), row.names = FALSE)
write.csv(gse18965_results |> filter(dataset == "GSE18965_CEL_sensitivity"),
          file.path(CONFIG$out_dir, "04_GSE18965_CEL_sensitivity", "GSE18965_CEL_external_ECM_results.csv"), row.names = FALSE)
write.csv(logistic_results |> filter(grepl("GSE118761", dataset)),
          file.path(CONFIG$out_dir, "05_GSE118761_nasal", "GSE118761_external_ECM_results_combined.csv"), row.names = FALSE)
write.csv(logistic_results |> filter(dataset == "GSE152004"),
          file.path(CONFIG$out_dir, "07_GSE152004", "GSE152004_external_ECM_results.csv"), row.names = FALSE)

log_msg("Comparing external ECM scores with original 32-gene score.")
compare_rows <- list()
for (cohort_name in names(cohorts)) {
  orig <- original32_lookup[[cohort_name]]
  orig_hi <- ifelse(orig >= median(orig, na.rm = TRUE), "high", "low")
  for (set_name in names(gene_sets)) {
    ext <- score_lookup[[paste(cohort_name, set_name, "original", "mean_z", sep = "|")]]
    ext_ex <- score_lookup[[paste(cohort_name, set_name, "excluding_original32", "mean_z", sep = "|")]]
    ext_hi <- ifelse(ext >= median(ext, na.rm = TRUE), "high", "low")
    kap <- kappa2(orig_hi, ext_hi)
    compare_rows[[length(compare_rows) + 1]] <- bind_rows(
      safe_cor(ext, orig, cohort_name, set_name, "external_vs_original32", "pearson"),
      safe_cor(ext, orig, cohort_name, set_name, "external_vs_original32", "spearman")
    ) |>
      mutate(overlap_genes = paste(sort(intersect(gene_sets[[set_name]], original32)), collapse = ";"),
             n_overlap = length(intersect(gene_sets[[set_name]], original32)),
             high_low_agreement = kap[["agreement"]], cohen_kappa = kap[["kappa"]],
             high_low_n = kap[["n"]],
             excluding_original32_direction_available = any(is.finite(ext_ex)))
  }
}
external_vs_original32 <- bind_rows(compare_rows)
external_vs_original32$BH_FDR_exploratory <- ave(external_vs_original32$p_value, external_vs_original32$dataset,
                                                 FUN = function(p) p.adjust(p, "BH"))
write.csv(external_vs_original32, file.path(CONFIG$out_dir, "08_original32_comparison", "external_ECM_vs_original32.csv"), row.names = FALSE)

log_msg("Running cross-cohort meta-analysis and direction summaries.")
direction_table <- bind_rows(
  gse18965_results |> filter(method == "mean_z_original") |>
    transmute(dataset, gene_set, endpoint = "AA_vs_HN", effect_type = "mean_difference_or_logOR",
              effect = mean_difference, p_value = welch_p, direction = ifelse(effect > 0, "higher_in_positive", "lower_in_positive")),
  logistic_results |> filter(method == "mean_z_original", estimable, model_label %in% c("age_sex_adjusted", "unadjusted")) |>
    transmute(dataset, gene_set, endpoint = outcome, effect_type = "mean_difference_or_logOR",
              effect = beta_per_sd, p_value, direction = ifelse(effect > 0, "OR_greater_than_1", "OR_less_than_1"))
)
meta_rows <- list()
for (set_name in names(gene_sets)) {
  ast <- logistic_results |> filter(gene_set == set_name, method == "mean_z_original", outcome == "asthma",
                                    dataset %in% c("GSE118761_nasal", "GSE152004"),
                                    model_label %in% c("age_sex_adjusted", "unadjusted"), estimable)
  if (nrow(ast) >= 2) {
    fit <- metafor::rma(yi = ast$beta_per_sd, sei = ast$se, method = "REML")
    meta_rows[[length(meta_rows) + 1]] <- data.frame(contrast = "asthma_GSE118761_nasal_GSE152004",
                                                     gene_set = set_name, k = fit$k, beta = fit$b[1],
                                                     CI_low = fit$ci.lb, CI_high = fit$ci.ub,
                                                     p_value = fit$pval, tau2 = fit$tau2, I2 = fit$I2)
  }
}
meta_analysis <- bind_rows(meta_rows)
if (nrow(meta_analysis)) meta_analysis$BH_FDR <- p.adjust(meta_analysis$p_value, "BH")
write.csv(direction_table, file.path(CONFIG$out_dir, "09_cross_cohort_meta", "cross_cohort_effect_direction_table.csv"), row.names = FALSE)
write.csv(meta_analysis, file.path(CONFIG$out_dir, "09_cross_cohort_meta", "cross_cohort_meta_analysis.csv"), row.names = FALSE)

log_msg("Running downstream cluster/projection sensitivity.")
reactome <- "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"
down_rows <- list()
quad_rows <- list()
meta118_nasal <- cohorts$GSE118761_nasal$meta
ecm118 <- score_lookup[[paste("GSE118761_nasal", reactome, "original", "mean_z", sep = "|")]]
ifn118 <- meta118_nasal$IFN_score
quad <- ifelse(ecm118 >= median(ecm118, na.rm = TRUE), "ECM_high", "ECM_low")
quad <- paste(quad, ifelse(ifn118 >= median(ifn118, na.rm = TRUE), "IFN_high", "IFN_low"), sep = "_")
for (q in sort(unique(quad))) {
  idx <- quad == q
  quad_rows[[length(quad_rows) + 1]] <- data.frame(dataset = "GSE118761_nasal", quadrant = q, n = sum(idx),
                                                   atopy_rate = mean(meta118_nasal$atopy[idx] == 1),
                                                   asthma_rate = mean(meta118_nasal$asthma[idx] == 1),
                                                   wheeze_rate = mean(meta118_nasal$wheeze[idx] == 1))
}

meta152_all <- cohorts$GSE152004$meta
ecm152 <- score_lookup[[paste("GSE152004", reactome, "original", "mean_z", sep = "|")]]
feat152 <- data.frame(T2 = meta152_all$T2_score, IFN = meta152_all$IFN_score, ECM = ecm152)
feat152_z <- scale(feat152)
set.seed(20260626)
km <- kmeans(feat152_z, centers = 3, nstart = 100)
new_cluster <- paste0("C", km$cluster)
sil <- cluster::silhouette(km$cluster, dist(feat152_z))
ari <- mclust::adjustedRandIndex(new_cluster, meta152_all$original_cluster)
nmi_val <- nmi(new_cluster, meta152_all$original_cluster)
cluster_summary <- data.frame(cluster = sort(unique(new_cluster))) |>
  rowwise() |>
  mutate(n = sum(new_cluster == cluster),
         mean_T2_z = mean(feat152_z[new_cluster == cluster, "T2"]),
         mean_IFN_z = mean(feat152_z[new_cluster == cluster, "IFN"]),
         mean_ECM_z = mean(feat152_z[new_cluster == cluster, "ECM"]),
         T2_high_rate = mean(meta152_all$T2_high[new_cluster == cluster] == 1),
         asthma_rate = mean(meta152_all$asthma[new_cluster == cluster] == 1),
         atopy_rate = NA_real_) |>
  ungroup()
cluster_metrics <- data.frame(dataset = "GSE152004", ARI_vs_original = ari, NMI_vs_original = nmi_val,
                              mean_silhouette = mean(sil[, "sil_width"]))
cluster_cross <- as.data.frame.matrix(table(new_cluster, meta152_all$original_cluster))
cluster_cross$new_cluster <- rownames(cluster_cross)
cluster_cross <- cluster_cross[, c("new_cluster", setdiff(names(cluster_cross), "new_cluster"))]

centers <- km$centers
mu152 <- attr(feat152_z, "scaled:center")
sd152 <- attr(feat152_z, "scaled:scale")
feat118 <- data.frame(T2 = meta118_nasal$T2_score, IFN = meta118_nasal$IFN_score, ECM = ecm118)
feat118_z <- sweep(feat118, 2, mu152, "-")
feat118_z <- sweep(feat118_z, 2, sd152, "/")
nearest <- apply(as.matrix(feat118_z), 1, function(v) {
  paste0("C", which.min(rowSums((centers - matrix(v, nrow = nrow(centers), ncol = ncol(centers), byrow = TRUE))^2)))
})
proj_summary <- data.frame(projected_cluster = sort(unique(nearest))) |>
  rowwise() |>
  mutate(n = sum(nearest == projected_cluster),
         atopy_rate = mean(meta118_nasal$atopy[nearest == projected_cluster] == 1),
         asthma_rate = mean(meta118_nasal$asthma[nearest == projected_cluster] == 1),
         wheeze_rate = mean(meta118_nasal$wheeze[nearest == projected_cluster] == 1)) |>
  ungroup()
proj_metrics <- data.frame(dataset = "GSE118761_nasal_projected_to_GSE152004_centroids",
                           ARI_vs_original_projected = mclust::adjustedRandIndex(nearest, meta118_nasal$original_projected_cluster),
                           NMI_vs_original_projected = nmi(nearest, meta118_nasal$original_projected_cluster))

downstream <- list(
  IFN_ECM_quadrants_GSE118761 = bind_rows(quad_rows),
  kmeans_cluster_summary_GSE152004 = cluster_summary,
  kmeans_metrics_GSE152004 = cluster_metrics,
  kmeans_cross_table_GSE152004 = cluster_cross,
  projection_summary_GSE118761 = proj_summary,
  projection_metrics_GSE118761 = proj_metrics
)
write.csv(downstream$IFN_ECM_quadrants_GSE118761, file.path(CONFIG$out_dir, "10_downstream_sensitivity", "GSE118761_IFN_ECM_quadrants.csv"), row.names = FALSE)
write.csv(downstream$kmeans_cluster_summary_GSE152004, file.path(CONFIG$out_dir, "10_downstream_sensitivity", "GSE152004_k3_cluster_summary.csv"), row.names = FALSE)
write.csv(downstream$projection_summary_GSE118761, file.path(CONFIG$out_dir, "10_downstream_sensitivity", "GSE118761_projection_summary.csv"), row.names = FALSE)

log_msg("Applying final support classification.")
get_pf <- function(ds, endpoint) primary_fam |> filter(dataset == ds, endpoint == endpoint)
pf18965 <- get_pf("GSE18965_series", "AA_vs_HN")
pf118 <- get_pf("GSE118761_nasal", "atopy")
pf152 <- get_pf("GSE152004", "T2_high")
reactome_series_positive <- nrow(pf18965) && is.finite(pf18965$effect[1]) && pf18965$effect[1] > 0
independent_sig <- (nrow(pf118) && pf118$effect[1] > 0 && pf118$BH_FDR[1] < 0.05) ||
  (nrow(pf152) && pf152$effect[1] > 0 && pf152$BH_FDR[1] < 0.05)
other_consistent <- (nrow(pf118) && pf118$effect[1] > 0) && (nrow(pf152) && pf152$effect[1] > 0)
ex_series <- gse18965_results |> filter(dataset == "GSE18965_series", gene_set == reactome, method == "mean_z_excluding_original32")
ex_118 <- logistic_results |> filter(dataset == "GSE118761_nasal", gene_set == reactome,
                                     method == "mean_z_excluding_original32", outcome == "atopy",
                                     model_label == "age_sex_adjusted", estimable)
ex_152 <- logistic_results |> filter(dataset == "GSE152004", gene_set == reactome,
                                     method == "mean_z_excluding_original32", outcome == "T2_high",
                                     model_label == "unadjusted", estimable)
excl_ok <- nrow(ex_series) && ex_series$mean_difference[1] > 0 &&
  nrow(ex_118) && ex_118$beta_per_sd[1] > 0 &&
  nrow(ex_152) && ex_152$beta_per_sd[1] > 0
ss_series <- gse18965_results |> filter(dataset == "GSE18965_series", gene_set == reactome, method == "ssGSEA_original")
ss_ok <- nrow(ss_series) && ss_series$mean_difference[1] > 0
cel <- gse18965_results |> filter(dataset == "GSE18965_CEL_sensitivity", gene_set == reactome, method == "mean_z_original")
cel_not_opposite <- nrow(cel) && cel$mean_difference[1] >= 0
if (reactome_series_positive && independent_sig && other_consistent && excl_ok && ss_ok && cel_not_opposite) {
  final_class <- "Strong support"
  manuscript_action <- "retain repair-ECM as a main axis"
} else if (reactome_series_positive &&
           ((nrow(pf118) && pf118$effect[1] > 0) || (nrow(pf152) && pf152$effect[1] > 0))) {
  final_class <- "Partial support"
  manuscript_action <- "retain repair-ECM only as exploratory/sensitivity axis"
} else {
  final_class <- "Not supported"
  manuscript_action <- "remove repair-ECM from title and primary framework"
}

classification <- data.frame(
  rule = c("GSE18965 Series Reactome AA higher", "At least one independent primary FDR<0.05 and positive",
           "Both independent primary directions positive", "Excluding original32 overlap remains positive",
           "ssGSEA Series direction positive",
           "CEL sensitivity not opposite", "Final evidence classification", "Consequence for manuscript"),
  value = c(reactome_series_positive, independent_sig, other_consistent, excl_ok, ss_ok, cel_not_opposite,
            final_class, manuscript_action)
)

log_msg("Writing workbooks.")
make_wb <- function(path, sheets) {
  wb <- createWorkbook()
  for (nm in names(sheets)) add_sheet(wb, nm, sheets[[nm]])
  saveWorkbook(wb, path, overwrite = TRUE)
}
make_wb(file.path(CONFIG$out_dir, "02_gene_sets", "external_ECM_gene_coverage_by_cohort.xlsx"), list(coverage = coverage))
make_wb(file.path(CONFIG$out_dir, "02_gene_sets", "sample_level_external_ECM_scores.xlsx"), list(scores = sample_scores))
make_wb(file.path(CONFIG$out_dir, "03_GSE18965_series", "GSE18965_external_ECM_results.xlsx"), list(GSE18965 = gse18965_results, FDR_families = fdr_families))
make_wb(file.path(CONFIG$out_dir, "05_GSE118761_nasal", "GSE118761_external_ECM_results.xlsx"), list(GSE118761 = logistic_results |> filter(grepl("GSE118761", dataset))))
make_wb(file.path(CONFIG$out_dir, "07_GSE152004", "GSE152004_external_ECM_results.xlsx"), list(GSE152004 = logistic_results |> filter(dataset == "GSE152004")))
make_wb(file.path(CONFIG$out_dir, "08_original32_comparison", "external_ECM_vs_original32.xlsx"), list(comparison = external_vs_original32))
make_wb(file.path(CONFIG$out_dir, "09_cross_cohort_meta", "cross_cohort_effect_direction_table.xlsx"), list(direction = direction_table))
make_wb(file.path(CONFIG$out_dir, "09_cross_cohort_meta", "cross_cohort_meta_analysis.xlsx"), list(meta_analysis = meta_analysis))
make_wb(file.path(CONFIG$out_dir, "10_downstream_sensitivity", "downstream_cluster_projection_sensitivity.xlsx"), downstream)
make_wb(file.path(CONFIG$out_dir, "12_report", "External_ECM_multicohort_full_results.xlsx"),
        list(summary = classification, FDR_families = fdr_families, coverage = coverage,
             sample_scores = sample_scores, GSE18965 = gse18965_results,
             logistic_results = logistic_results, original32_comparison = external_vs_original32,
             direction = direction_table, meta_analysis = meta_analysis,
             package_versions = data.frame(package = required, version = vapply(required, function(p) as.character(packageVersion(p)), character(1)))))

copy_map <- list(
  "02_gene_sets/external_ECM_gene_sets_original.xlsx" = "external_ECM_gene_sets_original.xlsx",
  "02_gene_sets/external_ECM_gene_coverage_by_cohort.xlsx" = "external_ECM_gene_coverage_by_cohort.xlsx",
  "02_gene_sets/sample_level_external_ECM_scores.xlsx" = "sample_level_external_ECM_scores.xlsx",
  "03_GSE18965_series/GSE18965_external_ECM_results.xlsx" = "GSE18965_external_ECM_results.xlsx",
  "05_GSE118761_nasal/GSE118761_external_ECM_results.xlsx" = "GSE118761_external_ECM_results.xlsx",
  "07_GSE152004/GSE152004_external_ECM_results.xlsx" = "GSE152004_external_ECM_results.xlsx",
  "08_original32_comparison/external_ECM_vs_original32.xlsx" = "external_ECM_vs_original32.xlsx",
  "09_cross_cohort_meta/cross_cohort_effect_direction_table.xlsx" = "cross_cohort_effect_direction_table.xlsx",
  "09_cross_cohort_meta/cross_cohort_meta_analysis.xlsx" = "cross_cohort_meta_analysis.xlsx",
  "10_downstream_sensitivity/downstream_cluster_projection_sensitivity.xlsx" = "downstream_cluster_projection_sensitivity.xlsx",
  "12_report/External_ECM_multicohort_full_results.xlsx" = "External_ECM_multicohort_full_results.xlsx"
)
for (src in names(copy_map)) file.copy(file.path(CONFIG$out_dir, src), file.path(CONFIG$out_dir, copy_map[[src]]), overwrite = TRUE)

summary_value <- function(df, rows, fields) {
  z <- df[rows, , drop = FALSE]
  if (!nrow(z)) return("not estimable")
  paste(paste(fields, signif(unlist(z[1, fields, drop = FALSE]), 4), sep = "="), collapse = ", ")
}
react_cov_series <- coverage |> filter(cohort == "GSE18965_series", gene_set == reactome, variant == "original")
report <- c(
  "# External ECM Multicohort Analysis Report",
  "",
  "## Fixed Opening Summary",
  "- MSigDB version: 2026.1.Hs; Reactome/Hallmark/KEGG from local v2026.1 GMT, NABA from local MSigDB 2026.1.Hs C2 RDS.",
  "- Primary external set: REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION.",
  "- Confirmatory external set: NABA_CORE_MATRISOME.",
  paste0("- GSE18965 Series score coverage and AA-HN effect: coverage=", signif(react_cov_series$coverage[1], 4), "; ",
         summary_value(gse18965_results, gse18965_results$dataset == "GSE18965_series" & gse18965_results$gene_set == reactome & gse18965_results$method == "mean_z_original", c("mean_difference", "ci_low", "ci_high", "welch_p"))),
  paste0("- GSE18965 CEL sensitivity effect: ",
         summary_value(gse18965_results, gse18965_results$dataset == "GSE18965_CEL_sensitivity" & gse18965_results$gene_set == reactome & gse18965_results$method == "mean_z_original", c("mean_difference", "ci_low", "ci_high", "welch_p"))),
  paste0("- GSE118761 nasal atopy OR per 1 SD: ",
         summary_value(logistic_results, logistic_results$dataset == "GSE118761_nasal" & logistic_results$gene_set == reactome & logistic_results$method == "mean_z_original" & logistic_results$outcome == "atopy" & logistic_results$model_label == "age_sex_adjusted", c("OR_per_sd", "CI_low", "CI_high", "p_value"))),
  paste0("- GSE118761 nasal asthma OR per 1 SD: ",
         summary_value(logistic_results, logistic_results$dataset == "GSE118761_nasal" & logistic_results$gene_set == reactome & logistic_results$method == "mean_z_original" & logistic_results$outcome == "asthma" & logistic_results$model_label == "age_sex_adjusted", c("OR_per_sd", "CI_low", "CI_high", "p_value"))),
  "- GSE118761 tracheal sensitivity: reported separately in GSE118761_external_ECM_results.xlsx; nasal and tracheal were not pooled.",
  paste0("- GSE152004 T2-high OR per 1 SD: ",
         summary_value(logistic_results, logistic_results$dataset == "GSE152004" & logistic_results$gene_set == reactome & logistic_results$method == "mean_z_original" & logistic_results$outcome == "T2_high", c("OR_per_sd", "CI_low", "CI_high", "p_value"))),
  paste0("- GSE152004 asthma OR per 1 SD: ",
         summary_value(logistic_results, logistic_results$dataset == "GSE152004" & logistic_results$gene_set == reactome & logistic_results$method == "mean_z_original" & logistic_results$outcome == "asthma", c("OR_per_sd", "CI_low", "CI_high", "p_value"))),
  paste0("- Primary-family BH-FDR values: ", paste(primary_fam$dataset, primary_fam$endpoint, signif(primary_fam$BH_FDR, 4), sep = "=", collapse = "; ")),
  "- Results after excluding original 32 overlapping genes: complete results are retained with method mean_z_excluding_original32.",
  "- Mean-z versus ssGSEA consistency: complete sensitivity results are retained with method ssGSEA_original.",
  "- Cross-cohort direction consistency: see cross_cohort_effect_direction_table.xlsx.",
  paste0("- Original k=3 clustering robustness: ARI=", signif(cluster_metrics$ARI_vs_original[1], 4), "; NMI=", signif(cluster_metrics$NMI_vs_original[1], 4), "; mean silhouette=", signif(cluster_metrics$mean_silhouette[1], 4), "."),
  paste0("- Projection robustness in GSE118761: ARI=", signif(proj_metrics$ARI_vs_original_projected[1], 4), "; NMI=", signif(proj_metrics$NMI_vs_original_projected[1], 4), "."),
  paste0("- Final evidence classification: ", final_class),
  paste0("- Consequence for manuscript: ", manuscript_action, "."),
  "",
  "## Notes",
  "- GSE152004 has no age, sex, batch, or atopy columns in the provided metadata; age/sex-adjusted and atopy analyses are therefore not estimated.",
  "- GSE118761 technical batch and treatment status are not available; age/sex adjusted models are provided.",
  "- NABA gene sets were read from local MSigDB 2026.1.Hs C2 RDS because c2.cgp/c2.all GMT was not present locally; exact gene set names and members are exported.",
  "- No manuscript, figure, title, or response letter was modified."
)
writeLines(report, file.path(CONFIG$out_dir, "12_report", "External_ECM_multicohort_report.md"))
file.copy(file.path(CONFIG$out_dir, "12_report", "External_ECM_multicohort_report.md"),
          file.path(CONFIG$out_dir, "External_ECM_multicohort_report.md"), overwrite = TRUE)

changelog <- c(
  "# CHANGELOG_external_ECM",
  "",
  "- Used R 4.2.3 with unified library user-managed R library via config/use_local_r_libs.R.",
  "- Installed/loaded required local packages from unified library; downloaded package archives were backed up under user-managed package download cache when provided.",
  "- Generated INPUT_MAPPING.md before modeling and recorded missing GSE152004 age/sex/atopy/batch as not estimable rather than guessed.",
  "- Resolved a single duplicated IGK row in the prior GSE18965 gene-level exports by retaining the highest-mean row; IGK is not part of the locked external ECM sets.",
  "- Locked the six prespecified external ECM gene sets before outcome modeling.",
  "- Used GSE18965 Series Matrix gene-level log2 expression as main analysis and CEL/GC-RMA as sensitivity.",
  "- Processed GSE118761 and GSE152004 raw counts through nonzero filtering, duplicate symbol summing, and DESeq2 VST(blind=TRUE).",
  "- Computed cohort-internal mean-z scores, ssGSEA sensitivity scores, and excluding-original32 scores.",
  "- Retained GSVA as the first ssGSEA backend; when GSVA 1.46.0 failed under the installed matrixStats interface, used a deterministic in-script sample-ranking ssGSEA fallback and logged this compatibility path.",
  "- Ran all prespecified group/logistic/meta/downstream analyses and retained non-significant, opposite-direction, and not-estimable results.",
  "- Did not modify manuscript, figures, title, response letter, original 32-gene list, or external gene sets."
)
writeLines(changelog, file.path(CONFIG$out_dir, "CHANGELOG_external_ECM.md"))

capture.output(sessionInfo(), file = file.path(CONFIG$out_dir, "sessionInfo.txt"))
file.copy("External_ECM_multicohort_analysis.R", file.path(CONFIG$out_dir, "External_ECM_multicohort_analysis.R"), overwrite = TRUE)

log_msg("Completed External ECM multicohort analysis.")






