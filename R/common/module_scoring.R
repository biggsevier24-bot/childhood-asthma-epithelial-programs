options(stringsAsFactors = FALSE)

read_gene_set <- function(path) {
  x <- readLines(path, warn = FALSE)
  unique(toupper(trimws(x[nzchar(trimws(x))])))
}

read_gene_expression_csv <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  genes <- toupper(trimws(as.character(x[[1]])))
  x[[1]] <- NULL
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  rownames(m) <- genes
  if (anyDuplicated(genes)) {
    mu <- rowMeans(m, na.rm = TRUE)
    keep <- unlist(tapply(seq_along(genes), genes, function(i) i[which.max(mu[i])]), use.names = FALSE)
    m <- m[sort(keep), , drop = FALSE]
  }
  m
}

align_expression_metadata <- function(expr, metadata) {
  stopifnot("sample_id" %in% names(metadata))
  common <- intersect(colnames(expr), as.character(metadata$sample_id))
  if (!length(common)) stop("No matching sample IDs between expression and metadata")
  expr <- expr[, common, drop = FALSE]
  metadata <- metadata[match(common, metadata$sample_id), , drop = FALSE]
  stopifnot(identical(colnames(expr), as.character(metadata$sample_id)))
  list(expr = expr, metadata = metadata)
}

score_mean_gene_z <- function(expr, genes, min_present = 1L) {
  present <- intersect(unique(toupper(genes)), rownames(expr))
  if (length(present) < min_present) stop("Insufficient genes present for score")
  m <- expr[present, , drop = FALSE]
  s <- apply(m, 1, sd, na.rm = TRUE)
  m <- m[is.finite(s) & s > 0, , drop = FALSE]
  z <- t(scale(t(m)))
  out <- colMeans(z, na.rm = TRUE)
  attr(out, "genes_present") <- rownames(z)
  out
}

score_mean_expression_per_sd <- function(expr, genes) {
  present <- intersect(unique(toupper(genes)), rownames(expr))
  if (!length(present)) stop("No genes present for alternative score")
  as.numeric(scale(colMeans(expr[present, , drop = FALSE], na.rm = TRUE)))
}

score_official_ssgsea <- function(expr, genes, min_size = 1L) {
  if (!requireNamespace("GSVA", quietly = TRUE)) stop("Official GSVA package is required")
  present <- intersect(unique(toupper(genes)), rownames(expr))
  if (length(present) < min_size) stop("Insufficient genes present for ssGSEA")
  sets <- list(program = present)
  exports <- getNamespaceExports("GSVA")
  if ("ssgseaParam" %in% exports) {
    param <- GSVA::ssgseaParam(exprData = expr, geneSets = sets, minSize = min_size,
                               maxSize = Inf, normalize = TRUE)
    return(as.numeric(GSVA::gsva(param, verbose = FALSE)[1, ]))
  }
  as.numeric(GSVA::gsva(expr, sets, method = "ssgsea", min.sz = min_size,
                        max.sz = Inf, ssgsea.norm = TRUE, verbose = FALSE)[1, ])
}

safe_z <- function(x) as.numeric(scale(as.numeric(x)))

binary_logistic_term <- function(data, formula, term) {
  fit <- glm(formula, data = data, family = binomial())
  b <- coef(summary(fit))[term, ]
  ci <- confint.default(fit, term)
  data.frame(term = term, beta = unname(b["Estimate"]), OR = exp(unname(b["Estimate"])),
             CI_low = exp(ci[1]), CI_high = exp(ci[2]), p_value = unname(b["Pr(>|z|)"]),
             n = nobs(fit), events = sum(model.response(model.frame(fit))), stringsAsFactors = FALSE)
}
