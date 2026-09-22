source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "00_projection_common.R"))
projection_log("START 06_summarize_projection_phenotypes.R")

d <- projection_read_tsv(file.path(projection_out, "GSE118761_projection_assignments_FINAL.tsv"))
summarize_tissue <- function(tissue) {
  x <- d[d$tissue == tissue, ]
  do.call(rbind, lapply(c("atopy", "asthma", "wheeze"), function(phenotype) {
    y <- as.integer(x[[phenotype]])
    tab <- table(factor(x$assigned_state, levels = c("E1", "E2", "E3")), factor(y, levels = c(0, 1)))
    overall_p <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
    do.call(rbind, lapply(c("E1", "E2", "E3"), function(state) {
      yy <- y[x$assigned_state == state & !is.na(y)]
      positive <- sum(yy == 1L)
      total <- length(yy)
      ci <- if (total) binom.test(positive, total)$conf.int else c(NA_real_, NA_real_)
      data.frame(tissue = tissue, state = state, phenotype = phenotype,
                 n_total = total, n_positive = positive,
                 percentage = if (total) 100 * positive / total else NA_real_,
                 CI_low = 100 * ci[1], CI_high = 100 * ci[2], P_value_overall = overall_p)
    }))
  }))
}
phenotypes <- do.call(rbind, lapply(c("nasal", "tracheal"), summarize_tissue))
projection_write_tsv(phenotypes[phenotypes$tissue == "nasal", ], file.path(projection_out, "GSE118761_nasal_projection_phenotypes_FINAL.tsv"))
projection_write_tsv(phenotypes[phenotypes$tissue == "tracheal", ], file.path(projection_out, "GSE118761_tracheal_projection_phenotypes_FINAL.tsv"))

confidence <- do.call(rbind, lapply(split(d, list(d$tissue, d$assigned_state), drop = TRUE), function(x) {
  data.frame(tissue = x$tissue[1], state = x$assigned_state[1], n = nrow(x),
             median_margin = median(x$margin), IQR_margin = IQR(x$margin),
             low_confidence_n = sum(x$assignment_confidence == "low"),
             low_confidence_percent = 100 * mean(x$assignment_confidence == "low"))
}))
confidence <- confidence[order(match(confidence$tissue, c("nasal", "tracheal")), match(confidence$state, c("E1", "E2", "E3"))), ]
projection_write_tsv(confidence, file.path(projection_out, "GSE118761_projection_confidence_FINAL.tsv"))

nasal <- phenotypes[phenotypes$tissue == "nasal", ]
claim_status <- function(phenotype) {
  z <- nasal[nasal$phenotype == phenotype, ]
  p <- unique(z$P_value_overall)
  e1 <- z$percentage[z$state == "E1"]
  other <- z$percentage[z$state != "E1"]
  if (e1 > max(other) && is.finite(p) && p < 0.05) return("SUPPORTED")
  if (e1 > max(other)) return("SUPPORTED_BUT_WEAKER")
  if (e1 < min(other)) return("CHANGED_DIRECTION")
  "NOT_SUPPORTED"
}
checks <- do.call(rbind, lapply(c("atopy", "asthma", "wheeze"), function(ph) {
  z <- nasal[nasal$phenotype == ph, ]
  data.frame(
    claim = paste0("E1 showed higher ", ph, " frequency after nasal projection"),
    old_support = if (ph == "atopy") "Locked aggregate: E1/E2/E3 76.7/16.7/28.6%" else "Locked aggregate: E1/E2/E3 53.3/33.3/71.4%",
    new_support = paste0("E1/E2/E3 ", paste(sprintf("%.1f%%", z$percentage[match(c("E1", "E2", "E3"), z$state)]), collapse = "/"),
                         "; Fisher P=", signif(unique(z$P_value_overall), 5)),
    status = claim_status(ph)
  )
}))
projection_write_tsv(checks, file.path(projection_audit_dir, "projection_conclusion_check.tsv"))
projection_log("PHENOTYPES recalculated from final sample assignments")
projection_log("END 06_summarize_projection_phenotypes.R")
