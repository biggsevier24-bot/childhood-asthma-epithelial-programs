options(stringsAsFactors = FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "../.."), winslash = "/", mustWork = TRUE)
read_csv <- function(path) read.csv(file.path(repo, path), check.names = FALSE)

verification <- read.delim(file.path(repo, "outputs/verification.tsv"), check.names = FALSE)
scores118 <- read_csv("outputs/02_GSE118761/GSE118761_all_sample_module_scores.csv")
quadrants <- read_csv("outputs/02_GSE118761/GSE118761_quadrant_summary.csv")
auc <- read_csv("outputs/02_GSE118761/GSE118761_wheezer_atopy_CV_ROC_summary.csv")
model152 <- read_csv("outputs/03_GSE152004/GSE152004_molecular_T2_high_logistic_model.csv")
states152 <- read_csv("outputs/03_GSE152004/GSE152004_final_state_profiles.csv")
s13 <- read_csv("outputs/03_GSE152004/Supplementary_Table_S13_multicollinearity_residualization.csv")
alt152 <- read_csv("outputs/03_GSE152004/GSE152004_alternative_repair_recomputed.csv")

q_counts <- setNames(quadrants$n, quadrants$quadrant)
s152 <- s13[s13$cohort == "GSE152004", ][1, ]
rows <- data.frame(
  result = c("GSE18965_exact32", "GSE118761_n", "GSE118761_nasal_n", "GSE118761_tracheal_n", "GSE118761_quadrants",
             "GSE118761_pooled_OOF_AUC", "GSE118761_projection", "GSE152004_n", "GSE152004_T2_split", "GSE152004_repair_OR",
             "GSE152004_IFN_OR", "GSE152004_k3", "GSE152004_repair_T2", "GSE152004_repair_IFN", "GSE152004_VIF_range",
             "GSE152004_condition_number", "GSE152004_repair_residual_variance", "GSE152004_alternative_repair_OR"),
  observed = c(as.character(verification$EXACT_MATCH_32), nrow(scores118), sum(scores118$tissue == "nasal"), sum(scores118$tissue == "tracheal"),
               paste(sort(q_counts), collapse = "/"), signif(auc$pooled_cv_auc, 12), "UNRESOLVED",
               unique(model152$n), paste(unique(model152$events), unique(model152$n - model152$events), sep = "/"),
               signif(model152$OR[model152$term == "Repair_z"], 12), signif(model152$OR[model152$term == "IFN_z"], 12),
               paste(sort(states152$n), collapse = "/"), signif(s152$Pearson_T2_repair, 12), signif(s152$Pearson_IFN_repair, 12),
               paste(signif(min(s152[c("VIF_T2", "VIF_IFN", "VIF_repair_ECM")]), 8), signif(max(s152[c("VIF_T2", "VIF_IFN", "VIF_repair_ECM")]), 8), sep = "-"),
               signif(s152$condition_number, 12), signif(s152$residual_variance_retained, 12), signif(alt152$OR, 12)),
  expected = c("TRUE", "104", "55", "49", "12/12/15/16", "0.719780219780", "30/18/7;15/29/5", "695", "348/347",
               "0.542309 approximately", "0.770618 approximately", "58/299/338", "-0.292 approximately", "0.206 approximately",
               "1.058-1.124 approximately", "1.434 approximately", "0.889 approximately", "1.16 approximately"),
  status = c("PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "FAIL_UNRESOLVED", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS")
)
write.table(rows, file.path(repo, "audit/final_full_clean_rerun/key_result_verification.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

html <- c("<!doctype html><html><head><meta charset='utf-8'><title>Final repository audit</title><style>body{font-family:Arial,sans-serif;max-width:1000px;margin:32px auto;line-height:1.45}table{border-collapse:collapse;width:100%}th,td{border:1px solid #bbb;padding:6px;text-align:left}.pass{color:#176b2c}.fail{color:#a21b1b}</style></head><body>",
  "<h1>Full manuscript GitHub final release audit</h1>",
  "<p><strong>PRELOCK DEVELOPMENT REPORT</strong><br>Superseded projection audit retained for history only.</p>",
  "<h2>Key-result gate</h2><table><tr><th>Result</th><th>Observed</th><th>Expected</th><th>Status</th></tr>",
  apply(rows, 1, function(x) sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>", x[[1]], x[[2]], x[[3]], x[[4]])),
  "</table><h2>Blocking issue</h2><p>The executable method that produced the manuscript GSE118761 projection counts was not recovered. No counts or assignments were hardcoded. See <code>audit/GSE118761_projection_unresolved_report.md</code>.</p></body></html>")
writeLines(html, file.path(repo, "FULL_MANUSCRIPT_GITHUB_FINAL_RELEASE_REPORT.html"), useBytes = TRUE)
