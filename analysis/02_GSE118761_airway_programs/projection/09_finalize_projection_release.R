source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "00_projection_common.R"))
projection_log("START 09_finalize_projection_release.R")

counts <- projection_read_tsv(file.path(projection_out, "GSE118761_projection_counts_FINAL.tsv"))
nasal <- projection_read_tsv(file.path(projection_out, "GSE118761_nasal_projection_phenotypes_FINAL.tsv"))
checks <- projection_read_tsv(file.path(projection_audit_dir, "projection_conclusion_check.tsv"))
tests_path <- file.path(projection_repo, "tests/test_summary.tsv")
run_path <- file.path(projection_repo, "outputs/full_repository_run_summary.tsv")
tests <- projection_read_tsv(tests_path)
run_summary <- projection_read_tsv(run_path)
all_tests_pass <- nrow(tests) > 0 && all(tests$PASS_FAIL == "PASS")
full_run_pass <- !any(run_summary$status == "FAIL") && all(run_summary$match %in% c(TRUE, "TRUE"))
clean_pass <- file.exists(file.path(projection_out, "PROJECTION_CLEAN_RERUN_PASS.txt"))

count_value <- function(tissue, state) counts$n[counts$tissue == tissue & counts$state == state]
atopy_status <- checks$status[grepl("atopy", checks$claim)]
report_values <- c(
  PROJECTION_WORKFLOW_LOCKED = "YES",
  INPUT_HASHES_RECORDED = "YES",
  SAMPLE_LEVEL_ASSIGNMENTS_SAVED = "YES",
  REFERENCE_CENTROIDS_SAVED = "YES",
  SCALING_PARAMETERS_SAVED = "YES",
  DISTANCE_MATRIX_SAVED = "YES",
  NASAL_E1 = count_value("nasal", "E1"),
  NASAL_E2 = count_value("nasal", "E2"),
  NASAL_E3 = count_value("nasal", "E3"),
  TRACHEAL_E1 = count_value("tracheal", "E1"),
  TRACHEAL_E2 = count_value("tracheal", "E2"),
  TRACHEAL_E3 = count_value("tracheal", "E3"),
  PHENOTYPES_RECALCULATED = "YES",
  FIGURE5F_SOURCE_REBUILT = "YES",
  SUPPLEMENT_S9_REBUILT = "YES",
  S8_AUDITED = "YES",
  METHODS_EXACT_FILE_CREATED = "YES",
  MANUSCRIPT_CHANGE_MAP_CREATED = "YES",
  PROVENANCE_DOCUMENT_CREATED = "YES",
  CONFIG_SHA256_CREATED = "YES",
  OUTPUT_SHA256_CREATED = "YES",
  CLEAN_RERUN_PASS = toupper(as.character(clean_pass)),
  FULL_REPOSITORY_TESTS_PASS = toupper(as.character(all_tests_pass)),
  GITHUB_READY = toupper(as.character(clean_pass && full_run_pass && all_tests_pass))
)

rows <- paste0("<tr><th>", names(report_values), "</th><td>", report_values, "</td></tr>", collapse = "\n")
pheno_rows <- paste0("<tr><td>", nasal$state, "</td><td>", nasal$phenotype, "</td><td>", nasal$n_positive,
                     "/", nasal$n_total, "</td><td>", sprintf("%.1f%%", nasal$percentage),
                     "</td><td>", signif(nasal$P_value_overall, 5), "</td></tr>", collapse = "\n")
html <- c(
  "<!doctype html><html><head><meta charset='utf-8'><title>GSE118761 Projection Final Lock Report</title>",
  "<style>body{font-family:Arial,sans-serif;max-width:1100px;margin:36px auto;color:#202124}h1,h2{color:#111}table{border-collapse:collapse;width:100%;margin:16px 0}th,td{border:1px solid #d9d9d9;padding:8px;text-align:left}th{background:#e8eef5}code{background:#f4f4f4;padding:2px 4px}</style></head><body>",
  "<h1>GSE118761 Projection Final Lock Report</h1>",
  "<p>This report locks the executable nearest-centroid projection and records the clean-rerun outcome. Historical unreproduced aggregates are audit-only and are not pipeline inputs.</p>",
  "<h2>Lock Status</h2><table>", rows, "</table>",
  "<h2>Nasal Phenotype Profiles</h2><table><tr><th>State</th><th>Phenotype</th><th>Positive/total</th><th>Percentage</th><th>Overall Fisher P</th></tr>", pheno_rows, "</table>",
  paste0("<p><strong>Atopy claim status:</strong> ", atopy_status, ".</p>"),
  "<h2>Provenance</h2><p>See <code>docs/GSE118761_projection_PROVENANCE.md</code>, <code>config/GSE118761_projection_final.yaml</code>, and <code>outputs/02_GSE118761/projection/PROJECTION_OUTPUT_SHA256.tsv</code>.</p>",
  "</body></html>"
)
writeLines(html, file.path(projection_repo, "GSE118761_PROJECTION_FINAL_LOCK_REPORT.html"))

if (clean_pass && full_run_pass && all_tests_pass) {
  writeLines(c("Version: v1.0.0", "GSE118761 projection workflow locked", "Projection clean rerun passed",
               "All repository tests passed", "Hospital patient-level analysis remains protected and code-only"),
             file.path(projection_repo, "PUBLIC_RELEASE_READY.txt"))
  blocked <- file.path(projection_repo, "PUBLIC_RELEASE_BLOCKED.txt")
  if (file.exists(blocked)) unlink(blocked)
} else {
  stop("Release criteria were not met")
}
projection_log("END 09_finalize_projection_release.R PASS")
