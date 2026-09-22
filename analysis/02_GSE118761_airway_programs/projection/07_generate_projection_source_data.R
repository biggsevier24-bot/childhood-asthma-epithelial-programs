source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "00_projection_common.R"))
projection_log("START 07_generate_projection_source_data.R")

assignments <- projection_read_tsv(file.path(projection_out, "GSE118761_projection_assignments_FINAL.tsv"))
counts <- projection_read_tsv(file.path(projection_out, "GSE118761_projection_counts_FINAL.tsv"))
nasal_pheno <- projection_read_tsv(file.path(projection_out, "GSE118761_nasal_projection_phenotypes_FINAL.tsv"))
trach_pheno <- projection_read_tsv(file.path(projection_out, "GSE118761_tracheal_projection_phenotypes_FINAL.tsv"))
phenotypes <- rbind(nasal_pheno, trach_pheno)
confidence <- projection_read_tsv(file.path(projection_out, "GSE118761_projection_confidence_FINAL.tsv"))

figure5f <- nasal_pheno[nasal_pheno$phenotype %in% c("atopy", "asthma"),
                       c("state", "phenotype", "n_positive", "n_total", "percentage", "CI_low", "CI_high", "P_value_overall")]
projection_write_tsv(figure5f, file.path(projection_figure_dir, "Figure5f_source_data.tsv"))

s9 <- merge(counts, confidence, by = c("tissue", "state", "n"), all.x = TRUE, sort = FALSE)
for (ph in c("atopy", "asthma", "wheeze")) {
  z <- phenotypes[phenotypes$phenotype == ph, c("tissue", "state", "n_positive", "percentage")]
  names(z)[3:4] <- paste0(ph, c("_positive_n", "_percentage"))
  s9 <- merge(s9, z, by = c("tissue", "state"), all.x = TRUE, sort = FALSE)
}
s9 <- s9[order(match(s9$tissue, c("nasal", "tracheal")), match(s9$state, c("E1", "E2", "E3"))), ]
projection_write_tsv(s9, file.path(projection_table_dir, "Supplementary_Table_S9_source.tsv"))

s8_audit <- data.frame(
  row = "Supplementary Table S8 external-ECM random-set projection assignment margins",
  depends_on_projection = "No - uses Reactome/NABA external-definition specificity pipeline, not final fixed32 assignments",
  old_value = "Retained external-specificity values",
  new_value = "Unchanged",
  changed_YN = "N",
  source_file = "outputs/05_external_ECM/corrected_projection_specificity.xlsx"
)
projection_write_tsv(s8_audit, file.path(projection_audit_dir, "S8_projection_dependency_audit.tsv"))

methods_lines <- c(
  "GSE118761 projection used the GSE152004 three-state solution as the reference.",
  "For GSE118761, T2, IFN, and fixed 32-gene repair-ECM raw scores were calculated as the mean of gene-wise z-standardized VST expression among available genes, with gene-wise standardization performed across all 104 GSE118761 samples.",
  "Each raw program score was then z-standardized separately within the 55 nasal and 49 tracheal samples.",
  "Reference coordinates were the standardized T2, IFN, and repair-ECM scores of the 695 GSE152004 samples assigned by the locked k=3 solution; each E1, E2, and E3 centroid was the arithmetic mean of these three coordinates within that state.",
  "Nasal and tracheal samples were projected separately by Euclidean distance in the ordered T2/IFN/repair-ECM feature space and assigned to the nearest centroid; exact ties, if any, were resolved in E1, E2, E3 order.",
  "Assignment margin was defined as the second-nearest minus nearest centroid distance, and margins less than or equal to 0.10 were classified as low confidence.",
  "Atopy, asthma, and wheeze summaries were calculated directly from the final sample-level assignments."
)
writeLines(methods_lines, file.path(projection_docs_dir, "GSE118761_projection_METHODS_EXACT.txt"))

fmt <- function(ph) {
  z <- nasal_pheno[nasal_pheno$phenotype == ph, ]
  paste(sprintf("%s %.1f%% (%d/%d)", z$state, z$percentage, z$n_positive, z$n_total), collapse = ", ")
}
result_lines <- c(
  "Using the fully executable nearest-centroid workflow, GSE118761 nasal samples were assigned to E1, E2, and E3 as reported in GSE118761_projection_counts_FINAL.tsv.",
  paste0("Nasal atopy: ", fmt("atopy"), "."),
  paste0("Nasal asthma: ", fmt("asthma"), "."),
  paste0("Nasal wheeze: ", fmt("wheeze"), "."),
  "These projection results are exploratory and do not constitute formal external validation of clinical endotypes."
)
writeLines(result_lines, file.path(projection_docs_dir, "GSE118761_projection_RESULTS_FINAL.txt"))
writeLines(c(
  "Figure 5f | Nearest-centroid projection of the GSE152004 states into GSE118761 nasal samples.",
  "Bars show the proportions with atopy and asthma calculated from the final sample-level assignments.",
  "Projection used tissue-standardized T2, IFN, and fixed 32-gene repair-ECM scores and Euclidean distance to the GSE152004 E1, E2, and E3 centroids.",
  "The projection is interpreted as exploratory cross-cohort recapitulation rather than formal validation of clinical endotypes."
), file.path(projection_docs_dir, "GSE118761_projection_FIGURE5_LEGEND_FINAL.txt"))

checks <- projection_read_tsv(file.path(projection_audit_dir, "projection_conclusion_check.tsv"))
atopy_status <- checks$status[grepl("atopy", checks$claim)]
change_map <- data.frame(
  location = c("Abstract Results", "Results cross-cohort projection", "Discussion", "Figure 5 legend", "Supplementary Table S8", "Supplementary Table S9"),
  current_text = c(
    "T2-high/IFN-low state showed higher atopy frequency after projection into GSE118761.",
    "Superseded aggregate projection summary.",
    "Projection interpreted as partial cross-cohort recapitulation.",
    "Superseded aggregate projection summary.",
    "External Reactome/NABA projection-margin specificity results.",
    "Historical aggregate fixed32 projection profiles."
  ),
  current_value = c("higher atopy frequency", "superseded aggregate projection summary", "partial recapitulation", "superseded aggregate projection summary", "external ECM values", "superseded projection profiles"),
  new_value = c(atopy_status, "See GSE118761_projection_RESULTS_FINAL.txt", atopy_status, "See Figure5f_source_data.tsv", "unchanged", "See Supplementary_Table_S9_source.tsv"),
  needs_change = c(ifelse(atopy_status == "SUPPORTED", "wording review", "YES"), "YES", "wording review", "YES", "NO", "YES"),
  reason = c("Reassess from reproducible sample assignments", "Historical aggregate replaced by executable workflow", "Align interpretation with reproducible result", "Panel source rebuilt", "S8 is not derived from fixed32 final assignments", "S9 rebuilt from final sample assignments"),
  source_output = c("audit/projection_conclusion_check.tsv", "docs/GSE118761_projection_RESULTS_FINAL.txt", "audit/projection_conclusion_check.tsv", "figures/Figure5/Figure5f_source_data.tsv", "audit/S8_projection_dependency_audit.tsv", "tables/Supplementary_Table_S9_source.tsv")
)
projection_write_tsv(change_map, file.path(projection_docs_dir, "GSE118761_projection_manuscript_change_map.tsv"))

provenance <- c(
  "# GSE118761 projection provenance",
  "",
  "1. Inputs: outputs/02_GSE118761/projection/provenance/input_manifest.tsv",
  "2. Expression source: public GSE118761 VST matrix and deposited metadata.",
  "3. Program scores: outputs/02_GSE118761/projection/GSE118761_projection_program_scores.tsv",
  "4. Standardization: outputs/02_GSE118761/projection/projection_scaling_parameters.tsv",
  "5. Reference states: outputs/02_GSE118761/projection/GSE152004_centroid_source_samples.tsv",
  "6. Centroids: outputs/02_GSE118761/projection/GSE152004_reference_centroids.tsv",
  "7. Distances: outputs/02_GSE118761/projection/GSE118761_projection_distances.tsv",
  "8. Assignments: outputs/02_GSE118761/projection/GSE118761_projection_assignments_FINAL.tsv",
  "9. Counts: outputs/02_GSE118761/projection/GSE118761_projection_counts_FINAL.tsv",
  "10. Phenotypes: outputs/02_GSE118761/projection/GSE118761_*_projection_phenotypes_FINAL.tsv",
  "11. Figure 5f: figures/Figure5/Figure5f_source_data.tsv",
  "12. Supplementary S8/S9: audit/S8_projection_dependency_audit.tsv and tables/Supplementary_Table_S9_source.tsv",
  "13. Manuscript-dependent text: docs/GSE118761_projection_manuscript_change_map.tsv",
  "",
  "Superseded aggregate projection values are audit-only and are not read by the formal pipeline."
)
writeLines(provenance, file.path(projection_docs_dir, "GSE118761_projection_PROVENANCE.md"))
projection_log("SOURCE_DATA Figure5f and Supplementary Table S9 rebuilt")
projection_log("END 07_generate_projection_source_data.R")
