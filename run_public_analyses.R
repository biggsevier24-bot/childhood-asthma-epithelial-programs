codex_r_helper <- Sys.getenv("CODEX_R_LIB_HELPER", "")
if (nzchar(codex_r_helper) && file.exists(codex_r_helper)) source(codex_r_helper)
options(stringsAsFactors=FALSE)
repo <- normalizePath(getwd(),winslash="/",mustWork=TRUE)
rscript <- file.path(R.home("bin"),if(.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
rows <- list()
run_step <- function(name, scripts, required) {
  missing <- required[!file.exists(required)]
  if(length(missing)) {
    cat(name,": FAIL: required repository file unavailable\n",sep="")
    rows[[length(rows)+1L]] <<- data.frame(analysis=name,status="FAIL",key_result="required repository file unavailable",expected="See data/README.md",match=FALSE)
    return(invisible(FALSE))
  }
  for(s in scripts) {
    z <- system2(rscript,shQuote(file.path(repo,s)),stdout=TRUE,stderr=TRUE)
    if(!is.null(attr(z,"status")) && attr(z,"status")!=0L) {
      rows[[length(rows)+1L]] <<- data.frame(analysis=name,status="FAIL",key_result=paste(tail(z,3),collapse=" | "),expected="successful execution",match=FALSE)
      return(invisible(FALSE))
    }
  }
  rows[[length(rows)+1L]] <<- data.frame(analysis=name,status="PASS",key_result="completed",expected="successful execution",match=TRUE)
  invisible(TRUE)
}

run_step("GSE18965 exact32",c("analysis/01_GSE18965_repair_ECM/06_run_full_derivation.R","analysis/01_GSE18965_repair_ECM/07_verify_final_module.R"),c("data/GSE18965/GSE18965_series_matrix.txt.gz","data/GSE18965/GPL96-57554.txt"))
gse118_standardized <- file.exists(Sys.getenv("GSE118761_EXPR", file.path(repo,"data/public/GSE118761/GSE118761_VST_gene_expression.csv"))) && file.exists(Sys.getenv("GSE118761_META", file.path(repo,"data/public/GSE118761/GSE118761_metadata.csv")))
gse118_scripts <- c(if(!gse118_standardized) "analysis/02_GSE118761_airway_programs/00_download_GSE118761.R",sprintf("analysis/02_GSE118761_airway_programs/%02d_%s.R",1:6,c("prepare_GSE118761","module_scores","PCA","nasal_tracheal_analysis","program_relationships","quadrant_analysis")))
run_step("GSE118761",gse118_scripts,character())
gse152_standardized <- file.exists(Sys.getenv("GSE152004_EXPR", file.path(repo,"data/public/GSE152004/GSE152004_VST_gene_expression.csv"))) && file.exists(Sys.getenv("GSE152004_META", file.path(repo,"data/public/GSE152004/GSE152004_metadata.csv")))
gse152_scripts <- c(if(!gse152_standardized) "analysis/03_GSE152004_primary_validation/00_download_GSE152004.R",sprintf("analysis/03_GSE152004_primary_validation/%02d_%s.R",1:7,c("prepare_GSE152004","define_T2_high","module_scores","logistic_regression","E1_E2_E3_clustering","multicollinearity_and_residualization_recovered","residualization")))
run_step("GSE152004",gse152_scripts,character())
gse152_raw <- file.path(repo, "data/public/GSE152004/raw/GSE152004_695_raw_counts.txt.gz")
gse152_cache <- Sys.getenv("GSE152004_COUNTS_CACHE", "")
if (file.exists(gse152_raw) || (nzchar(gse152_cache) && file.exists(gse152_cache))) {
  run_step(
    "GSE152004 score construction sensitivity",
    c(
      "analysis/03_GSE152004_primary_validation/00_download_raw_counts.R",
      "analysis/03_GSE152004_primary_validation/08_score_construction_sensitivity.R"
    ),
    character()
  )
} else {
  verified <- file.path(repo, "outputs/03_GSE152004/GSE152004_alternative_repair_verified_output.csv")
  rows[[length(rows) + 1L]] <- data.frame(
    analysis = "GSE152004 score construction sensitivity",
    status = "VERIFIED_OUTPUT_INPUT_NOT_REDISTRIBUTED",
    key_result = if (file.exists(verified)) "retained machine-readable output verified" else "verified output missing",
    expected = "secondary sensitivity; official raw-count URL and download script provided",
    match = file.exists(verified)
  )
}
run_step("GSE118761 projection",sprintf("analysis/02_GSE118761_airway_programs/projection/%02d_%s.R",1:8,c("prepare_projection_inputs","calculate_projection_features","build_reference_centroids","assign_nearest_centroid","summarize_projection_counts","summarize_projection_phenotypes","generate_projection_source_data","verify_projection_results")),c("config/GSE118761_projection_final.yaml"))
run_step("GSE118761 projection deliverables",c("figures/Figure5/make_Figure5f.R","tables/make_Supplementary_Table_S9.R"),c("figures/Figure5/Figure5f_source_data.tsv","tables/Supplementary_Table_S9_source.tsv"))
run_step("External ECM",sprintf("analysis/05_external_ECM_sensitivity/%02d_%s.R",1:5,c("Reactome_ECM","NABA_core_matrisome","overlap_exclusion","random_gene_set_audit","ssGSEA_sensitivity")),c("data/reference_gene_sets/Reactome_ECM_manifest.csv","data/reference_gene_sets/NABA_core_matrisome_manifest.csv"))
rows[[length(rows)+1L]] <- data.frame(analysis="Hospital cohort",status="CODE_ONLY",key_result="protected patient-level data not public",expected="code and schema only",match=TRUE)
out <- do.call(rbind,rows)
dir.create(file.path(repo,"outputs"),showWarnings=FALSE)
write.table(out,file.path(repo,"outputs/full_repository_run_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
print(out)
if(any(out$status=="FAIL")) quit(status=1)
