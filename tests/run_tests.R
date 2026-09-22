codex_r_helper <- Sys.getenv("CODEX_R_LIB_HELPER", "")
if (nzchar(codex_r_helper) && file.exists(codex_r_helper)) source(codex_r_helper)
options(stringsAsFactors=FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])),".."),winslash="/",mustWork=TRUE)
rscript <- file.path(R.home("bin"),if(.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
tests <- c(
  "test_GSE18965_exact32.R",
  "test_no_hardcoding.R",
  "test_config_code_consistency.R",
  "test_external_ECM_not_used_for_derivation.R",
  "test_T2_high_definition.R",
  "test_module_scores.R",
  "test_public_input_rebuild.R",
  "test_key_manuscript_results.R",
  "test_repository_structure.R",
  "test_GSE118761_projection.R",
  "test_GSE118761_sample_counts.R",
  "test_GSE118761_Figure4_PCA.R",
  "test_archive_independence.R",
  "test_security_scan.R"
)
rows <- lapply(tests,function(t){z<-system2(rscript,shQuote(file.path(repo,"tests",t)),stdout=TRUE,stderr=TRUE); ok<-is.null(attr(z,"status"))||attr(z,"status")==0L; data.frame(test=t,PASS_FAIL=if(ok)"PASS" else "FAIL",details=paste(tail(z,4),collapse=" | "),stringsAsFactors=FALSE)})
out<-do.call(rbind,rows); write.table(out,file.path(repo,"tests/test_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE); print(out)
if(any(out$PASS_FAIL=="FAIL")) quit(status=1)
