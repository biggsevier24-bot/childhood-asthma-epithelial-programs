options(stringsAsFactors=FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])),".."),winslash="/",mustWork=TRUE)
m <- read.csv(file.path(repo,"outputs/03_GSE152004/GSE152004_molecular_T2_high_logistic_model.csv"),check.names=FALSE)
stopifnot(unique(m$n)==695L,unique(m$events)==348L)
stopifnot(abs(m$OR[m$term=="Repair_z"]-0.5423094732)<1e-6,abs(m$OR[m$term=="IFN_z"]-0.7706175706)<1e-6)
c3 <- read.csv(file.path(repo,"outputs/03_GSE152004/GSE152004_final_state_profiles.csv"),check.names=FALSE)
stopifnot(identical(sort(c3$n),sort(c(338L,299L,58L))))
alt_path <- file.path(repo,"outputs/03_GSE152004/GSE152004_alternative_repair_recomputed.csv")
stopifnot(file.exists(alt_path))
alt <- read.csv(alt_path,check.names=FALSE)
stopifnot(
  nrow(alt)==1L,
  identical(alt$branch,"canonical_public_rebuild"),
  alt$N==695L, alt$events==348L,
  alt$repair_genes_expected==32L, alt$repair_genes_used==32L,
  abs(alt$repair_score_SD-0.124834159728788)<1e-10,
  abs(alt$beta-0.13923965972929803)<1e-10,
  abs(alt$SE-0.0781492347094933)<1e-10,
  abs(alt$OR-1.149399531777976)<1e-10,
  abs(alt$CI_low-0.986166548142085)<1e-10,
  abs(alt$CI_high-1.339651285211802)<1e-10,
  abs(alt$p_value-0.0747957278590982)<1e-10,
  abs(alt$IFN_beta-(-0.36197941614975432))<1e-10,
  abs(alt$IFN_SE-0.0851713005946017)<1e-10,
  abs(alt$IFN_OR-0.696296700161266)<1e-10,
  abs(alt$IFN_CI_low-0.589245662743692)<1e-10,
  abs(alt$IFN_CI_high-0.822796204214671)<1e-10,
  abs(alt$IFN_p_value-2.13754956058918e-05)<1e-10
)
s118 <- read.csv(file.path(repo,"outputs/02_GSE118761/GSE118761_all_sample_module_scores.csv"),check.names=FALSE)
stopifnot(nrow(s118)==104L,sum(tolower(s118$tissue)=="nasal")==55L,sum(tolower(s118$tissue)=="tracheal")==49L)
cat("PASS: reproducible key manuscript values; projection is verified separately from sample-level distances\n")
