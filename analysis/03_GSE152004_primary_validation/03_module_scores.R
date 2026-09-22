options(stringsAsFactors = FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "../.."), winslash="/", mustWork=TRUE)
source(file.path(repo,"R/common/module_scoring.R"))
a <- readRDS(file.path(repo,"outputs/03_GSE152004/GSE152004_aligned_input.rds"))
t2def <- read.csv(file.path(repo,"outputs/03_GSE152004/GSE152004_T2_definition_recomputed.csv"))
sets <- list(T2=read_gene_set(file.path(repo,"data/reference_gene_sets/T2_fixed3.txt")), IFN=read_gene_set(file.path(repo,"data/reference_gene_sets/IFN_fixed50.txt")), repair_ECM=read_gene_set(file.path(repo,"data/reference_gene_sets/repair_ECM_fixed32.txt")))
out <- a$metadata
out$T2_z <- safe_z(score_mean_gene_z(a$expr,sets$T2))
out$IFN_z <- safe_z(score_mean_gene_z(a$expr,sets$IFN))
out$repair_ECM_z <- safe_z(score_mean_gene_z(a$expr,sets$repair_ECM))
out$T2_high <- t2def$T2_high[match(out$sample_id,t2def$sample_id)]
write.csv(out,file.path(repo,"outputs/03_GSE152004/GSE152004_module_scores_recomputed.csv"),row.names=FALSE)
