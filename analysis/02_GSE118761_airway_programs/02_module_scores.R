options(stringsAsFactors = FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "../.."), winslash="/", mustWork=TRUE)
source(file.path(repo, "R/common/module_scoring.R"))
a <- readRDS(file.path(repo, "outputs/02_GSE118761/GSE118761_aligned_input.rds"))
sets <- list(T2=read_gene_set(file.path(repo,"data/reference_gene_sets/T2_fixed3.txt")), IFN=read_gene_set(file.path(repo,"data/reference_gene_sets/IFN_fixed50.txt")), repair_ECM=read_gene_set(file.path(repo,"data/reference_gene_sets/repair_ECM_fixed32.txt")))
out <- a$metadata
for (nm in names(sets)) out[[paste0(nm,"_z")]] <- score_mean_gene_z(a$expr, sets[[nm]])
write.csv(out, file.path(repo,"outputs/02_GSE118761/GSE118761_recomputed_module_scores.csv"), row.names=FALSE)
formal <- out
names(formal)[names(formal) == "repair_ECM_z"] <- "Repair_ECM_32_score"
names(formal)[names(formal) == "IFN_z"] <- "IFN_score"
names(formal)[names(formal) == "T2_z"] <- "T2_score"
write.csv(formal, file.path(repo,"outputs/02_GSE118761/GSE118761_all_sample_module_scores.csv"), row.names=FALSE)
