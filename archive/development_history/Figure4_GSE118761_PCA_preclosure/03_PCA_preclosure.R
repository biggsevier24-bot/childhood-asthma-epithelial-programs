options(stringsAsFactors = FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "../.."), winslash="/", mustWork=TRUE)
a <- readRDS(file.path(repo,"outputs/02_GSE118761/GSE118761_aligned_input.rds"))
for (tissue in c("nasal","tracheal")) {
  idx <- which(tolower(a$metadata$tissue)==tissue)
  p <- prcomp(t(a$expr[,idx,drop=FALSE]), center=TRUE, scale.=FALSE)
  scores <- data.frame(sample_id=a$metadata$sample_id[idx], tissue=tissue, p$x[,1:min(10,ncol(p$x)),drop=FALSE], check.names=FALSE)
  write.csv(scores, file.path(repo,paste0("outputs/02_GSE118761/GSE118761_",tissue,"_PCA_scores_recomputed.csv")), row.names=FALSE)
  modules <- read.csv(file.path(repo,"outputs/02_GSE118761/GSE118761_recomputed_module_scores.csv"),check.names=FALSE)
  modules <- modules[match(scores$sample_id,modules$sample_id),]
  variables <- c(repair_ECM="repair_ECM_z",atopy="atopy",asthma="asthma")
  rows <- do.call(rbind,lapply(seq_len(min(10,ncol(p$x))),function(j) do.call(rbind,lapply(names(variables),function(v){
    data.frame(PC=paste0("PC",j),variable=v,rho=cor(p$x[,j],modules[[variables[[v]]]],method="spearman",use="pairwise.complete.obs"))
  }))))
  write.csv(rows,file.path(repo,paste0("outputs/02_GSE118761/GSE118761_",tissue,"_PCA_correlations.csv")),row.names=FALSE)
}
