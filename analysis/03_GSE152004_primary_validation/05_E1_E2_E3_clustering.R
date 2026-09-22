options(stringsAsFactors = FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "../.."), winslash="/", mustWork=TRUE)
set.seed(20260626)
d <- read.csv(file.path(repo,"outputs/03_GSE152004/GSE152004_module_scores_recomputed.csv"),check.names=FALSE)
x <- scale(d[,c("T2_z","IFN_z","repair_ECM_z")])
km <- kmeans(x,centers=3,nstart=500)
cent <- aggregate(x,list(cluster=km$cluster),mean)
e1 <- cent$cluster[which.max(cent$T2_z)]
e3 <- cent$cluster[which.max(cent$IFN_z+cent$repair_ECM_z)]
e2 <- setdiff(cent$cluster,c(e1,e3))[1]
map <- setNames(c("E1_T2_high_IFN_low","E2_T2_low","E3_IFN_high_repair_high"),c(e1,e2,e3))
out <- data.frame(sample_id=d$sample_id,state_label=unname(map[as.character(km$cluster)]),d[,c("T2_z","IFN_z","repair_ECM_z")])
write.csv(out,file.path(repo,"outputs/03_GSE152004/GSE152004_k3_assignments_recomputed.csv"),row.names=FALSE)
counts <- as.data.frame(table(out$state_label)); names(counts) <- c("state_label", "n")
write.csv(counts,file.path(repo,"outputs/03_GSE152004/GSE152004_k3_counts_recomputed.csv"),row.names=FALSE)
write.csv(out,file.path(repo,"outputs/03_GSE152004/GSE152004_final_state_assignments.csv"),row.names=FALSE)
write.csv(counts,file.path(repo,"outputs/03_GSE152004/GSE152004_final_state_profiles.csv"),row.names=FALSE)
