options(stringsAsFactors=FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "../.."), winslash="/", mustWork=TRUE)
fixed <- toupper(readLines(file.path(repo,"data/reference_gene_sets/repair_ECM_fixed32.txt")))
for(f in c("Reactome_ECM_manifest.csv","NABA_core_matrisome_manifest.csv")){
  x<-read.csv(file.path(repo,"data/reference_gene_sets",f),check.names=FALSE); genes<-unique(toupper(unlist(x))); genes<-genes[nzchar(genes)&!is.na(genes)]
  cat(f,": total=",length(genes)," overlap=",length(intersect(genes,fixed))," excluded=",length(setdiff(genes,fixed)),"\n",sep="")
}
