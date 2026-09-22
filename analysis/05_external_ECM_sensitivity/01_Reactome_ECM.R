repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "../.."), winslash="/", mustWork=TRUE)
cat("Reactome ECM Organization manifest: ",file.path(repo,"data/reference_gene_sets/Reactome_ECM_manifest.csv"),"\n",sep="")
