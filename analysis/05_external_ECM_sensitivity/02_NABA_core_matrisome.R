repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "../.."), winslash="/", mustWork=TRUE)
cat("NABA Core Matrisome manifest: ",file.path(repo,"data/reference_gene_sets/NABA_core_matrisome_manifest.csv"),"\n",sep="")
