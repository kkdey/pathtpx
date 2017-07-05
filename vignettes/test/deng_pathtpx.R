
####################  deng et al pathway analysis  #######################

library(singleCellRNASeqMouseDeng2014)
deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

pathways <- read.delim("CPDB_pathways_genes.tab")
pathways1 <- pathways[which(pathways$source == "KEGG"),]


pathways1$hgnc_symbol_ids

pathway_names <- list()
for (l in 1:dim(pathways1)[1]){
  temp <- strsplit(as.character(pathways1[l,4]), "[,]")[[1]]
  index1 <- match(temp, as.character(deng.gene_names))
  index2 <- index1[!is.na(index1)]
  pathway_names[[l]] <- as.character(deng.gene_names)[index2]
}

counts <- t(deng.counts)

topic_clus <- pathtopics(counts[,1:5000], pathway_names, ord=FALSE, tol = 10)
