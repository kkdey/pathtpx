
digamma(10^(-7)) - digamma(10^(-6))

sigma <- trigamma(10^(-7)) + trigamma(10^(-6))

library(devtools)
library(singleCellRNASeqMouseDeng2014)
library(GTExV6Brain)
gtex.counts <- Biobase::exprs(GTExV6Brain)
gtex.meta_data <- Biobase::pData(GTExV6Brain)
gtex.gene_names <- rownames(gtex.counts)

library(maptpx)

out <- maptpx::topics(t(gtex.counts),  K=3, tol=100)
