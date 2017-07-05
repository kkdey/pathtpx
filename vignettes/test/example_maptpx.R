

###############  maptpx solutions ########################


library(devtools)
library(CountClust)
library(singleCellRNASeqMouseDeng2014)
library(GTExV6Brain)
gtex.counts <- Biobase::exprs(GTExV6Brain)
gtex.meta_data <- Biobase::pData(GTExV6Brain)
gtex.gene_names <- rownames(gtex.counts)

library(maptpx)

out <- maptpx::topics(t(gtex.counts),  K=4, tol=0.1)

library(CountClust)
FitGoM(t(gtex.counts),
       K=4, tol=0.1,
       path_rda="GTExV6Brain.FitGoM.rda")

deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

fitgom_gtex <- get(load("../../CountClust/data/GTExV6Brain.FitGoM.rda"))
fitgom_deng <- get(load("../../CountClust/data/MouseDeng2014.FitGoM.rda"))

theta <- as.matrix(fitgom_deng$clust_3$theta)
plot(theta[,2])

mu <- log(theta[,2]) - log(theta[which.min(theta[,2]),2])
mu1 <- digamma(1000*theta[,2]) - digamma(1000*theta[which.min(theta[,2]),2])
sigma <- trigamma(1000*theta[,2]) + trigamma(1000*theta[which.min(theta[,2]),2])

out <- ashr::ash(mu1, sqrt(sigma), mixcompdist = "normal", mode="estimate")
res <- out$result$PosteriorMean
theta_ash <- exp(out$result$PosteriorMean)*theta[which.min(theta[,2]),2]
