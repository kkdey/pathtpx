

################  Simulation of pathtpx models  ########################

G= 1000
N = 100
L = 40

pathways <- list()
for(l in 1:L){
  pathways[[l]] <- sample(1:G, 20, replace = FALSE)
}

pathways_mat <- matrix(0, G, length(pathways))
for(l in 1:length(pathways)){
  pathways_mat[pathways[[l]], l] <- rep(1, length(pathways[[l]])) 
}

library(slam)
library(maptpx)

pathnormalizetpx <- function(x, byrow=TRUE){
  if(byrow){ s <- row_sums(x)
  s[s==0] <- 1
  return( x/s ) }
  else{
    s <- col_sums(x)
    s[s==0] <- 1
    return(t(t(x)/s)) }
}

freq <- t(pathnormalizetpx(pathways_mat, byrow = FALSE))

omega_sim <- matrix(0, N, L)
for(n in 1:N){
  num <- floor(sample(1:30))
  omega_sim[n, sample(L, num, replace = FALSE)] <- 1
}

omega_sim <- pathnormalizetpx(omega_sim, byrow = TRUE)

counts <- t( do.call(cbind,
                     lapply(1:dim(omega_sim)[1], 
                            function(x) 
                              rmultinom(1,10000,prob=omega_sim[x,]%*%freq))))


topic_clus <- pathtopics(counts, pathways, ord=FALSE, tol = 0.1, tmax=100)

omega_out <- topic_clus$omega
theta_out <- topic_clus$theta

plot(theta_out[,5])
plot(freq[5,])
plot(omega_out[20,], col="blue", pch=20)
points(omega_sim[20,], col="red", pch=20)

plot(omega_sim[90,], col="blue", pch=20)
points(omega_out[90,], col="red", pch=20)

