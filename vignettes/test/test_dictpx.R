
#######################  test dictpx  #######################################

library(singleCellRNASeqMouseDeng2014)
deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

deng.counts_red <- deng.counts[1:500,]

counts <- t(deng.counts_red)
X <- CheckCounts(counts)
p <- ncol(X)

pathway <- list()
pathway[[1]] <- 1:10
pathway[[2]] <- 50:60
pathway[[3]] <- 300:330
pathway[[4]] <- 90:120
pathway[[5]] <- 200:250
pathway[[6]] <- 40:70

pathway_mat <- matrix(0, dim(counts)[2], length(pathway))
for(l in 1:length(pathway)){
  pathway_mat[pathway[[l]], l] <- rep(1, length(pathway[[l]])) 
}



initopics <- ashnormalizetpx(pathway_mat, byrow = FALSE)


sample_init=TRUE
shape=NULL
initopics=NULL
tol=0.1
bf=FALSE
kill=2
ord=TRUE
verb=1
admix=TRUE
nbundles=1
use_squarem=FALSE
init.adapt=FALSE
type="full"
ind_model_indices = NULL
signatures=NULL
light=1
method_admix=1
sample_init=TRUE

theta <- initopics
alpha <- shape
K <- ncol(theta)
n <- nrow(X)
p <- ncol(X)
m <- row_sums(X)
if(is.null(alpha)){ alpha <- 1/(K*p) }
if(is.matrix(alpha)){ if(nrow(alpha)!=p || ncol(alpha)!=K){ stop("bad matrix alpha dimensions") }}

## recycle these in tpcweights to save time
xvo <- X$v[order(X$i)]
wrd <- X$j[order(X$i)]-1
doc <- c(0,cumsum(as.double(table(factor(X$i, levels=c(1:nrow(X)))))))

theta2 <- theta
theta2[theta2 == 1] <- 1 - 1e-14
theta2[theta2 == 0] <- 1e-14

system.time(omega <- tpxweights(n=n, p=p, xvo=xvo, wrd=wrd, doc=doc, start=tpxOmegaStart(X, theta2), theta=theta2))

iter <- 0
dif <- tol+1+qn
update <- TRUE
if(verb>0){
  cat("log posterior increase: " )
  digits <- max(1, -floor(log(tol, base=10))) }

Y <- NULL # only used for qn > 0
Q0 <- col_sums(X)/sum(X)
L <- tpxlpost(X=X, theta=theta2, omega=omega, alpha=alpha, admix=admix, grp=grp)
if(is.infinite(L)){ L <- sum( (log(Q0)*col_sums(X))[Q0>0] ) }

iter <- 1;
wtol=10^{-4}
qn=100
nonzero=FALSE
dcut=-10
top_genes=150
burn_in=5
while( update  && iter < tmax ){
  if(admix && wtol > 0 && (iter-1)%%nbundles==0){
    Wfit <- tpxweights(n=nrow(X), p=ncol(X), xvo=xvo, wrd=wrd, doc=doc,
                       start=omega, theta=theta2,  verb=0, nef=TRUE, wtol=wtol, tmax=20)
  }else{
    Wfit <- omega;
  }
  
  lambda <- omega%*%t(theta);
  counts2 <- as.matrix(X);
  temp <- counts2/lambda;
  t_matrix <- (t(temp) %*% omega)*theta;
  t_matrix[t_matrix == "NaN"] = 0
  
  lambda2 <- omega%*%t(theta2);
  temp2 <- counts2/lambda2;
  w_matrix <- (temp2 %*% theta)*omega;
  
  for(l in 1:length(pathway)){
    t_matrix[pathway[[l]], l] <- t_matrix[pathway[[l]], l] + alpha
  }
  
  omega <- ashnormalizetpx(w_matrix+(1/(n*K)), byrow=TRUE)
  full_indices <- which(omega==1, arr.ind=T)
  full_indices_rows <- unique(full_indices[,1]);
  omega[full_indices_rows,] <- omega[full_indices_rows,] + (1/(n*K));
  omega <- ashnormalizetpx(omega, byrow=TRUE)
  
  theta <- ashnormalizetpx(t_matrix, byrow=FALSE)
  
  move <- list(theta=theta, omega=omega)

  dif <- tol+1+qn
  update <- TRUE
  QNup <- tpxQN(move=move, Y=Y, X=X, alpha=alpha, verb=verb, admix=admix, grp=grp, doqn=qn-dif)
  move <- list(theta = move$theta, omega = QNup$move$omega)
  
  Y <- QNup$Y
  
  
  
  
