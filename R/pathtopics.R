##### Estimation for Topic Models ######

## intended main function; provides defaults and selects K via marginal lhd
pathtopics <- function(counts,
                   pathways,
                   shape=NULL, initopics=NULL, tol=0.1,
                   bf=FALSE, kill=2, verb=1, admix=TRUE,
                   nbundles=1, ord = TRUE, tmax=10000,...)
  ## tpxselect defaults: tmax=10000, wtol=10^(-4), qn=100, grp=NULL,
  ## nonzero=FALSE, dcut=-10, top_genes=100, burn_in=5
{
  if(is.null(pathways)){
    stop("User must provide a list of pathway variable names")
  }
  if(is.null(colnames(counts))){
    colnames(counts) <- 1:dim(counts)[2]
  }
  
  X <- CheckCounts(counts)
  p <- ncol(X)
  if(verb>0)
    cat(sprintf("\nEstimating on a %d document collection.\n", nrow(X)))

  ## check the prior parameters for theta
  if(prod(shape>0) != 1){ stop("use shape > 0\n") }


  pathways_indices <- list()
  l <- 1
  while(l <= length(pathways)){
     temp   <- match(pathways[[l]], colnames(counts))
     temp2 <- temp[!is.na(temp)]
     if(length(temp2) != 0){
       pathways_indices[[l]] <- temp2
     }
     l = l+1
  }
  
  if(length(pathways_indices) == 0){
    stop("The pathway names must match with column names of the counts matrix")
  }
  
  pathways_mat <- matrix(0, dim(counts)[2], length(pathways_indices))
  for(l in 1:length(pathways_indices)){
    pathways_mat[pathways_indices[[l]], l] <- rep(1, length(pathways_indices[[l]])) 
  }
  
  initopics <- pathnormalizetpx(pathways_mat, byrow = FALSE)
  K <- length(pathways_indices)
  
  tpx <- pathtpxfit(X=X, theta=initopics, pathways = pathways_indices,
                alpha=shape, tol=tol, verb=verb,
                admix=TRUE, method_admix=1, grp=NULL, tmax = tmax, 
                nbundles=1, light = 1,
                wtol=10^{-4}, qn=100,
                nonzero=FALSE, dcut=-10,
                top_genes=150, burn_in=5)

  K <- tpx$K

  ## clean up and out
  if(ord){ worder <- order(col_sums(tpx$omega), decreasing=TRUE) } # order by decreasing usage
  else{ worder <- 1:K }
  ## Main parameters
  theta=matrix(tpx$theta[,worder], ncol=K, dimnames=list(phrase=dimnames(X)[[2]], topic=paste(1:K)) )
  omega=matrix(tpx$omega[,worder], ncol=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  if(nrow(omega)==nrow(X)){ dimnames(omega)[[1]] <- dimnames(X)[[1]] }

  ## topic object
  out <- list(K=K, theta=theta, omega=omega, BF=tpx$BF, D=tpx$D, X=X)
  class(out) <- "topics"
  invisible(out) }
