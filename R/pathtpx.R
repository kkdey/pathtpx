
CheckCounts <- function(counts){
  if(class(counts)[1] == "TermDocumentMatrix"){ counts <- t(counts) }
  if(is.null(dimnames(counts)[[1]])){ dimnames(counts)[[1]] <- paste("doc",1:nrow(counts)) }
  if(is.null(dimnames(counts)[[2]])){ dimnames(counts)[[2]] <- paste("wrd",1:ncol(counts)) }
  empty <- row_sums(counts) == 0
  if(sum(empty) != 0){
    counts <- counts[!empty,]
    cat(paste("Removed", sum(empty), "blank documents.\n")) }
  return(as.simple_triplet_matrix(counts))
}


pathtpxfit <- function(X, theta, pathways, alpha, tol, verb,
                   admix, method_admix, grp, tmax, nbundles,
                   type, ind_model_indices, signatures, light,
                   wtol=10^{-4}, qn=100,
                   nonzero=FALSE, dcut=-10,
                   top_genes=150, burn_in=5)
{
  ## inputs and dimensions
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix") }
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
  
  system.time(omega <- pathtpxweights(n=n, p=p, xvo=xvo, wrd=wrd, doc=doc, start=pathtpxOmegaStart(X,theta2), theta=theta2))
  if(!admix){ omega <- matrix(apply(omega,2, function(w) tapply(w,grp,mean)), ncol=K) }
  
  iter <- 0
  dif <- tol+1+qn
  update <- TRUE
  if(verb>0){
    cat("log posterior increase: " )
    digits <- max(1, -floor(log(tol, base=10))) }
  
  Y <- NULL # only used for qn > 0
  Q0 <- col_sums(X)/sum(X)
  L <- pathtpxlpost(X=X, theta=theta2, omega=omega, alpha=alpha, admix=admix, grp=grp)
  if(is.infinite(L)){ L <- sum( (log(Q0)*col_sums(X))[Q0>0] ) }
  
  iter <- 1;
  
  ## Iterate towards MAP
  while( update  && iter < tmax ){
    
    if(light==2){
      if(admix && wtol > 0 && (iter-1)%%nbundles==0){
        if(iter <= burn_in){
          Wfit <- pathtpxweights(n=nrow(X), p=ncol(X), xvo=xvo, wrd=wrd, doc=doc,
                             start=omega, theta=theta,  verb=0, nef=TRUE, wtol=wtol, tmax=20)
        }
        if(iter > burn_in){
          # system.time(suppressWarnings(select_genes <- as.numeric(na.omit(as.vector(ExtractTopFeatures(theta, top_genes))))))
          ptm <- proc.time()
          
          topic_index <- apply(theta,1, which.max);
          select_genes <- vector();
          for(k in 1:K){
            topic_labels <- which(topic_index==k);
            select_genes <- c(select_genes,
                              topic_labels[order(apply(theta[which(topic_index==k),],1, sd),
                                                 decreasing=TRUE)[1:round(min(top_genes,(length(topic_labels)*0.9)))]]);
          }
          select_genes <- as.numeric(select_genes);
          
          counts <- as.matrix(X);
          counts.mod <- counts[,select_genes];
          Xmod <- CheckCounts(counts.mod)
          xvo2 <- Xmod$v[order(Xmod$i)]
          wrd2 <- Xmod$j[order(Xmod$i)]-1
          doc2 <- c(0,cumsum(as.double(table(factor(Xmod$i, levels=c(1:nrow(Xmod)))))))
          
          if(length(which(rowSums(counts.mod)==0))==0){
            Wfit <- pathtpxweights(n=nrow(Xmod), p=ncol(Xmod), xvo=xvo2, wrd=wrd2, doc=doc2,
                               start=omega, theta=theta[select_genes,],  verb=0, nef=TRUE, wtol=wtol, tmax=20)
          }else{
            Wfit <- matrix(0, nrow=nrow(X), ncol=ncol(X));
            Wfit[which(rowSums(counts.mod)==0),] <- omega[which(rowSums(counts.mod)==0),];
            Wfit[which(rowSums(counts.mod)!=0),] <- pathtpxweights(n=nrow(Xmod), p=ncol(Xmod), xvo=xvo2, wrd=wrd2, doc=doc2,
                                                               start=omega, theta=theta[select_genes,],  verb=0, nef=TRUE, wtol=wtol, tmax=20)
          }
          
          proc.time() - ptm
        }
      }else{
        Wfit <- omega;
      }}
    
    if(light==1){
      if(admix && wtol > 0 && (iter-1)%%nbundles==0){
        Wfit <- pathtpxweights(n=nrow(X), p=ncol(X), xvo=xvo, wrd=wrd, doc=doc,
                           start=omega, theta=theta2,  verb=0, nef=TRUE, wtol=wtol, tmax=20)
      }else{
        Wfit <- omega;
      }}
    
    if(light==0){
      Wfit <- omega;
    }
    
    move <- pathtpxEM(X=X, m=m, theta=theta, omega=Wfit, 
                  pathways = pathways, alpha=alpha, admix=admix,
                  method_admix=method_admix, grp=grp)
    
    QNup <- pathtpxQN(move=move, Y=Y, X=X, alpha=alpha, verb=verb, admix=admix, grp=grp, doqn=qn-dif)
    
    move <- list(theta = move$theta, omega = QNup$move$omega)
    
    Y <- QNup$Y
    
    if(QNup$L < L){  # happens on bad Wfit, so fully reverse
      if(verb > 10){ cat("_reversing a step_") }
      move <- pathtpxEM(X=X, m=m, theta=theta, omega=omega, pathways = pathways,
                    alpha=alpha, admix=admix,
                    method_admix=method_admix,grp=grp)
      QNup$L <-  pathtpxlpost(X=X, theta=move$theta, omega=move$omega, alpha=alpha, admix=admix, grp=grp) 
      
      move <- list(theta = move$theta, omega = QNup$move$omega)

    }
    
    dif <- (QNup$L-L)
    
    L <- QNup$L
    
    
    ## check convergence
    if(abs(dif) < tol){
      if(sum(abs(theta-move$theta)) < tol){ update = FALSE } }
    
    ## print
    if(verb>0 && (iter-1)%%ceiling(10/verb)==0 && iter>0){
      ##if(verb>0 && iter>0){
      cat( paste( round(abs(dif),digits), #" (", sum(abs(theta-move$theta)),")",
                  ", ", sep="") ) }
    
    ## heartbeat for long jobs
    if(((iter+1)%%1000)==0){
      cat(sprintf("p %d iter %d diff %g\n",
                  nrow(theta), iter+1,round(dif))) }
    
    ## iterate
    iter <- iter+1
    theta <- move$theta
    omega <- move$omega
    
  }
  
  L <- pathtpxlpost(X=X, theta=theta, omega=omega, alpha=alpha, admix=admix, grp=grp)
  
  if(verb>0){
    cat("done.")
    if(verb>1) { cat(paste(" (L = ", round(L,digits), ")", sep="")) }
    cat("\n")
  }
  
  out <- list(theta=theta, omega=omega, K=K, alpha=alpha, L=L, iter=iter)
  invisible(out) 
  
}

pathtpxweights <- function(n, p, xvo, wrd, doc, start, theta, verb=FALSE, nef=TRUE, wtol=10^{-5}, tmax=1000)
{
  K <- ncol(theta)
  start[start == 0] <- 0.1/K
  start <- start/rowSums(start)
  omega <- .C("Romega",
              n = as.integer(n),
              p = as.integer(p),
              K = as.integer(K),
              doc = as.integer(doc),
              wrd = as.integer(wrd),
              X = as.double(xvo),
              theta = as.double(theta),
              W = as.double(t(start)),
              nef = as.integer(nef),
              tol = as.double(wtol),
              tmax = as.integer(tmax),
              verb = as.integer(verb),
              PACKAGE="maptpx")
  return(t(matrix(omega$W, nrow=ncol(theta), ncol=n))) }

pathtpxEM <- function(X, m, theta, omega, pathways, alpha, admix, method_admix, grp)
{
  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(theta)
  
  if(admix==TRUE & method_admix==1){
    
    theta2 <- theta
    theta2[theta2 == 1] <- 1 - 1e-14
    theta2[theta2 == 0] <- 1e-14
    
    lambda <- omega%*%t(theta);
    counts2 <- as.matrix(X);
    temp <- counts2/lambda;
    t_matrix <- (t(temp) %*% omega)*theta;
    t_matrix[t_matrix == "NaN"] = 0
    
    lambda2 <- omega%*%t(theta2);
    temp2 <- counts2/lambda2;
    w_matrix <- (temp2 %*% theta)*omega;
    
    for(l in 1:length(pathways)){
      t_matrix[pathways[[l]], l] <- t_matrix[pathways[[l]], l] + alpha
    }
    
    omega <- pathnormalizetpx(w_matrix+(1/(n*K)), byrow=TRUE)
    full_indices <- which(omega==1, arr.ind=T)
    full_indices_rows <- unique(full_indices[,1]);
    omega[full_indices_rows,] <- omega[full_indices_rows,] + (1/(n*K));
    omega <- pathnormalizetpx(omega, byrow=TRUE)
    theta <- pathnormalizetpx(t_matrix, byrow=FALSE)
  }
  
  if(admix==TRUE & method_admix==2){ Xhat <- (X$v/pathtpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j))*(omega[X$i,]*theta[X$j,])
  Zhat <- .C("Rzhat", n=as.integer(n), p=as.integer(p), K=as.integer(K), N=as.integer(nrow(Xhat)),
             Xhat=as.double(Xhat), doc=as.integer(X$i-1), wrd=as.integer(X$j-1),
             zj = as.double(rep(0,K*p)), zi = as.double(rep(0,K*n)), PACKAGE="maptpx")
  theta <- pathnormalizetpx(matrix(Zhat$zj+alpha, ncol=K), byrow=FALSE)
  omega <- pathnormalizetpx(matrix(Zhat$zi+1/K, ncol=K)) }
  if(!admix){
    qhat <- pathtpxMixQ(X, omega, theta, grp, qhat=TRUE)$qhat
    ## EM update
    theta <- pathnormalizetpx(tcrossprod_simple_triplet_matrix( t(X), t(qhat) ) + alpha, byrow=FALSE)
    omega <- pathnormalizetpx(matrix(apply(qhat*m,2, function(x) tapply(x,grp,sum)), ncol=K)+1/K )  }
  
  return(list(theta=theta, omega=omega)) }

pathtpxQN <- function(move, Y, X, alpha, verb, admix, grp, doqn)
{
  move$theta[move$theta==1] <- 1 - 1e-14;
  move$omega[move$omega==1] <- 1 - 1e-14;
  move$omega[move$omega==0] <- 1e-14;
  move$theta[move$theta==0] <- 1e-14;
  move$theta <- pathnormalizetpx(move$theta, byrow = FALSE)
  move$omega <- pathnormalizetpx(move$omega, byrow = TRUE)
  
  ## always check likelihood
  L <- pathtpxlpost(X=X, theta=move$theta, omega=move$omega,
                alpha=alpha, admix=admix, grp=grp)
  
  if(doqn < 0){ return(list(move=move, L=L, Y=Y)) }
  
  ## update Y accounting
  Y <- cbind(Y, pathtpxToNEF(theta=move$theta, omega=move$omega))
  if(ncol(Y) < 3){ return(list(Y=Y, move=move, L=L)) }
  if(ncol(Y) > 3){ warning("mis-specification in quasi-newton update; please report this bug.") }
  
  ## Check quasinewton secant conditions and solve F(x) - x = 0.
  U <- as.matrix(Y[,2]-Y[,1])
  V <- as.matrix(Y[,3]-Y[,2])
  sUU <- sum(U^2)
  sVU <- sum(V*U)
  Ynew <- Y[,3] + V*(sVU/(sUU-sVU))
  qnup <- pathtpxFromNEF(Ynew, n=nrow(move$omega),
                     p=nrow(move$theta), K=ncol(move$theta))
  
  ## check for a likelihood improvement
  Lqnup <- try(pathtpxlpost(X=X, theta=qnup$theta, omega=qnup$omega,
                        alpha=alpha, admix=admix, grp=grp), silent=TRUE)
  
  if(inherits(Lqnup, "try-error")){
    if(verb>10){ cat("(QN: try error) ") }
    return(list(Y=Y[,-1], move=move, L=L)) }
  
  if(verb>10){ cat(paste("(QN diff ", round(Lqnup-L,3), ")\n", sep="")) }
  
  if(Lqnup < L){
    return(list(Y=Y[,-1], move=move, L=L)) }
  else{
    L <- Lqnup
    Y <- cbind(Y[,2],Ynew)
    return( list(Y=Y, move=qnup, L=L) )
  }
}

pathtpxlpost <- function(X, theta, omega, alpha, admix=TRUE, grp=NULL)
{
  theta[theta==1] <- 1 - 1e-10;
  omega[omega==1] <- 1 - 1e-10;
  omega[omega==0] <- 1e-10;
  theta[theta==0] <- 1e-10;
  theta <- pathnormalizetpx(theta, byrow = FALSE)
  omega <- pathnormalizetpx(omega, byrow = TRUE)
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }
  K <- ncol(theta)
  
  if(admix){ L <- sum( X$v*log(pathtpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)) ) }
  else{ L <- sum(pathtpxMixQ(X, omega, theta, grp)$lqlhd) }
  if(is.null(nrow(alpha))){ if(alpha != 0){ L <- L + sum(alpha*log(theta))  } } # unnormalized prior
  L <- L + sum(log(omega))/K
  
  return(L) }

pathtpxQ <- function(theta, omega, doc, wrd){
  
  theta[theta==1] <- 1 - 1e-14;
  theta[theta==0] <- 1e-14;
  theta <- pathnormalizetpx(theta, byrow = FALSE)
  
  omega[omega==1] <- 1 - 1e-14;
  omega[omega==0] <- 1e-14;
  omega <- pathnormalizetpx(omega, byrow = TRUE)
  
  if(length(wrd)!=length(doc)){stop("index mis-match in tpxQ") }
  if(ncol(omega)!=ncol(theta)){stop("theta/omega mis-match in tpxQ") }
  
  out <- .C("RcalcQ",
            n = as.integer(nrow(omega)),
            p = as.integer(nrow(theta)),
            K = as.integer(ncol(theta)),
            doc = as.integer(doc-1),
            wrd = as.integer(wrd-1),
            N = as.integer(length(wrd)),
            omega = as.double(omega),
            theta = as.double(theta),
            q = double(length(wrd)),
            PACKAGE="maptpx" )
  
  return( out$q ) }


pathtpxMixQ <- function(X, omega, theta, grp=NULL, qhat=FALSE){
  if(is.null(grp)){ grp <- rep(1, nrow(X)) }
  
  theta[theta==1] <- 1 - 1e-14;
  theta[theta==0] <- 1e-14;
  theta <- pathnormalizetpx(theta, byrow = FALSE)
  
  omega[omega==1] <- 1 - 1e-14;
  omega[omega==0] <- 1e-14;
  omega <- pathnormalizetpx(omega, byrow = TRUE)
  
  K <- ncol(omega)
  n <- nrow(X)
  mixhat <- .C("RmixQ",
               n = as.integer(nrow(X)),
               p = as.integer(ncol(X)),
               K = as.integer(K),
               N = as.integer(length(X$v)),
               B = as.integer(nrow(omega)),
               cnt = as.double(X$v),
               doc = as.integer(X$i-1),
               wrd = as.integer(X$j-1),
               grp = as.integer(as.numeric(grp)-1),
               omega = as.double(omega),
               theta = as.double(theta),
               Q = double(K*n),
               PACKAGE="maptpx")
  ## model and component likelihoods
  lQ <- matrix(mixhat$Q, ncol=K)
  lqlhd <- log(row_sums(exp(lQ)))
  lqlhd[is.infinite(lqlhd)] <- -600 # remove infs
  if(qhat){
    qhat <- exp(lQ-lqlhd)
    ## deal with numerical overload
    infq <- row_sums(qhat) < .999
    if(sum(infq)>0){
      qhat[infq,] <- 0
      qhat[n*(apply(matrix(lQ[infq,],ncol=K),1,which.max)-1) + (1:n)[infq]] <- 1 }
  }
  return(list(lQ=lQ, lqlhd=lqlhd, qhat=qhat)) }

pathtpxToNEF <- function(theta, omega){
  n <- nrow(omega)
  p <- nrow(theta)
  K <- ncol(omega)
  return(.C("RtoNEF",
            n=as.integer(n), p=as.integer(p), K=as.integer(K),
            Y=double((p-1)*K + n*(K-1)),
            theta=as.double(theta), tomega=as.double(t(omega)),
            PACKAGE="maptpx")$Y)
}

## 'From' NEF representation back to probabilities
pathtpxFromNEF <- function(Y, n, p, K){
  bck <- .C("RfromNEF",
            n=as.integer(n), p=as.integer(p), K=as.integer(K),
            Y=as.double(Y), theta=double(K*p), tomega=double(K*n),
            PACKAGE="maptpx")
  return(list(omega=t( matrix(bck$tomega, nrow=K) ), theta=matrix(bck$theta, ncol=K)))
}


pathtpxOmegaStart <- function(X, theta)
{
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }
  omega <- try(tcrossprod_simple_triplet_matrix(X, solve(t(theta)%*%theta)%*%t(theta)), silent=TRUE )
  if(inherits(omega,"try-error")){ return( matrix( 1/ncol(theta), nrow=nrow(X), ncol=ncol(theta) ) ) }
  omega[omega <= 0] <- .5
  return( pathnormalizetpx(omega, byrow=TRUE) )
}





