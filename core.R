library(transport)
library(nloptr)
library(EMCluster)
library(MGMM)
library(LaplacesDemon)
library(mvtnorm)


F.sqdW2 <- function( X ) {
  ot <- transport( pp(X), pp(Y), p=2, method="primaldual" )
  F.W2 <- (wasserstein( pp(X), pp(Y), p=2, tplan=ot))^2
}

F.symKL <- function(X) {
  
  # grid of locations for computing KL divergence
  XX <- as.matrix( expand.grid(seq(-5,15,length.out=ng),seq(-10,10,length.out=ng)) )

  # clustering the data of the point cloud X via EM
  em <- emcluster(X,emobj=init.EM(X))
  
  means <- list()
  for( i in 1:nrow(em$Mu) )
    means[[i]] <- em$Mu[i,]
  
  gmm <- FitGMM(data=X, k=em$nclass, init_means=means )
  
  # probability density in XX for the point cloud X (as the Gaussian Mixture obtained via EM)  
  if( class(gmm)=="mix" ) {
    qx <- 0
    for( i in 1:em$nclass )
      qx <- qx + gmm@Proportions[i] * mvtnorm::dmvnorm(XX,mean=gmm@Means[[i]],sigma=gmm@Covariances[[i]])
  } else {
    qx <- mvtnorm::dmvnorm(XX,mean=gmm@Mean,sigma=gmm@Covariance)
  }
  
  # probability density in XX for the target distribution
  if( is.list(target.mean) ) {
    px <- 0
    for( i in 1:length(target.mean) )
      px <- px + (1/length(target.mean)) * mvtnorm::dmvnorm( XX, mean=target.mean[[i]], sigma=target.sigma[[i]] )  
  } else {
    px <- mvtnorm::dmvnorm( XX, mean=target.mean, sigma=target.sigma )  
  }
  
  tmp <- KLD(px,qx)
  F.symKL <- 2 * tmp$mean.sum.KLD
}

# Extended to deal with d>2 [RM computation has been changed; alpha is now a matrix with size N * (d-1) )
push_fwd <- function( X, lambda, alpha ) {
  X_ <- X
  d <- ncol(X)
  N <- nrow(X)
  stopifnot( is.matrix(X) && length(alpha)==N*(d-1) )
  alpha <- matrix( alpha, N, d-1, byrow=T )
  for( i in 1:N ) {
    for( j in 1:(d-1) ) {
      R0 <- diag(d)
      R0[ j:(j+1), j:(j+1) ] <- matrix( c(cos(alpha[i,j]), -sin(alpha[i,j]), sin(alpha[i,j]), cos(alpha[i,j])), nrow=2, byrow=T )
      if( j==1 ) {
        RM <- R0
      } else {
        RM <- RM %*% R0
      }
    }
    X_[i,] <- X[i,] - (lambda[i] * t(rep(sqrt(d),d))) %*% RM
  }
  push_fwd <- X_
}

# old version: only for d=2
# push_fwd <- function( X, lambda, alpha ) {
#   stopifnot( is.matrix(X) && ncol(X)==2 )
#   X_ <- X
#   for( i in 1:nrow(X) ) {
#     RM <- matrix( c(cos(alpha[i]), -sin(alpha[i]), sin(alpha[i]), cos(alpha[i])), nrow=ncol(X), byrow=T )
#     X_[i,] <- X[i,] - (lambda[i] * rep(sqrt(ncol(X)),ncol(X))) %*% RM
#   }
#   push_fwd <- X_
# }

jko.obj <- function( x, Ff, Xk, h ) {
  lambda <- x[1:nrow(Xk)]
  alpha <- x[-c(1:nrow(Xk))]
  X_ <- push_fwd(X=Xk,lambda=lambda,alpha=alpha)
  W2 <- sqrt(sum((Xk - X_)^2)/nrow(Xk))
  term1 <- (W2^2)/(2*h)
  term2 <- Ff(X_)
  jko.obj <- term1 + term2
}

JKO <- function( X0, Ff, h, K, patience=10 ) {
  d <- ncol(X0)
  m <- nrow(X0)
  Xks <- list()
  X <- X0 
  objs <- NULL
  Lambdas <- Alphas <- list()
  y.best <- Inf
  notImproved <- 0

  k <- 1
  runTimes <- NULL
  while( k <=K && notImproved<patience ) {
    if( k %% 10 == 0 )
      cat(".")
    if( k %% 500 == 0 )
      cat("\n")
    
    elapsed <- Sys.time()
    
    lambda0 <- runif( n=m(X), min=0, max=sqrt(d) )
    alpha0 <- runif( n=(d-1)*m, min=0, max=2*pi )

    res <- optim( par=c(lambda0,alpha0), fn=jko.obj, gr=NULL,
                  method="L-BFGS-B",
                  lower=numeric(m+(d-1)*m),
                  upper=c(rep(Inf,m),rep(2*pi,(d-1)*m)),
                  control=list(trace=0),
                  Ff=Ff, Xk=X, h=h )

    elapsed <- difftime(Sys.time(),elapsed,units="secs")
    runTimes <- c(runTimes,elapsed)
    
    lambda <- res$par[1:m]
    alpha <- res$par[-c(1:m)]
    objs <- c(objs,res$value)
    
    if( y.best - res$value >= 10^-6 ) {
      y.best <- res$value
      notImproved <- 0
    } else {
      notImproved <- notImproved+1
    }

    
    Lambdas[[k]] <- lambda
    Alphas[[k]] <- alpha
  
    X <- push_fwd(X=X,lambda=lambda,alpha=alpha)
    Xks[[k]] <- X
    
    k <- k+1
  }
  
  return( list(Xks=Xks, Lambdas=Lambdas, Alphas=Alphas, objs=objs, runTimes=runTimes) )
}


cjko.obj <- function( x, Ff, Xk, eps ) {
  lambda <- x[1:nrow(Xk)]
  alpha <- x[-c(1:nrow(Xk))]
  X_ <- push_fwd( X=Xk, lambda=lambda, alpha=alpha )
  cjko.obj <- Ff( X_ )
}

cjko.constr <- function( x, Ff, Xk, eps ) {
  lambda <- x[1:nrow(Xk)]
  cjko.constr <- sum(sqrt(lambda))/nrow(Xk) - eps
}

cJKO <- function( X0, Ff, eps, K, patience=10 ) {
  d <- ncol(X0)
  m <- nrow(X0)
  Xks <- list()
  X <- X0 
  objs <- W2s <- NULL
  Lambdas <- Alphas <- list()
  y.best <- Inf
  notImproved <- 0
  
  k <- 1
  runTimes <- NULL
  while( k <=K && notImproved<patience ) {
    
    if( k %% 10 == 0 )
      cat(".")
    if( k %% 500 == 0 )
      cat("\n")
    
    elapsed <- Sys.time()

    res <- NULL
    while( is.null(res) ) {
      lambda0 <<- runif( n=m, min=0, max=(m*eps)^2 )
      if( sum(sqrt(lambda0))>m*eps )
        lambda0 <<- ((m*eps) * (sqrt(lambda0) / sum(sqrt(lambda0))))^2
      alpha0 <<- runif( n=(d-1)*m, min=0, max=2*pi )
      
      try( expr = (res <- nloptr( x0=c(lambda0,alpha0), eval_f=cjko.obj, eval_g_ineq=cjko.constr,
                                  lb=numeric(m+(d-1)*m),
                                  ub=c(rep((m*eps)^2,m),rep(2*pi,(d-1)*m)),
                                  # opts=list(algorithm="NLOPT_GN_ISRES",xtol_rel=10^-6), [-]
                                  opts=list(algorithm="NLOPT_LN_COBYLA",xtol_rel=10^-6), # [+]
                                  Ff=Ff, Xk=X, eps=eps ) ),
           silent = T )
      if( is.null(res) )
        cat("(!)")
    }

    elapsed <- difftime(Sys.time(),elapsed,units="secs")
    runTimes <- c(runTimes,elapsed)
    
    lambda <- res$solution[1:nrow(X)]
    alpha <- res$solution[-c(1:nrow(X))]
    objs <- c(objs,res$objective)
    W2s <- c(W2s,sum(sqrt(lambda)))
    
    if( y.best-res$objective >= 10^-6 ) {
      y.best <- res$objective
      notImproved <- 0
    } else {
      notImproved <- notImproved+1
    }
    
    Lambdas[[k]] <- lambda
    Alphas[[k]] <- alpha
    
    X <- push_fwd(X=X,lambda=lambda,alpha=alpha)
    
    Xks[[k]] <- X
    
    k <- k+1
  }

  
  return( list(Xks=Xks, Lambdas=Lambdas, Alphas=Alphas, objs=objs, W2s=W2s, runTimes=runTimes) )
}