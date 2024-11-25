rm(list=ls()); graphics.off(); cat("\014")

source("core.R")

d <- 2
N <- 10

K.cjko <- 5000
patience.cjko <- 10

all.eps <- c(0.05, 0.1, 0.25) # 0.5, 0.25, 0.1, 0.05, 0.01
all.res <- list()


Ff <- F.sqdW2
set.seed(42)

# example for Panos-70
Y <- rmvnorm( n=N, mean=c(5,-4), sigma=2.5*diag(x=1,nrow=d) )
X0 <- rmvnorm( n=N, mean=c(-7,3), sigma=1*diag(1,nrow=d) )



for( eps in all.eps ) {
  
  cat("> Starting constrained-JKO with eps =",eps,"\n")
  cjko.res <- cJKO( X=X0, Ff=Ff, eps=eps, K=K.cjko, patience=patience.cjko )
  cat("\n> Solved in:",sum(cjko.res$runTimes),"[secs]\n\n")
  all.res[[length(all.res)+1]] <- cjko.res
  
}

clrs <- c("deepskyblue","green4","red")
clrs2 <- c("lightblue","lightgreen","pink")

par(mfrow=c(1,4))

par(mar=c(4.1,4.6,1.1,1.1))
for( j in 1:length(all.res) ) {
  res <- all.res[[j]]
  K <- length(res$Xks)
  plot( X0[,1], X0[,2], pch=19, col=clrs[j],
        xlim=c(-10,10), ylim=c(-10,10),
        xlab=expression(x[1]), ylab=expression(x[2]),
        cex.lab=2, cex.axis=2 )
  for( i in 1:nrow(X0) ) {
    lines( c(X0[i,1],res$Xks[[1]][i,1]), c(X0[i,2],res$Xks[[1]][i,2]), col=clrs2[j] )
    for( k in 2:K ) {
      lines( c(res$Xks[[k-1]][i,1],res$Xks[[k]][i,1]),
             c(res$Xks[[k-1]][i,2],res$Xks[[k]][i,2]),
             col=clrs2[j], lwd=2 )
    }
  }
  points( X0[,1], X0[,2], pch=4, cex=1.5, col=clrs[j], lwd=2 )
  points( res$Xks[[K]][,1], res$Xks[[K]][,2], pch=19, col=clrs[j] ) 
  points( Y[,1], Y[,2], pch=1, lwd=2, cex=1.7 ) 
  
  if( j==1 ) {
    legend( "topright", legend=c("source","transported","target"),
            col=c("deepskyblue","deepskyblue","black"), pch=c(4,19,1), cex=1.7,
            lwd=2, lty=0, pt.cex=c(1.7,1.7,2)  )
  }
}


par(mar=c(5.1,5.6,2.1,1.1))
for( j in 1:length(all.res) ) {
  res <- all.res[[j]]
  if( j==1 ) {
    plot( res$objs, type="l", lwd=3, col=clrs[j],
          ylab=expression(F(P)),
          xlab="cJKO iterations",
          cex.lab=2, cex.axis=2 )
  } else {
    lines( res$objs, lwd=3, col=clrs[j] )
  }
}
legend( "topright", cex=2,
        title=expression(epsilon),
        legend=all.eps, col=clrs, lwd=3 )
