rm(list=ls()); graphics.off(); cat("\014")

library(mvtnorm)

N <- 30 # must be a multiple of 3
ng <- 300
XX <- as.matrix(expand.grid( x1=seq(-5,15,length.out=ng), x2=seq(-10,10,length.out=ng) ))

par(mfrow=c(1,3))
par(mar=c(3.1,3.1,2.1,1.1))

# G2G

set.seed(42)
target.mean <- c(10,5)
target.sigma <- 2*matrix( c(1,-0.75,-0.75,1), nrow=ncol(XX) )
Y <- rmvnorm( n=N, mean=target.mean, sigma=target.sigma )
source.mean <- c(0,0)
source.sigma <- diag(x=1,nrow=ncol(XX))
X0 <- rmvnorm( n=N, mean=source.mean, sigma=source.sigma )

d.src <- dmvnorm(XX,mean=source.mean,sigma=source.sigma)
d.tgt <- dmvnorm(XX,mean=target.mean,sigma=target.sigma)

contour( x=unique(XX[,1]), y=unique(XX[,2]), z=matrix(d.src,ng,ng),
         nlevels=7, col="blue", lwd=2, main="MVN to MVN",
         cex.axis=2, cex.lab=2, cex.main=2, drawlabels=F  )
contour( x=unique(XX[,1]), y=unique(XX[,2]), z=matrix(d.tgt,ng,ng),
         nlevels=7, col="firebrick", lwd=2, drawlabels=F, add=T )

# points( X0[,1], X0[,2], pch=19, col="blue")
# points( Y[,1], Y[,2], pch=19, col="red")

legend( "bottomleft", legend=c("source","target"),
        col=c("blue","firebrick"), lwd=3, cex=2 )

# G2Mx

d1 <- dmvnorm(XX,mean=c(10,7),sigma=1.5*diag(2))
d2 <- dmvnorm(XX,mean=c(5,5),sigma=2*matrix(c(2,0.5,0.5,1),2,2))
d3 <- dmvnorm(XX,mean=c(10,-2.5),sigma=matrix(c(2.5,-1.5,-1.5,2.5),2,2))
d.tgt <- (d1+d2+d3)/3

Y <- rmvnorm(N/3,mean=c(10,7),sigma=1.5*diag(2))
Y <- rbind(Y,rmvnorm(N/3,mean=c(5,5),sigma=2*matrix(c(2,0.5,0.5,1),2,2)))
Y <- rbind(Y,rmvnorm(N/3,mean=c(10,-2.5),sigma=matrix(c(2.5,-1.5,-1.5,2.5),2,2)))


contour( x=unique(XX[,1]), y=unique(XX[,2]), z=matrix(d.src,ng,ng),
         nlevels=7, col="blue", lwd=2, drawlabels= F,
         main="MVN to Gaussian Mixture", cex.lab=2, cex.axis=2, cex.main=2 )
contour( x=unique(XX[,1]), y=unique(XX[,2]), z=matrix(d.tgt,ng,ng),
         nlevels=7, col="firebrick", lwd=2, drawlabels=F, add=T )

# points( X0[,1], X0[,2], pch=19, col="blue")
# points( Y[,1], Y[,2], pch=19, col="red")



# Mx2G

set.seed(42)
target.mean <- c(5,-5)
target.sigma <- diag(2)
Y <- rmvnorm( n=N, mean=target.mean, sigma=target.sigma )
d.tgt <- dmvnorm(XX,mean=target.mean,sigma=target.sigma)

source.mean <- c(0,0)
source.sigma <- diag(x=1,nrow=ncol(XX))

d1 <- dmvnorm(XX,mean=c(0,5),sigma=matrix(c(2,0.5,0.5,2),2,2))
d2 <- dmvnorm(XX,mean=c(5,5),sigma=2*matrix(c(2,-0.5,-0.5,2),2,2))
d3 <- dmvnorm(XX,mean=c(12,0),sigma=matrix(c(2.5,1.5,1.5,2.5),2,2))
d.src <- (d1+d2+d3)/3

X0 <- rmvnorm(N/3,mean=c(0,5),sigma=matrix(c(2,0.5,0.5,2),2,2))
X0 <- rbind(X0,rmvnorm(N/3,mean=c(5,5),sigma=2*matrix(c(2,-0.5,-0.5,1),2,2)))
X0 <- rbind(X0,rmvnorm(N/3,mean=c(12,0),sigma=matrix(c(2,1,1,2),2,2)))


contour( x=unique(XX[,1]), y=unique(XX[,2]), z=matrix(d.src,ng,ng),
         nlevels=7, col="blue", main="Gaussian Mixture to MVN",
         cex.axis=2, cex.main=2, cex.lab=2, lwd=2, drawlabels=F )
contour( x=unique(XX[,1]), y=unique(XX[,2]), z=matrix(d.tgt,ng,ng),
         nlevels=7, col="firebrick", lwd=2, drawlabels=F, add=T )

# points( X0[,1], X0[,2], pch=19, col="blue")
# points( Y[,1], Y[,2], pch=19, col="red")
