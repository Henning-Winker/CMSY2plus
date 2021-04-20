#' mvn()
#'
#' Function to create multivariate-normal distribution for r-k for plotting
#' @param n number of r-k Monte-Carlo iteratition  
#' @param mean.log.r expected mean
#' @param sd.log.r sd on log-scale
#' @param mean.log.k expected mean
#' @param sd.log.k sd on log-scale
#' @return posterior of mvn.log.rk
#' @export 
mvn   <- function(n,mean.log.r,sd.log.r,mean.log.k,sd.log.k) {
  cov.log.rk <- cor.log.rk*sd.log.r*sd.log.k # covariance with empirical correlation and prior variances  covar.log.rk = matrix(NA, ncol=2,nrow=2)   # contract covariance matrix
  covar.log.rk      <- matrix(NA, ncol=2,nrow=2) # covariance matrix
  covar.log.rk[1,1] <- sd.log.r^2                # position [1,1] is variance of log.r
  covar.log.rk[2,2] <- sd.log.k^2               # position [2,2] is variance of log.k
  covar.log.rk[1,2] = covar.log.rk[2,1] = cov.log.rk     # positions [1,2] and [2,1] are correlations
  mu.log.rk  <- (c(mean.log.r,mean.log.k))      # vector of log.means
  mvn.log.rk <- rmvnorm(n,mean=mu.log.rk,sigma=covar.log.rk,method="svd")
  return(mvn.log.rk)
}



#' beta prior for relative biomass b/k
#'
#' convert b.prior ranges into beta priors
#' @param b.prior lower and upper range
#' @return vector of beta a and b parameters (shape and scale)
#' @export
beta.prior = function(b.prior){
  bk.beta = matrix(0,nrow = 3,ncol=3)
  for(i in 1:3){
    sd.bk = (b.prior[2,i]-b.prior[1,i])/(4*0.98)
    mu.bk = mean(b.prior[1:2,i])
    cv.bk = sd.bk/mu.bk
    bk.beta[1:2,i] = get_beta(mu.bk,cv.bk)
  }
  return(bk.beta)
}

#' traceEllipse()
#'
#' Fits an ellipse around the CMSY r-k cloud and estimates the rightmost focus
#' @param rs posterior of r 
#' @param ks posterior of k 
#' @return COM bias correct ellipse space for r-k
#' @export

traceEllipse<-function(rs,ks){
  log.rs<-log(rs)
  log.ks<-log(ks)
  

  #prepare data for ellipse fitting
  cloud.data <- as.matrix(data.frame(x = log.rs, y = log.ks))
  ellip <- EllipseDirectFit(cloud.data)
  #estimate ellipse characteristics
  atog<-AtoG(ellip)
  ellipG <- atog$ParG
  ell.center.x<-ellipG[1]
  ell.center.y<-ellipG[2]
  ell.axis.a<-ellipG[3]
  ell.axis.b<-ellipG[4]
  ell.tilt.angle.deg<-180/pi*ellipG[5]
  ell.slope<-tan(ellipG[5])
  xy.ell<-calculateEllipse(ell.center.x,
                           ell.center.y,
                           ell.axis.a,
                           ell.axis.b,
                           ell.tilt.angle.deg)
  ell.intercept.1 = ell.center.y-ell.center.x*ell.slope
  ell.demiaxis.c.sqr<-(0.25*ell.axis.a*ell.axis.a)-(0.25*ell.axis.b*ell.axis.b)
  if (ell.demiaxis.c.sqr<0)
    ell.demiaxis.c.sqr<-ell.axis.a/2
  else
    ell.demiaxis.c<-sqrt(ell.demiaxis.c.sqr)
  sin.c<-ell.demiaxis.c*sin(ellipG[5])
  cos.c<-ell.demiaxis.c*cos(ellipG[5])
  ell.foc.y<-ell.center.y-sin.c
  ell.foc.x<-ell.center.x-cos.c
  
  return (c(exp(ell.foc.x),exp(ell.foc.y)))
}