### Functions for data process
#' Generate initialization from specific space
#'
#' @param p p*p orthodox matrix
#' @param d p*d orthodox matrix
#' @param runs the number of orthodox matrix
#' @param orth.method orthodox method
#'
#' @return w.list
#' @export
#'
#' @examples gen.inits(2,3,3,'svd')
gen.inits <- function(p,d,runs,orth.method=c('svd','givens')) {
  orth.method=match.arg(orth.method) # the first value in orth.metod #
  W.list = list()
  for(i in 1:runs) {
    if(orth.method=='givens') {
      W.list[[i]] <- as.matrix(theta2W(runif(n=choose(p,2),min=0,max=2*pi)))[,1:d]  # convert vector into orthodox matirx and get p*d #
    } else {
      temp = matrix(rnorm(p*d),p,d)
      W.list[[i]] <- svd(temp)$u     # svd left matrix p*d #
    }
  }
  W.list
}

#' Convert angle vector into  orthodox matrix
#'
#' @param theta vector of angles theta
#'
#' @return an orthodox matrix
#' @export
#'
#'
theta2W = function(theta)
{
  #<author hidden>
  # For a vector of angles theta, returns W, a d x d Givens rotation matrix:
  # W = Q.1,d %*% ... %*% Q.d-1,d %*% Q.1,d-1 %*% ... %*% Q.1,3 %*% Q.2,3 %*% Q.1,2
  ##  if(theta < 0  || pi < theta){stop("theta must be in the interval [0,pi]")}
  d = (sqrt(8*length(theta)+1)+1)/2
  if(d - floor(d) != 0){stop("theta must have length: d(d-1)/2")}
  W = diag(d)
  index = 1
  for(j in 1:(d-1)){
    for(i in (j+1):d){
      Q.ij = givens.rotation(theta[index], d, c(i,j))
      W = Q.ij %*% W
      index = index + 1
    }
  }
  W
}


#-----------------------
#' Orthogonalization of matrix
#'
#' @param W arbitrary matrix
#'
#' @return orthogonalized matrix
#' @export
#'
orthogonalize = function (W) {
  ##For arbitrary W, this is equivalent to (WW^T)^{-1/2} W
  temp <- svd(W)
  tcrossprod(temp$u,temp$v)
}




#' Whitening Function
#'
#' @param X dataset
#' @param n.comp the number of components
#' @param center.row whether center the row of data
#'
#' @return a whitener matrix
#' @export
#'
#' @import MASS
whitener <- function(X,n.comp=ncol(X),center.row=FALSE) {

  #X must be n x d
  if(ncol(X)>nrow(X)) warning('X is whitened with respect to columns')
  #Corresponds to model of the form X.center = S A, where S are orthogonal with covariance = identity.
  x.center=scale(X,center=TRUE,scale=FALSE)
  if(center.row==TRUE) x.center = x.center - rowMeans(x.center)
  n.rep=dim(x.center)[1]
  svd.x=svd(x.center,nu=n.comp,nv=n.comp)

  #RETURNS PARAMETERIZATION AS IN fastICA (i.e., X is n x d)
  #NOTE: For whitened X, re-whitening leads to different X
  #The output for square A is equivalent to solve(K)
  return(list(whitener=t(ginv(svd.x$v%*%diag(svd.x$d[1:n.comp])/sqrt(n.rep-1))),Z=sqrt(n.rep-1)*svd.x$u,mean=apply(X,2,mean)))
}

###############################################################
# Auxillary functions taken from other people
##############################################################

# Original function by Yunfeng Zhang - calculate the power of a matrix
#' Calculate the power of a square matrix
#'
#' returns a matrix composed of eigenvector x diag(eigenvalue ^ power) x eigenvector'
#' @param S a square matrix
#' @param power the times of power
#'
#' @return a matrix after power calculation that eigenvector x diag(eigenvalue ^ power) x eigenvector'
#' @export
#'
#' @examples a <- matrix(1:9,3,3)
#' a %^% 2
#'
"%^%" <- function(S, power){
  out = eigen(S)
  nonzero = which(out$values > 1e-8)
  out$vectors[,nonzero] %*% (out$values[nonzero]^power * t(out$vectors[,nonzero]))  # which equals to out$vectors[,nonzero] %*% (diag(out$values[nonzero]^power) %*% t(out$vectors[,nonzero]))
}

#######################################################
# Objective functions - to monitor the convergence
####################################################
# Function that calculates JB
#' JB score calculate
#'
#' @param U U matrix for matched columns rj x n
#' @param X whitened data matrix n x px, data = whitenerXA \%*\% dXcentered
#' @param S the variable loadings r x px.
#' @param alpha default = 0.8
#'
#' @return the sum of JB score for each component, the row of Sj.
#' @export
#'
calculateJB <- function(U=NULL, X=NULL, S=NULL, alpha = 0.8){
  if((is.null(U)|is.null(X))&is.null(S)){stop("At least input S matrix or both U and X matrices.")}
  if(is.null(S)) {# Calculate u^{top}X for each u_l and Xj, this is UX
  S = U %*% X # Sx matrix with rx x p, Sx = Ux * Lx * Xc, in which X = Lx * Xc.
  #p = ncol(X)
  # Calculate all individual components
  }
  gamma = rowMeans(S^3) # length r
  kappa = rowMeans(S^4 - 3) # length r

  # TU must be r by n
  JB = sum(alpha*gamma^2 + (1-alpha)*kappa^2)
  return(JB)

}


# Separate function for chordal distance of two vectors
chordal <- function(x, y){
  sum((tcrossprod(x)/sum(x^2) - tcrossprod(y)/sum(y^2))^2)
}

############
# BRisk: returns 4/5 of Jin's function
calculateJBofS <- function(S, alpha = 0.8){
  p = ncol(S)
  if (p < nrow(S)) warning("S should be r x p")
  # Calculate all individual components
  gamma = rowMeans(S^3) # length r
  kappa = rowMeans(S^4 - 3) # length r

  JB = sum(alpha*gamma^2 + (1-alpha)*kappa^2)

  return(JB)
}


#' Sign change for S matrix to image
#'
#' @param S t(S) px x r.
#' @param M another matrix.
#'
#' @return newS the original S * signskew, which is mean(S.column^3)
#' @export
#'
signchange = function(S,M=NULL) {
  # S: 59,412 x r_J
  signskew = sign(apply(S,2,function(x) mean(x^3)))
  newS = signskew*S          # t(signskew*t(S))
  if(!is.null(M)) newM = t(signskew*t(M))
  ifelse(is.null(M),return(newS),return(list(S=newS,M=newM)))
}

###############
# logistic <- function(xData, scale=sqrt(3)/pi, df=0) {
#   #maximizes likelihood given s then calculates gradient w.r.t. w.hat
#   #df is not used
#   xData = as.vector(xData)
#   list(Gs = -xData/scale - log(scale) - 2*log(1+exp(-xData/scale)), gs = -1/scale + 2*exp(-xData/scale)/(scale*(1+exp(-xData/scale))), gps = (2*exp(-2*xData/scale) - 2*exp(-xData/scale)*(1+exp(-xData/scale))) / (scale^2*(1+exp(-xData/scale))^2))
# }


logistic <- function(xData, scale=sqrt(3)/pi, df=0) {
  #maximizes likelihood given s then calculates gradient w.r.t. w.hat
  #df is not used
  xData = as.vector(xData)
  list(Gs = -xData/scale - log(scale) - 2*log(1+exp(-xData/scale)), gs = -1/scale + 2*exp(-xData/scale)/(scale*(1+exp(-xData/scale))), gps = - 2*exp(-xData/scale) / (scale^2*(1+exp(-xData/scale))^2))
}


jb.stat <- function(x, df=0) {
  n <- length(x)
  s <- sum(x^3)
  k <- sum(x^4)
  Gs <- x^3 * s / n^2 + (x^4 * k / n^2 + 9 / n - 6 * x^4 / n) / 4
  gs <- 6 * x^2 * s / n^2 + (8 * x^3 * (k / n - 3) / n) / 4
  gps <- 6 * (3 * x^4 + 2 * x * s) / n^2 + (24 * x^2 * (k / n - 3) / n + 32 * x^6 / n^2) / 4
  list(Gs = Gs, gs = gs, gps = gps)
}



#' tiltedgaussian
#'
#' @param xData input data
#' @param df degree freedom
#' @param B default value=100
#' @param ...
#'
#' @import gam
tiltedgaussian = function (xData, df = 8, B = 100, ...) {
  #This function is based on ProDenICA::GPois by Trevor Hastie
  #NOTE: Assumes data are zero mean.

  n <- length(xData)
  sd.x = sd(xData)
  rx <- c(min(xData)-0.1*sd.x, max(xData)+0.1*sd.x)
  xg <- seq(from = rx[1], to = rx[2], length = B)
  gaps <- diff(rx)/(B - 1)
  xcuts <- c(rx[1] - gaps/2, xg[-B] + gaps/2, rx[2] + gaps/2)
  #NOTE: I use the response variable that corresponds to the LCA paper.
  #This differs from the GPois algorithm in ProDenICA
  ys <- as.vector(table(cut(xData, xcuts)))/(gaps*n)
  pois.fit <- suppressWarnings(gam(ys ~ s(xg, df)+offset(dnorm(xg,log=TRUE)), family = poisson, ...))
  Gs <- predict(pois.fit) #log tilt function predicted at grid locations (note: predict on gam object can not be used to obtain derivatives)
  # the gam object with the predict function can not be used directly to obtain the derivatives
  # of the smoothing spline.
  # Here, we refit another iteration of the IRWLS algorithm used in gam:
  # Note: working residuals = (y - mu0)/mu0
  # weights = mu0
  # fitted(pois.fit) = mu0
  # predict(pois.fit) = eta0 = log(mu0)
  sGs = Gs #+ log(sum(dnorm(xg))/sum(fitted(pois.fit)))
  z0 <- sGs + residuals(pois.fit, type='working')
  pois.refit <- smooth.spline(x=xg, y=z0, w=fitted(pois.fit),df=df) #obtain the log tilt function in an object that can be used to obtain derivatives
  Gs <- predict(pois.refit, xData, deriv = 0)$y
  gs <- predict(pois.refit, xData, deriv = 1)$y
  gps <- predict(pois.refit, xData, deriv = 2)$y
  fGs <- function(x) predict(pois.refit,x,deriv=0)$y
  fgs <- function(x) predict(pois.refit,x,deriv=1)$y
  fgps <- function(x) predict(pois.refit,x,deriv=2)$y
  list(Gs = Gs, gs = gs, gps = gps, fGs = fGs, fgs=fgs, fgps=fgps)
}
#---------------------------------------------
#---------------------------------------
standardizeM <- function(M) {
  #M is d x p
  diag(diag(M%*%t(M))^(-1/2))%*%M
}



rtwonorm <- function(n, mean=c(0,5), sd=c(2/3,1), prob=0.5) {
  k <- rbinom(n,1,prob=prob)
  k*rnorm(n,mean[1],sd[1])+(1-k)*rnorm(n,mean[2],sd[2])
}

rmixnorm <- function(n, pars = list(mean=c(0,5), sd = c(2/3,1), prob=c(0.25,0.75))) {
  probs = pars[['prob']]
  means = pars[['mean']]
  sigmas = pars[['sd']]
  if(sum(probs)!=1) stop('Probabilities must sum to one')
  z = rmultinom(n=n, size=1, prob = probs)
  k = length(probs)
  #use rnorm recycling:
  x = rnorm(k*n,means,sigmas)
  dim(x) = c(k,n)
  apply(x*z,2,sum)
}

SimTwoNorms <- function(n.samples, distribution=c('mix-sub','mix-super'),snr,noisyICA=FALSE) {
  distribution = match.arg(distribution)
  if(distribution=='mix-sub') {
    mean=c(-1.7,1.7); sd=c(1,1); prob=0.75
  }
  if(distribution=='mix-super') {
    mean=c(0,5); sd=c(2/3,1); prob=0.95
  }
  sim.S <- rtwonorm(n=2*n.samples, mean=mean, sd=sd, prob=prob)
  dim(sim.S) <- c(n.samples,2)
  sim.M = myMixmat(5)
  sim.W = solve(sim.M)
  if(noisyICA) {
    sim.N <- matrix(rnorm(n=5*n.samples,mean=0,sd=1),nrow=n.samples,ncol=5)
  } else {
    sim.N <- matrix(rnorm(n=3*n.samples,mean=0,sd=1),nrow=n.samples,ncol=3)
  }

  sim.Ms = sim.M[1:2,]
  sim.Xs = sim.S%*%sim.Ms
  if(noisyICA) {
    sim.Mn = NULL
    sim.Xn <- sim.N
  } else {
    sim.Mn <- sim.M[3:5,]
    sim.Xn <- sim.N%*%sim.Mn
  }
  #alpha = 1/sqrt(mean(sim.Xs^2))
  alpha = 1/sd(as.vector(sim.Xs))
  sim.Xs = sim.Xs*alpha
  mean.S = apply(sim.S,2,mean)
  temp.S = scale(sim.S,center=TRUE,scale=FALSE)
  scalingMat = diag(apply(temp.S,2,sd))
  scaled.sim.S = temp.S%*%solve(scalingMat)
  scaled.sim.Ms = sqrt(snr)*alpha*scalingMat%*%sim.Ms
  #sim.Xn = sim.Xn/sqrt(mean(sim.Xn^2))
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn))
  sim.Xs = sqrt(snr)*sim.Xs #equivalent to scaled.sim.S%*%(alpha*sqrt(snr)*scalingMat%*%sim.Ms)+alpha*sqrt(snr)*matrix(mean.S,nrow=n.samples,ncol=2,byrow=TRUE)%*%sim.Ms
  #since we only recover scaled.sim.S, "true Ms" for, e.g., IFA, is defined as in scaled.sim.Ms
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)
  return(list(sim.S = sim.S, sim.N = sim.N, sim.Ms = sim.Ms, sim.Mn = sim.Mn, sim.X=sim.X, scaled.sim.S = scale(sim.S),scaled.sim.Ms = scaled.sim.Ms,scaled.sim.X = scale(sim.X), whitened.sim.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
}


#' Standardization with double centered and column scaling
#'
#' @param data input matrix with n x px.
#' @param dif.tol the value for the threshold of scaling
#' @param max.iter default value = 10
#'
#' @return standardized matrix with n x px.
#' @export
#'
#' @examples
#' spmwm = 3*matrix(rnorm(100000),nrow=100)+1
#' dim(spmwm)
#' apply(spmwm,1,mean) # we want these to be 0
#' apply(spmwm,2,mean) # we want these to be 0
#' apply(spmwm,2,sd) # we want each of these variances to be 1
#'
#' spmwm_cp=standard(spmwm)
#' max(abs(apply(spmwm_cp,1,mean)))
#' max(abs(apply(spmwm_cp,2,mean)))
#' max(abs(apply(spmwm_cp,2,sd)-1))
standard <- function(data,dif.tol=1e-03,max.iter=10){
  row_mean_max = max(abs(apply(data,1,mean)))
  col_mean_max = max(abs(apply(data,2,mean)))
  col_sd_max= max(abs(apply(data,2,sd)-1))
  n=0
  while(n<=max.iter & max(row_mean_max,col_mean_max,col_mean_max)>= dif.tol) {
    data=scale(data) # centering and scaling for each column
    data=t(scale(t(data),center = T,scale = F))
    n=n+1
  }
  return(data)
}

