#----------------------------------
# Benjamin Risk and Irina Gaynanova
# Contact: brisk@emory.edu
# Functions supporting SING
#------------------------------

# IGAY: edited on June 5th, 2019 to correct ordering of Ws
# BRisk: edited on 17 July 2019 to estimate M for xData, not for whitened data
#' Decomposite the original data through LNGCA method
#'
#' @param xData the original dataset for decomposition, matrix of px x n
#' @param n.comp the number of components need to be estimated.
#' @param W.list list
#' @param whiten whiten method with default
#' @param maxit max iteration, defalut = 1000
#' @param eps default = 1e-06
#' @param verbose default = FALSE
#' @param restarts.pbyd default = 0
#' @param restarts.dbyd default = 0
#' @param distribution distribution methods with default
#' @param density default = FALSE
#' @param out.all default = FALSE
#' @param orth.method default = 'svd'
#' @param reinit.max.comp default = FALSE
#' @param max.comp logical variable that estiamtes the max number of non-gaussian components.
#' @param df default = 0
#' @param stand whether to standard the data
#' @param ... ellipsis
#'
#' @return Function outputs a list including the following:
#' \describe{
#'       \item{\code{U}}{ matrix rx x n, part of the expression that Ax = Ux x Lx and Ax x Xc = Sx, where Lx is the whitener matrix.}
#'       \item{\code{loglik}}{the value of log-likelihood in the lngca method.}
#'       \item{\code{S}}{the variable loading matrix r x px, each row is a component, which can be used to measure nongaussianity}
#'       \item{\code{df}}{egree of freedom.}
#'       \item{\code{distribution}}{the method used for data decomposition.}
#'       \item{\code{whitener}}{A symmetric whitening matrix n x n from dX, the same with  whitenerXA = est.sigmaXA \%^\% -0.5}
#'       \item{\code{M}}{Mx Mtrix with n x rx.}
#'       \item{\code{nongaussianity}}{the nongaussianity score for each component saved in S matrix.}
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' #get simulation data
#' data(exampledata)
#' data=exampledata
#'
#' # To get n.comp value, we can use NG_number function.
#'
#' # use JB statistic as the measure of nongaussianity to run lngca with df=0
#' estX_JB = lngca(xData = t(data$dX), n.comp = 4, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB',df=0)
#'
#' # use the tiltedgaussian distribution to run lngca with df=3
#' estX_tilt = lngca(xData = t(data$dX), n.comp = 4, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='tiltedgaussian',df=3)
#'
#' # true non-gaussian component of Sx, includ individual and joint components
#' trueSx = rbind(data$sjX,data$siX)
#'
#' # use frobICA to compare the difference of the two methods
#' frobICA(S1 = t(trueSx),S2=t(estX_JB$S),standardize = TRUE) #0.2341096
#' frobICA(S1 = t(trueSx),S2=t(estX_tilt$S),standardize = TRUE) #0.1824922
#'
#' # the lngca using tiltedgaussian is more accurate with smaller frobICA value. But it will cost more time.
#'}
#'
lngca <- function(xData, n.comp = ncol(xData), W.list = NULL, whiten = c('eigenvec','sqrtprec','none'), maxit = 1000, eps = 1e-06, verbose = FALSE, restarts.pbyd = 0, restarts.dbyd = 0, distribution=c('tiltedgaussian','logistic','JB'), density=FALSE, out.all=FALSE, orth.method=c('svd','givens'), reinit.max.comp = FALSE, max.comp = FALSE, df=0,stand=FALSE,...) {

  #note: small changes from mlcaFP from the JASA paper:
  # 1) output Mhat.
  # 2) order by skewness, with option for n.comp=1

  #former option:
  #alg.typ = c('symmetric','deflation'),
  #alg.typ = match.arg(alg.typ)
  alg.typ = 'symmetric'


  distribution = match.arg(distribution)
  whiten=match.arg(whiten)

  if(restarts.dbyd>0 && whiten!='eigenvec') stop('Use whiten=eigenvec with restarts.dbyd')
  ## whiten:

  if(max.comp) { #if statement evaluates to true for all max.comp!=0
    s.comp = n.comp
    n.comp = max.comp
  }
  if(max.comp=='TRUE') stop('max.comp should be an integer or FALSE')
  if(reinit.max.comp && max.comp==FALSE) stop('Can not reinitialize from max.comp solution if max.comp==FALSE')
  if(reinit.max.comp && alg.typ=='deflation') stop('reinit.max.comp not yet written for deflation algorithm')

  #require(multidcov)
  if(distribution=='tiltedgaussian') {
    Gfunc = tiltedgaussian
    #require(ProDenICA)
  }


  if(distribution=='tiltedgaussian' && df==0) stop('df must be greater than 0 for tiltedgaussian')
  if(distribution=='logistic'  && df>0) stop('df should be set to zero when using logistic')
  if(distribution=='logistic') Gfunc = logistic
  if(distribution=='JB'  && df>0) stop('df should be set to zero when using JB')
  if(distribution=='JB') Gfunc = jb.stat
  if(!is.null(W.list) & class(W.list)!='list') stop('W.list must be a list')
  if(length(W.list) && (restarts.pbyd || restarts.dbyd)) stop('restarts.pbyd and restarts.dbyd must be equal to zero when supplying W.list')

  orth.method= match.arg(orth.method)
  p = ncol(xData)
  nRow = nrow(xData)
  d = n.comp

  # center xData such that ones%*%xData = 0
  # Liangkang fix: use the standard function to double center the data.
  if(stand){
    xData = t(standard(t(xData))) #standard input N x px matrix.
  } else {
    xData = scale(xData, center=TRUE, scale=FALSE)# minus the mean to center the xData
  }



  if (d > p) stop('d must be less than or equal to p')
  if (whiten=='eigenvec') {
    # Use whitener=='eigenvec' so that restarts.dbyd initiates from the
    # span of the first d eigenvectors.
    temp = whitener(X = xData,n.comp = p)
    zData = temp$Z
    whitener = temp$whitener
    rm(temp)
  }  else if (whiten=='sqrtprec') {
    est.sigma = stats::cov(xData)  ## Use eigenvalue decomposition rather than SVD.
    evd.sigma = svd(est.sigma)
    whitener = evd.sigma$u%*%diag(evd.sigma$d^(-1/2))%*%t(evd.sigma$u)
    zData = xData%*%whitener
  }
  else {
    # no whitening:
    zData = xData
    whitener = diag(p)
  }
  # warning('TO DO: Look at whitening methods and check for inconsistent options')
  if (is.null(W.list)) {
    if(restarts.pbyd) W.list = gen.inits(p=p,d=d,runs=restarts.pbyd,orth.method=orth.method)
    if(restarts.dbyd) {
      W.temp = gen.inits(p=d,d=d,runs=restarts.dbyd,orth.method=orth.method)
      #pad with zeros:
      zeros = matrix(0,p-d,d)
      W.temp = lapply(W.temp,FUN = function(x) rbind(x,zeros))
      W.list = c(W.list,W.temp)
    }
  }
  ## If restarts.pbyd and restarts.dbyd both equal zero:
  if (is.null(W.list)) W.list = gen.inits(p=p,d=d,runs=1,orth.method=orth.method)
  runs = length(W.list)
  out.list = NULL
  loglik.v = numeric(runs)
  for(k in 1:runs) {
    W0 = as.matrix(W.list[[k]])
    if(alg.typ == 'symmetric') {
      out.list[[k]] = lca.par(xData=zData,W0=W0,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=n.comp,df=df, ...)
      out.list[[k]]$df = df
    }
    #if(alg.typ == 'deflation') {
     # out.list[[k]] = lca.def(xData=zData,W0=W0,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=n.comp,df=df,...)
     #out.list[[k]]$df = df
    #}
    if(max.comp) {
      flist0 = list()
      for (j in 1:d) flist0[[j]] <- Gfunc(out.list[[k]]$S[, j], ...)
      loglik.S <- apply(sapply(flist0, "[[", "Gs"),2,sum)
      orderedLL = order(loglik.S,decreasing=TRUE)
      out.list[[k]]$S = out.list[[k]]$S[,orderedLL[1:s.comp]]
      out.list[[k]]$Ws = out.list[[k]]$Ws[,orderedLL[1:s.comp]]
      out.list[[k]]$loglik = sum(loglik.S[orderedLL[1:s.comp]])
      loglik.v[k] = out.list[[k]]$loglik
    } else {
      loglik.v[k] = out.list[[k]]$loglik
    }
  }
  for(i in 1:runs){
    out.list[[i]]$distribution=distribution
    out.list[[i]]$whitener = whitener
  }
  out = out.list[[which.max(loglik.v)]]

  maxS = out$S
  loglik.maxS = numeric(n.comp)

  #order by negentropy:
  densities.maxS = apply(maxS,2,Gfunc)

  # note: sum here leads to correct jb statistic
  # IGAY: this is the first place where ordering is applied to maxS (out$S)
  for(i in 1:n.comp) loglik.maxS[i] = sum(densities.maxS[[i]][['Gs']])
  o.maxS = maxS[,order(loglik.maxS,decreasing=TRUE)]
  rm(maxS)

  #force skewness to be positive:
  # added option for one component on 6 May 2019
  if(n.comp>1) {
    Shat = rightskew(S = o.maxS,order.skew = FALSE)
  } else {
    Shat = sign(mean(o.maxS^3))*o.maxS
  }
  # equivalent to est.M.ols for since 0 has been centered above
  # BRisk note: This gives the mixing matrix for the singular vectors, not for the data
  # so it is M t(V) D^{-1},

  Mhat <- t(Shat)%*%xData

  # IGAY: this is where out is reassigned for S and M, but not Ws with right ordering
  out$S = Shat
  out$M = Mhat
  # IGAY: fix the ordering of Ws based on loglik as well, fix on June 5th, 2019
  out$Ws = out$Ws[, order(loglik.maxS,decreasing=TRUE)]
  # Liangkang: change the out subset name from Ws to U, and transpose to r x n dimension
  out$U=t(out$Ws)
  out=within(out,rm(Ws))
  # Liangkang: transpose the M to n x r dimension.
  out$M=t(out$M)
  # Liangkang: transpose the S to r x px dimension.
  out$S=t(out$S)

  out$nongaussianity = sort(loglik.maxS,decreasing = TRUE)

  if(out.all==TRUE) {
    out.list[[which.max(loglik.v)]] = out
    out.list
  }  else {
    out
  }
}
#---------------------------------
#----------------------------
#
#' Estimate mixing matrix from estimates of components
#'
#' @param sData t(S) px x rx
#' @param xData t(dX) px x n
#' @param intercept default = TRUE
#'
#' @return a matrix t(Mx) n.comp x n
#' @export
#'
est.M.ols <- function(sData,xData,intercept=TRUE) {
  if(intercept) stats::coef(stats::lm(xData~sData))[-1,] else stats::coef(stats::lm(xData~sData-1))
}
# NOTE: for centered X, equivalent to xData %*% sData/(px-1)
# Correced by Liankang
# NOTE: for centered X, equivalent to t(sData)%*%xData/(V-1)



#' tiltedgaussian
#'
#' @param xData input data
#' @param df degree freedom
#' @param B default value=100
#' @param ... ellipsis
#'
#' @import gam
tiltedgaussian = function (xData, df = 8, B = 100, ...) {
  #This function is based on ProDenICA::GPois by Trevor Hastie
  #NOTE: Assumes data are zero mean.

  n <- length(xData)
  sd.x = stats::sd(xData)
  rx <- c(min(xData)-0.1*sd.x, max(xData)+0.1*sd.x)
  xg <- seq(from = rx[1], to = rx[2], length = B)
  gaps <- diff(rx)/(B - 1)
  xcuts <- c(rx[1] - gaps/2, xg[-B] + gaps/2, rx[2] + gaps/2)
  #NOTE: I use the response variable that corresponds to the LCA paper.
  #This differs from the GPois algorithm in ProDenICA
  ys <- as.vector(table(cut(xData, xcuts)))/(gaps*n)
  pois.fit <- suppressWarnings(gam(ys ~ s(xg, df)+offset(dnorm(xg,log=TRUE)), family = poisson, ...))
  Gs <- stats::predict(pois.fit) #log tilt function predicted at grid locations (note: predict on gam object can not be used to obtain derivatives)
  # the gam object with the predict function can not be used directly to obtain the derivatives
  # of the smoothing spline.
  # Here, we refit another iteration of the IRWLS algorithm used in gam:
  # Note: working residuals = (y - mu0)/mu0
  # weights = mu0
  # fitted(pois.fit) = mu0
  # predict(pois.fit) = eta0 = log(mu0)
  sGs = Gs #+ log(sum(dnorm(xg))/sum(fitted(pois.fit)))
  z0 <- sGs + stats::residuals(pois.fit, type='working')
  pois.refit <- stats::smooth.spline(x=xg, y=z0, w=stats::fitted(pois.fit),df=df) #obtain the log tilt function in an object that can be used to obtain derivatives
  Gs <- stats::predict(pois.refit, xData, deriv = 0)$y
  gs <- stats::predict(pois.refit, xData, deriv = 1)$y
  gps <- stats::predict(pois.refit, xData, deriv = 2)$y
  fGs <- function(x) stats::predict(pois.refit,x,deriv=0)$y
  fgs <- function(x) stats::predict(pois.refit,x,deriv=1)$y
  fgps <- function(x) stats::predict(pois.refit,x,deriv=2)$y
  list(Gs = Gs, gs = gs, gps = gps, fGs = fGs, fgs=fgs, fgps=fgps)
}
#---------------------------------------------
