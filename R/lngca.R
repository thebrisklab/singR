#----------------------------------
# Benjamin Risk and Irina Gaynanova
# Contact: brisk@emory.edu
# Functions supporting SING
#------------------------------

#' Decompose the original data through LNGCA method.
#'
#' Implements the methods of linear non-Gaussian component analysis (LNGCA) and likelihood component analysis (when using a density, e.g., tilted Gaussian) from the \href{https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1407772}{LNGCA paper}
#'
#' @param xData the original dataset for decomposition, matrix of n x px.
#' @param n.comp the number of components to be estimated.
#' @param Ux.list list of user specified initial values for Ux. If null, will generate random orthogonal matrices. See restarts.pbyd and restarts.dbyd
#' @param whiten whitening method. Defaults to "svd" which uses the n left eigenvectors divided by sqrt(px-1). Optionally uses the square root of the n x n "precision" matrix.
#' @param maxit max iteration, defalut = 1000
#' @param eps default = 1e-06
#' @param verbose default = FALSE
#' @param restarts.pbyd default = 0. Generates p x d random orthogonal matrices. Use a large number for large datasets. Note: it is recommended that you run lngca twice with different seeds and compare the results, which should be similar when a sufficient number of restarts is used. In practice, stability with large datasets and a large number of components can be challenging.
#' @param restarts.dbyd default = 0. These are d x d initial matrices padded with zeros, which results in initializations from the principal subspace. Can speed up convergence but may miss low variance non-Gaussian components.
#' @param distribution distribution methods with default to tilted Gaussian. "logistic" is similar to infomax ICA, JB is capable of capture super and sub Gaussian distribution while being faster than tilted Gaussian. (tilted Gaussian tends to be most accurate, but computationally much slower.)
#' @param density return the estimated tilted Gaussian density? default = FALSE
#' @param out.all default = FALSE
#' @param orth.method default = 'svd'. Method to generate random initial matrices. See [gen.inits()]
#' @param df default = 0, df of the spline used in fitting the non-parametric density. use df=8 or so for tilted gaussian. set df=0 for JB and logistic.
#' @param stand whether to standardize the data to have row and column means equal to 0 and the row standard deviation equal to 1 (i.e., all variables on same scale). Often used when combined with singR for data integration.
#' @param ... other arguments to tiltedgaussian estimation
#'
#' @return Function outputs a list including the following:
#' \describe{
#'       \item{\code{U}}{matrix rx x n, part of the expression that Ax = Ux x Lx and Ax x Xc = Sx, where Lx is the whitener matrix.}
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
#' \donttest{
#' #get simulation data
#' data(exampledata)
#' data=exampledata
#'
#' # To get n.comp value, we can use NG_number function.
#'
#' # use JB statistic as the measure of nongaussianity to run lngca with df=0
#' estX_JB = lngca(xData = data$dX, n.comp = 4,
#'  whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB',df=0)
#'
#' # use the tiltedgaussian distribution to run lngca with df=8. This takes a long time:
#' estX_tilt = lngca(xData = data$dX, n.comp = 4,
#'  whiten = 'sqrtprec', restarts.pbyd = 20, distribution='tiltedgaussian',df=8)
#'
#' # true non-gaussian component of Sx, include individual and joint components
#' trueSx = rbind(data$sjX,data$siX)
#'
#' # use pmse to compare the difference of the two methods
#' pmse(S1 = t(trueSx),S2=t(estX_JB$S),standardize = TRUE)
#' pmse(S1 = t(trueSx),S2=t(estX_tilt$S),standardize = TRUE)
#'
#' # the lngca using tiltedgaussian tends to be more accurate
#' # with smaller pmse value, but takes longer to run.
#'}
#'
lngca <- function(xData, n.comp = NULL, Ux.list = NULL, whiten = c('sqrtprec','eigenvec','none'), maxit = 1000, eps = 1e-06, verbose = FALSE, restarts.pbyd = 0, restarts.dbyd = 0, distribution=c('JB','tiltedgaussian','logistic'), density=FALSE, out.all=FALSE, orth.method=c('svd','givens'),df=0,stand=FALSE,...) {



  #note: small changes from mlcaFP from the JASA paper:
  # 1) output Mhat.
  # 2) order by skewness, with option for n.comp=1

  #former option: deflation optimization, not currently implemented
  #alg.typ = c('symmetric','deflation'),
  #alg.typ = match.arg(alg.typ)


  whiten=match.arg(arg = NULL,choices = whiten)
  if (nrow(xData)>ncol(xData)) warning("LNGCA may produce unexpected results with n>p.")
  if (nrow(xData)>ncol(xData) & whiten != "sqrtprec") stop("with n > p, use whiten = 'sqrtprec'. note: LNGCA maximizes non-Gaussianity across p_x.")

  # Liangkang: Change the input name Ux.list to W.list
  W.list=Ux.list

  alg.typ = 'symmetric'


  distribution = match.arg(arg = NULL,choices = distribution)



  if(distribution=="JB" || distribution == "logistic"){
    df=0
    if(distribution=='logistic') Gfunc = logistic
    if(distribution=='JB') Gfunc = jb.stat
    }

  if(restarts.dbyd>0 && whiten!='eigenvec') stop('Use whiten=eigenvec with restarts.dbyd')
  ## whiten:

  #require(multidcov)
  if(distribution=='tiltedgaussian') {
    Gfunc = tiltedgaussian
    if (df==0) {df = 8 }
    #require(ProDenICA)
  }


  #if(distribution=='tiltedgaussian' && df==0) stop('df must be greater than 0 for tiltedgaussian')
  #if(distribution=='logistic'  && df>0) stop('df should be set to zero when using logistic')

  #if(distribution=='JB'  && df>0) stop('df should be set to zero when using JB')

  #if(!is.null(W.list) && class(W.list)!='list') stop('W.list must be a list')
  if(!is.null(W.list) && !is.list(W.list)) stop('W.list must be a list')
  if(length(W.list) && (restarts.pbyd || restarts.dbyd)) stop('restarts.pbyd and restarts.dbyd must be equal to zero when supplying W.list')

  xData=t(xData) # transform the input dimension from n x px to px x n.

  if(is.null(n.comp)) {n.comp=ncol(xData)}

  orth.method= match.arg(arg = NULL,choices = orth.method)
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
      loglik.v[k] = out.list[[k]]$loglik
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
  # IGAY: re-order Ws based on loglik as well
  out$Ws = out$Ws[, order(loglik.maxS,decreasing=TRUE)]
  # Liangkang: change the out subset name from Ws to U, and transpose to r x n dimension
  out$U=t(out$Ws)
  if("Ws" %in% names(out)) out <- out[ - which(names(out) == "Ws")]
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
#' @param sData S rx x px
#' @param xData dX n x px
#' @param intercept default = TRUE
#'
#' @return a matrix Mx, dimension n x rx.
#' @export
#'
est.M.ols <- function(sData,xData,intercept=TRUE) {
  sData=t(sData)
  xData=t(xData)
  if(intercept) t(stats::coef(stats::lm(xData~sData))[-1,]) else t(stats::coef(stats::lm(xData~sData-1)))
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
  pois.fit <- suppressWarnings(gam(ys ~ s(xg, df)+offset(dnorm(xg,log=TRUE)), family = stats::poisson, ...))
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
