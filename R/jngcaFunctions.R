#----------------------------------
# Benjamin Risk and Irina Gaynanova
# Contact: brisk@emory.edu
# Functions supporting SING
#------------------------------


# This function for easy parallelization needs to be revised to work better on other systems. Right now it involves my user-specific directory structure;
# some issues with passing variables to cluster

lngca_multicore <- function(xData, n.comp = ncol(xData), ncores=NULL, W.list = NULL, whiten = c('eigenvec','none'), restarts=NULL, distribution=c('tiltedgaussian','logistic','JB'), df=0, initFOBI = FALSE, scratch="~/temp",keepall = FALSE,sourcefile='~/Dropbox/JINGCA/Programs/Functions/jngcaFunctions.R',...) {
  warning("If running mlcaFP_parallelized, you need to manually change the path that contains jngca_functions in order to load the functions onto the workers. Sorry for the inconvenience. See comments in the program.")


# require(parallel)
  if (is.null(ncores)) {
    ncores = detectCores()
  }

  system(paste('mkdir',scratch))  # create a temporary directory in the scratch adress#

  cl = makeSOCKcluster(rep("localhost",ncores)) # what is the make SOCKcluster function #

  if(is.null(restarts)) restarts=ncores  # default restars = ncores #

  source(sourcefile)    # load functions from source file #
  #clusterExport(cl,"sourcefile")
  #clusterEvalQ(cl,sourcefile)

  # can't get this to work with path name variable!!!
  # Need to edit this manually to get it to work:
  invisible(clusterEvalQ(cl,library(singR)))  # change the print mode #

  invisible(clusterEvalQ(cl,.libPaths('~/Rlibs')))   # what is the clusterEvalQ function #
  .libPaths('~/Rlibs')

  n.Voxels = nrow(xData)
  xData <- scale(xData, center=TRUE, scale=FALSE)  # normalization with sample mean and sample std #
  whiten=match.arg(whiten) # get the first value in whiten #

  distribution = match.arg(distribution)  # first value in distribution #
  # perform whitening once:
  if (whiten=='eigenvec') {
    zdata = svd(xData)$u*sqrt(n.Voxels) # x$u means left matrix in single value decomposition #
  } else {
    zdata = xData
  }

  n.Time = ncol(zdata)

  #Generate initializations from the principal subspace:
  W.temp = gen.inits(p=n.comp,d=n.comp,runs=ceiling(restarts/2)-1,orth.method='svd') # generate runs p*d orthodox matrix with svd method in w.list #
  zeros = matrix(0,n.Time-n.comp,n.comp)
  W.temp = lapply(W.temp,FUN = function(x) rbind(x,zeros)) #pad with zeros # the W.temp is n.Time*n.comp #
  W.temp[[1]] = diag(n.Time)[,1:n.comp] #one initialization from PCA solution. # n.Time * n.comp matrix #

  #Generate initializations from the column space of the data:
  W.list=gen.inits(p=n.Time,d=n.comp,runs=ceiling(restarts/2),orth.method='svd')
  W.list = c(W.temp,W.list)

  #Generate an initialization from the FOBI solution:
  if (initFOBI) {
    require(JADE)
    estFOBI = FOBI(zdata)
    W.temp = list(t(estFOBI[['W']])[,1:n.comp])
    W.list = c(W.temp,W.list)
    rm(estFOBI)  # delete estFOBI #
  }

  newruns = length(W.list)
  mymlcaFP <- function(x,W.list,xData,n.comp,maxit=maxit,distribution,df) {
    est = mlcaFP(xData = xData, n.comp = n.comp, W.list=list(W.list[[x]]), whiten="none", maxit=300, distribution=distribution, df=df)

    save(est,file=paste0(scratch,'/est',x,'.RData'))
  }

  # estimate on whitened data:
  clusterApply(cl=cl,x=as.list(1:newruns),fun=mymlcaFP,W.list=W.list,xData=zdata,n.comp=n.comp,distribution=distribution,df=df)
  stopCluster(cl)

  ###################
  # 2. Find argmax
  ##Get list of estimates:
  setwd(scratch)
  list.est = NULL
  for (x in 1:newruns) {
    list.est = c(list.est,paste0("est",x,'.RData'))
  }

  loglik.v = numeric(newruns)

  ##Find best likelihood:
  for(i in 1:newruns) {
    load(list.est[i])
    loglik.v[i] = est[['loglik']]
  }

  load(list.est[which.max(loglik.v)])
  # # clean up intermediate files:
   if (!keepall) {
     for (i in 1:newruns) {
     tmpfile = paste0("est",i,".RData")
     file.remove(file.path(scratch,tmpfile))
     }
   }
  #

  # equivalent to est.M.ols for since Xdata has been centered above
   Mhat <- t(est$S)%*%xData

  # 19 April 2019: changed named of Shat to S and Mhat to M
  # 16 July 2019: Deleted ordering since this was done in mlcaFP
#  return(list(S=Shat,M=Mhat,loglik.maxS=sort(loglik.maxS,decreasing=TRUE),nongaussianity=nongaussianity,loglik.allinits=loglik.v))
  return(list(S=est$S,M=Mhat,Ws=est$Ws,nongaussianity=est$nongaussianity,loglik.allinits=loglik.v,whitener=NULL))

  }



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
#' @param irlba whether require irlba package.
#' @param ... ellipsis
#'
#' @return Function outputs a list including the following:
#' \describe{
#'       \item{\code{Ws}}{t(Ux) matrix n x n.comp, part of the expression that Ax = Ux x Lx and Ax x Xc = Sx, where Lx is the whitener matrix.}
#'       \item{\code{loglik}}{the value of log-likelihood in the lngca method.}
#'       \item{\code{S}}{the variable loading matrix px x n.comp, each column is a component, which can be used to measure nongaussianity}
#'       \item{\code{df}}{egree of freedom.}
#'       \item{\code{distribution}}{the method used for data decomposition.}
#'       \item{\code{whitener}}{A symmetric whitening matrix n x n from dX, the same with  whitenerXA = est.sigmaXA \%^\% -0.5}
#'       \item{\code{M}}{Mtrix with n.comp x n.}
#'       \item{\code{nongaussianity}}{the nongaussianity score for each component saved in S matrix.}
#' }
#'
#' @export
#'
#' @import irlba
#' @import ProDenICA
#'
#'
lngca <- function(xData, n.comp = ncol(xData), W.list = NULL, whiten = c('eigenvec','sqrtprec','none'), maxit = 1000, eps = 1e-06, verbose = FALSE, restarts.pbyd = 0, restarts.dbyd = 0, distribution=c('tiltedgaussian','logistic','JB'), density=FALSE, out.all=FALSE, orth.method=c('svd','givens'), reinit.max.comp = FALSE, max.comp = FALSE, df=0, irlba=FALSE,...) {

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
    xData <- scale(xData, center=TRUE, scale=FALSE)  # minus the mean to center the xData #
    if (d > p) stop('d must be less than or equal to p')
    if (whiten=='eigenvec') {
      # Use whitener=='eigenvec' so that restarts.dbyd initiates from the
      # span of the first d eigenvectors.
      temp = whitener(X = xData,n.comp = p,irlba=irlba)
      zData = temp$Z
      whitener = temp$whitener
      rm(temp)
      }  else if (whiten=='sqrtprec') {
         est.sigma = cov(xData)  ## Use eigenvalue decomposition rather than SVD.
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
    if(alg.typ == 'deflation') {
      out.list[[k]] = lca.def(xData=zData,W0=W0,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=n.comp,df=df,...)
    out.list[[k]]$df = df
    }
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

  out$nongaussianity = sort(loglik.maxS,decreasing = TRUE)

  if(out.all==TRUE) {
    out.list[[which.max(loglik.v)]] = out
    out.list
  }  else {
    out
}
}
#---------------------------------

#-------------
orthonormalize <- function(xk,X,k) {
  #Gram-Schmidt orthogonalization
  #assumes columns of X have length equal to 1
  if(k!=1) {
    t <- numeric(length(xk))
    for (u in 1:(k-1)) {
      a <- sum(xk * X[,u])
      t <- t + a * X[,u]
    }
    xk <- xk - t
    }
  xk / sqrt(crossprod(xk))
}

genWfromWs <- function(Ws) {
  d = ncol(Ws)
  p = nrow(Ws)
  tempW = cbind(Ws,diag(p)[,(d+1):p])
  for(k in (d+1):p) {
    oldWk = tempW[,k]
    tempWk = tempW[,k]
    for(j in 1:(k-1)) {
      tempWj = tempW[,j]
      tempWk = tempWk - tempWj * crossprod(tempWj,oldWk)/crossprod(tempWj,tempWj)
  }
  tempW[,k] = tempWk/sqrt(crossprod(tempWk))
  }
  tempW
}

temp.orthogonalize <- function(V,W) {
  #orthogonalizes the vector V to all columns in W
  #and returns cbind(W,orthV)
  oldWk = V
  tempWk = V
  tempW=cbind(W,V)
  k=ncol(W)+1
  for(j in 1:(k-1)) {
      tempWj = tempW[,j]
      tempWk = tempWk - tempWj * crossprod(tempWj,oldWk)/crossprod(tempWj)
    }
  tempW[,k] = tempWk/sqrt(crossprod(tempWk))
  tempW
}


#------------------------------------------------
# symmetric algorithm:
lca.par <- function(xData,W0,Gfunc,maxit,verbose,density,eps,n.comp,df,...) {
  W0 = as.matrix(W0)
  d = ncol(W0)
  if(n.comp!=d) stop('W0 needs to be p x d')
  p = ncol(xData)
  nRow = nrow(xData)
  s <- xData %*% W0
  flist <- as.list(1:d)
  ##  Densities of likelihood components:
  for (j in 1:d) flist[[j]] <- Gfunc(s[, j], df=df,...)
  flist0 <- flist
  crit0 <- mean(sapply(flist0, "[[", "Gs"))
  nit <- 0
  nw <- 10
  repeat {
    nit <- nit + 1
    gS <- sapply(flist0, "[[", "gs")
    gpS <- sapply(flist0, "[[", "gps")
    #t1 <- t(xData) %*% gS/nRow
    t1 <- crossprod(xData,gS)/nRow
    t2 <- apply(gpS, 2, mean)
    if(d>1) W1 <- t1 - W0%*%diag(t2) else W1 <- t1 - W0*t2
    W1 <- orthogonalize(W1)
    if(d>1) nw <- frobICA(t(W0), t(W1))^2 else nw <- mean((W0-W1)^2) #Uses a measure that works for non-square matrices -- MSE. The measure is defined for M so here we use transpose of W.
    W0 <- W1
    s <- xData %*% W0
    for (j in 1:d) flist0[[j]] <- Gfunc(s[, j], df=df, ...)
    crit0 <- mean(sapply(flist0, "[[", "Gs"))
    if (verbose) cat("Iter", nit, "G", crit0, "Delta", nw, "\n")
    if ((nit >= maxit)) {
      warning('Max iter reached')
      break
    }
    if (nw < eps) break
  }
  out = list(Ws = W0, loglik = d*nRow*crit0, S = s)
  if(density) out$density = lapply(flist0, "[[", "density")
  out
}

#--------------------------------------
myMixmat <-  function (p = 2) {
  a <- matrix(rnorm(p * p), p, p)
  sa <- svd(a)
  d <- sort(runif(p,min=1,max=10))
  mat <- sa$u %*% (sa$v * d)
  attr(mat, "condition") <- d[p]/d[1]
  mat
}

#------------------------
#---------------------------------------
standardizeM <- function(M) {
  #M is d x p
  diag(diag(M%*%t(M))^(-1/2))%*%M
}



#----------------------------
#
#' Estimate mixing matrix from estimates of components
#'
#' @param sData t(S) px x n.comp
#' @param xData t(dX) px x n
#' @param intercept default = TRUE
#'
#' @return a matrix t(Mx) n.comp x n
#' @export
#'
est.M.ols <- function(sData,xData,intercept=TRUE) {
  if(intercept) coef(lm(xData~sData))[-1,] else coef(lm(xData~sData-1))
}
# NOTE: for centered X, equivalent to xData %*% sData/(px-1)
# Correced by Liankang
# NOTE: for centered X, equivalent to t(sData)%*%xData/(V-1)

#-----------------------------------
# order by likelihood
# option for positive skewness
order.likelihood <- function(S,positive.skew=TRUE,distribution=c('logistic','tiltedgaussian','logcosh','JB'),out.loglik=FALSE,...) {
  distribution = match.arg(distribution)
  nObs = nrow(S)
  d = ncol(S)
  if(distribution=='tiltedgaussian') Gfunc = tiltedgaussian
  if(distribution=='logistic') Gfunc = logistic
  if(distribution=='logcosh') Gfunc = ProDenICA::G1
  if(distribution=='JB') Gfunc = jb.stat
  if(positive.skew) {
    skewness <- function(x, n = nObs) (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
    skew = apply(S, 2, skewness)
    sign = -1 * (skew < 0) + 1 * (skew > 0)
    S = S %*% diag(sign)
  }
  flist0 = list()
  for (j in 1:d) flist0[[j]] <- Gfunc(S[, j], ...)
  loglik.S <- apply(sapply(flist0, "[[", "Gs"),2,sum)
  orderedLL = order(loglik.S,decreasing=TRUE)
  S = S[,orderedLL]
  if(out.loglik) return(list(S=S,loglik=sort(loglik.S,decreasing=TRUE))) else S
}

marginal.likelihoods <- function(S,distribution=c('logistic','tiltedgaussian','logcosh','GPois','JB'),...)
{
  distribution = match.arg(distribution)
  if(distribution=='tiltedgaussian') Gfunc = tiltedgaussian
  if(distribution=='logistic') Gfunc = logistic
  if(distribution=='logcosh') Gfunc = ProDenICA::G1
  if(distribution=='GPois') Gfunc = ProDenICA::GPois
  if(distribution=='JB') Gfunc = jb.stat
  d = ncol(S)
  flist0 = list()
  for (j in 1:d) flist0[[j]] <- Gfunc(S[, j], ...)
  apply(sapply(flist0, "[[", "Gs"),2,sum)
}




#-------------------------------------
# Match mixing matrices:
# This function does not require M to be square:
#' Match mixing matrices
#'
#' @param M1 Subject score 1
#' @param M2 Subject score 2
#' @param S1 Loading 1
#' @param S2 Loading 2
#' @param standardize whether to standardize
#'
#' @export
#' @import clue
frobICA<-function(M1=NULL,M2=NULL,S1=NULL,S2=NULL,standardize=FALSE) {
  #MODEL: X = S M + E, so M is d x p
  #standardize: if standardize==TRUE, then standardizes rows of M1 and M2
  #to have unit norm; if using S1 and S2, standardizes columns to have unit variance.
  #standardize=TRUE makes the measure scale invariant.


  tfun = function(x) all(x==0)
  if(is.null(M1) && is.null(M2) && is.null(S1) && is.null(S2)) stop("need to supply either M1 and M2 or S1 and S2")
  if(!is.null(M1) && !is.null(M2) && !is.null(S1) && !is.null(S2)) {
    stop("provide either (M1 and M2) or (S1 and S2) but not both (M1,M2) and (S1,S2)")
  }
  if(!is.null(M1) && nrow(M1) > ncol(M1)) stop("The input appears to be S1 and S2, but the arguments were not specified; re-run with S1=<object> and S2=<object>")

  if(is.null(M1)) {
    nS = nrow(S1)
    if(nS!=nrow(S2)) stop('S1 and S2 must have the same number of rows')
    if(sum(apply(S1,2,tfun)) + sum(apply(S2,2,tfun))) stop('frobICA not defined when S1 or S2 has a column of all zeros')
    if(standardize) {
      S1 = scale(S1)
      S2 = scale(S2)
    }
    p = ncol(S1)
    q = ncol(S2)
    if(p < q) {
      S1 = cbind(S1,matrix(0,nS,(q-p)))
    }
    if(q < p) {
      S2 = cbind(S2,matrix(0,nS,(p-q)))
    }
    Stemp = matchICA(S=S1,template=S2)
    n.comp = max(q,p)
    indices = c(1:n.comp)[!(apply(Stemp,2,tfun) | apply(S2,2,tfun))]
    return(sqrt(sum((Stemp[,indices] - S2[,indices])^2))/sqrt(nS*min(p,q)))
  }

  else {
    if(sum(apply(M1,1,tfun)) + sum(apply(M2,1,tfun))) stop('frobICA not defined when M1 or M2 has a row of all zeros')
    if(standardize) {
      temp = diag((diag(M1%*%t(M1)))^(-1/2))
      M1 = temp%*%M1
      temp = diag((diag(M2%*%t(M2)))^(-1/2))
      M2 = temp%*%M2
    }
    p = ncol(M1)
    if(p!=ncol(M2)) stop("M1 and M2 must have the same number of columns")
    d = nrow(M1)
    q = nrow(M2)
    n.comp=max(d,q)
    if(n.comp > p) warning("M should be d x p")
    if(d<q) {
      M1 = rbind(M1,matrix(0,(q-d),p))
    }
    if(q<d) {
      M2 = rbind(M2,matrix(0,(d-q),p))
    }
    l2.mat1=l2.mat2=matrix(NA,nrow=n.comp,ncol=n.comp)
    for (j in 1:n.comp) {
      for (i in 1:n.comp) {
        #since signs are arbitrary, take min of plus and minus:
        l2.mat1[i,j]=sum((M2[i,]-M1[j,])^2)
        l2.mat2[i,j]=sum((M2[i,]+M1[j,])^2)
      }
    }
    l2.mat1=sqrt(l2.mat1)
    l2.mat2=sqrt(l2.mat2)
    #take the min of plus/min l2 distances. This is okay because solve_LSAP is one to one
    l2.mat=l2.mat1*(l2.mat1<=l2.mat2)+l2.mat2*(l2.mat2<l2.mat1)
    map=as.vector(solve_LSAP(l2.mat))
    #retain relevant l2 distances:
    l2.1=diag(l2.mat1[,map])
    l2.2=diag(l2.mat2[,map])
    #sign.change is for re-ordered matrix 2
    sign.change=-1*(l2.2<l2.1)+1*(l2.1<=l2.2)
    perm=diag(n.comp)[,map]%*%diag(sign.change)
    M.perm=t(perm)%*%M1
    indices = c(1:n.comp)[!(apply(M.perm,1,tfun) | apply(M2,1,tfun))]
    return(sqrt(sum((M.perm[indices,]-M2[indices,])^2))/sqrt(p*min(d,q)))
  }
}


#----------------
#----------------------------------------
#Function to make most extreme values for the skewed tail positive, i.e., force all distributions to be right skewed, and order ICs by skewness.
rightskew=function(S,M=NULL,order.skew=TRUE) {
  #S: n x d matrix
  #A: d x d` corresponding to X = S A
  #If order = TRUE, then ICs are organized from HIGHEST to LOWEST skewness where skewness is forced to be positive for all ICs.
  nObs <- nrow(S)
  if(ncol(S)>nObs) stop('S must be n x d')
  skewness<-function(x,n = nObs) (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  skew=apply(S,2,skewness)
  sign=-1*(skew<0)+1*(skew>0)
  S.new=S%*%diag(sign)
  if(!is.null(M)) M.new=t(diag(sign))%*%M
  if(order.skew==TRUE) {
    skew.new=apply(S.new,2,skewness)
    perm=order(skew.new,decreasing=TRUE)
    S.new=S.new[,perm]
    if(!is.null(M)) M.new=M.new[perm,]
  }
  if(is.null(M)) {
    S.new
  } else
    return(list(S=S.new,M=M.new))
}
#---------------------------------------------

#' multiple ProdenICA
#'
#' @param X Original data
#' @param n.comp the max components
#' @param restarts the number of restart
#' @param tol tol amount
#' @param maxit the max iteration number
#' @param G G method
#' @param verbose default = false
#' @param whiten whether to whiten
#' @param ... ellipsis
#'
#' @return
#' @export
#' @import ProDenICA
mProDenICA <- function(X, n.comp = ncol(X), restarts=0, tol=1e-07,maxit=100,G = c('GPois','G0','G1'),verbose=FALSE,whiten=FALSE,...) {
  ##NOTE: the restarts in ProDenICA evaluate the likelihood at a sample of orthogonal matrices, identifies the random matrix associated with highest likelihood, and then estimates ICs for this single initialization. Here, I initiate from the entire set of random matrices.
  ##NOTE: Restarts defined differently here than in ProDenICA. ProDenICA is initiatialized from restarts+1 initial values.
  ##NOTE: G defined differently from ProDenICA's Gfunc; here it is a string

  G = match.arg(G)
  if(G=='G0') Gfunc=G0
  if(G=='G1') Gfunc=G1
  if(G=='GPois') Gfunc=GPois
  est.list = list()
  runs = restarts+1
  obj.v = numeric(runs)
  theta.list = lapply(rep(choose(n.comp, 2), runs), runif, min = 0, max = 2 * pi)
  W.list = lapply(theta.list, theta2W)
  if(whiten) {
    a<- whitener(X=X,n.comp=n.comp)
    zData <- a$Z
  } else {
    zData <- X[,1:n.comp]
  }

  for(i in 1:runs) {
    est.list[[i]] = ProDenICA(x=zData, k=n.comp, W0=W.list[[i]], whiten=FALSE, maxit = maxit, thresh = tol, trace=verbose, restarts=0, Gfunc=Gfunc,...)
    if (G=='G1') obj.v[i] = calc.negent.hyvarinen(s=est.list[[i]]$s) else obj.v[i] = est.list[[i]]$negentropy
  }
  out = est.list[[which.max(obj.v)]]
  if(G=='G1') out$negentropy=obj.v[which.max(obj.v)]
  if(whiten) out$whitener=a$whitener
  out
}


givens.rotation <- function(theta=0, d=2, which=c(1,2))
{
  # For a given angle theta, returns a d x d Givens rotation matrix
  #
  # Ex: for i < j , d = 2:  (c -s)
  #                         (s  c)
  c = cos(theta)
  s = sin(theta)
  M = diag(d)
  a = which[1]
  b = which[2]
  M[a,a] =  c
  M[b,b] =  c
  M[a,b] = -s
  M[b,a] =  s
  M
}


# Returns square root of the precision matrix for whitening:
#' Returns square root of the precision matrix for whitening
#'
#' @param X Matrix
#' @param n.comp the number of components
#' @param center.row whether to center
#' @export
#' @import MASS
covwhitener <- function(X,n.comp=ncol(X),center.row=FALSE) {

  #X must be n x d
  if(ncol(X)>nrow(X)) warning('X is whitened with respect to columns')
  #Creates model of the form X.center = S A, where S are orthogonal with covariance = identity.
  x.center=scale(X,center=TRUE,scale=FALSE)
  if(center.row==TRUE) x.center = x.center - rowMeans(x.center)
  n.rep=dim(x.center)[1]
  covmat = cov(x.center)
  evdcov = eigen(covmat,symmetric = TRUE)
  whitener = evdcov$vectors%*%diag(1/sqrt(evdcov$values))%*%t(evdcov$vectors)
  #RETURNS PARAMETERIZATION AS IN fastICA (i.e., X is n x d)
  #NOTE: For whitened X, re-whitening leads to different X
  #The output for square A is equivalent to solve(K)
  return(list(whitener=whitener,Z=x.center%*%whitener,mean=apply(X,2,mean)))
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


SimLCA <- function(n.samples, distribution = c('logistic','t','gumbel'), nu.vector = c(sqrt(3)/pi,sqrt(3)/pi), snr, k = NULL, noisyICA=FALSE) {
  ##k = number of noise components
  require(MASS)
  if(is.null(k)) {
    if (noisyICA) k = 5 else k = 3
  }
  d = length(nu.vector)
  p = ifelse(noisyICA==TRUE,k,k+d)
  distribution=match.arg(distribution)
  sim.S = NULL
  if(distribution == 'logistic') for (i in 1:d) sim.S = cbind(sim.S,rlogis(n=n.samples,scale=nu.vector[i]))
  if(distribution == 't') for (i in 1:d) sim.S = cbind(sim.S,rt(n=n.samples,df=nu.vector[i]))
  if(distribution == 'gumbel') {
    require(evd)
    for (i in 1:d) sim.S = cbind(sim.S,rgumbel(n=n.samples,scale=nu.vector[i]))
  }
  sim.N = matrix(rnorm(n.samples*k,sd=1),nrow=n.samples)
  sim.M = myMixmat(p)
  sim.Ms = sim.M[1:d,]
  sim.Xs = sim.S%*%sim.Ms
  if(noisyICA) {
    sim.Mn = NULL
    sim.Xn = sim.N
  } else {
    sim.Mn = sim.M[(d+1):p,]
    sim.Xn = sim.N%*%sim.Mn
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
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn)) #eigenvalues sum to p
  sim.Xs = sqrt(snr)*sim.Xs #equivalent to scaled.sim.S%*%(alpha*sqrt(snr)*scalingMat%*%sim.Ms)+alpha*sqrt(snr)*matrix(mean.S,nrow=n.samples,ncol=2,byrow=TRUE)%*%sim.Ms
  #since we only recover scaled.sim.S, "true Ms" for, e.g., IFA, is defined as in scaled.sim.Ms
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)

  return(list(sim.S = sim.S, sim.N = sim.N, sim.Ms = sim.Ms, scaled.sim.Ms = scaled.sim.Ms, sim.Mn = sim.Mn, sim.X=sim.X, scaled.sim.S = scale(sim.S), scaled.sim.X = scale(sim.X), whitened.sim.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
}




#--------------------------------------------

# function to create network matrices from vectorized lower diagonals:

#' Create network matrices from vectorized lower diagonals
#' \code{vec2net} transfer the matrix vectorized lower diagonals into net to show the component image.
#' @param invector vectorized lower diagonals.
#' @param make.diag default value = 1.
#'
#' @return a net matrx
#' @export
#'
#' @examples net = vec2net(1:10)
vec2net = function(invector,make.diag=1) {
  #invector: choose(p,2) x 1, where p is the number of nodes
  nNode = (1 + sqrt(1+8*length(invector)))/2
  outNet = matrix(0,nNode,nNode)
  outNet[lower.tri(outNet)] = invector
  dim(outNet) = c(nNode,nNode)
  outNet = outNet + t(outNet)
  diag(outNet) = make.diag
  outNet
}


##############



whitener.evd = function(xData) {
  est.sigma = cov(xData)  ## Use eigenvalue decomposition rather than SVD.
  evd.sigma = svd(est.sigma)
  whitener = evd.sigma$u%*%diag(evd.sigma$d^(-1/2))%*%t(evd.sigma$u)
  list(whitener=whitener,zData = xData%*%whitener)
}

# here subjects are by rows, columns correspond to components
aveM = function(mjX,mjY) {

  mjX = t(t(mjX) / sqrt(apply(mjX^2,2,sum))) # each column has eucledean norm 1
  mjY = t(t(mjY) / sqrt(apply(mjY^2,2,sum)))
  n = nrow(mjX) # number of subjects
  rj = ncol(mjX) # number of components
  aveMj = matrix(0,n,rj)
  for (j in 1:rj) {
    # Need to take sign into account
    signXY = sign(sum(mjX[, j] * mjY[, j])) # sign
    temp = (mjX[,j] + mjY[,j] * signXY)/2
    aveMj[,j] = temp/sqrt(sum(temp^2))
  }
  aveMj
}




# BRisk 27 January 2020: added comments, did not change any calculations
computeRho <- function(JBvalX, JBvalY, Mxall, Myall, rjoint){
  # Compute the minimal value of rho so that the joint components have lower starting objective than individual
  # JBvalX, JBvalY - JB values for all components, first are assumed to be joint
  # Mxall, Myall - corresponding mixing matrices
  # rjoint - number of joint components
  rho = 0
  nx = length(JBvalX)
  ny = length(JBvalY)
  if ((nx == rjoint)|(ny == rjoint)) {return(0)}

  Mxscaled = Mxall %*% diag(sqrt(1/colSums(Mxall^2)))
  Myscaled = Myall %*% diag(sqrt(1/colSums(Myall^2)))

  Mxy = crossprod(Mxscaled, Myscaled)

  chordalJ = rep(0, rjoint)
  for (k in 1:rjoint){
    chordalJ[k] = 2 - 2 * Mxy[k, k]^2
  }

  for (i in (rjoint + 1):nx){ # ith individual from X
    for (j in (rjoint + 1):ny){ # jth individual from Y

      ratio = rep(NA, rjoint) # rhs/lhs for each joint
      #print(paste('pairs',i,j))
      for (k in 1:rjoint){ #kth joint
        rhs = max(JBvalX[i] - JBvalX[k] + JBvalY[j] - JBvalY[k], 0)
        # Brisk comment: 2 - 2*Mxy[i,j]^2 is chordal norm between individual components in X and Y
        lhs = 2 - 2 * Mxy[i, j]^2 - chordalJ[k]
        if (lhs <= 0){
          ratio[k] = NA
        }else{
          ratio[k] = rhs/lhs
        }
      }
      rho = max(rho, max(ratio, na.rm = T))
    }
  }
  return(rho)
}




#IGAY: adjust function so that the path is not hard-coded
#' PlotNetwork for components
#'
#' @param component component for the plot, px x rx
#' @param title title for plot
#' @param qmin default value = 0.005
#' @param qmax default value = 0.995
#' @param path default path
#' @param make.diag default value = NA
#'
#' @return a list as followed:
#' \describe{
#'       \item{\code{netmatfig}}{component loadings from rs correlation}
#'       \item{\code{loadingsfig}}{the sum of the absolute values of the rows of netmatfig}
#'       \item{\code{netmat}}{a matrix can be used with image() function}
#'       \item{\code{loadingsummary}}{loadings summary for each row of netmatfig}
#' }
#' @export
#' @import ggplot2
#' @import grid
#' @import scales
plotNetwork = function(component,title='',qmin=0.005, qmax=0.995, path = '~/Dropbox/JINGCA/Data/community_affiliation_mmpplus.csv',make.diag=NA) {
  # component:
  # vectorized network of length choose(n,2)


  # load communities for plotting:
  mmp_modules = read.csv(path)
  mmp_order = order(mmp_modules$Community_Vector)

  #check community labels:
  #table(mmp_modules$Community_Label)
  #table(mmp_modules$Community_Label,mmp_modules$Community_Vector)

  labels = c('VI','SM','DS','VS','DM','CE','SC')
  coords = c(0,70.5,124.5,148.5,197.5,293.5,360.5)


  zmin = quantile(component,qmin)
  zmax = quantile(component,qmax)

  netmat = vec2net(component,make.diag)

  meltsub = create.graph.long(netmat,mmp_order)
  #g2 = ggplot(meltsub, aes(X1,X2,fill=value))+ geom_tile()+ scale_fill_gradient2(low = "blue",  high = "red",limits=c(zmin,zmax),oob=squish)+labs(title = paste0("Component ",component), x = "Node 1", y = "Node 2")+coord_cartesian(clip='off',xlim=c(-0,390))

  g2 = ggplot(meltsub, aes(X1,X2,fill=value))+
    geom_tile()+
    scale_fill_gradient2(low = "blue",  high = "red",limits=c(zmin,zmax),oob=squish)+
    labs(title = title, x = "Node 1", y = "Node 2")+
    coord_cartesian(clip='off',xlim=c(-0,390))

  for (i in 1:7) {
    if (i!=3) {
      g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+10), xmax = (coords[i]+10), ymin = -7, ymax = -7)
    } else{
      g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+1), xmax = (coords[i]+1), ymin = -7, ymax = -7)
    }
  }
  # which nodes are prominent:
  loadingsummary = apply(abs(netmat),1,sum,na.rm=TRUE)
  loadingsum2 = loadingsummary[mmp_order]

  Community = factor(mmp_modules$Community_Label)[mmp_order]

  g3 = qplot(c(1:379),loadingsum2,col=Community,size=I(3))+xlab('MMP Index')+ylab('L1 Norm of the Rows')

  return(list(netmatfig = g2, loadingsfig = g3, netmat=netmat, loadingsummary = loadingsummary))
}
#########
#' Plot for the changed example of Y component.
#'
#' @param component component for the plot, px x rx
#' @param title title for plot
#' @param qmin default value = 0.005
#' @param qmax default value = 0.995
#' @param path default path
#' @param make.diag default value = NA
#'
#' @return a list as followed:
#' \describe{
#'       \item{\code{netmatfig}}{component loadings from rs correlation}
#'       \item{\code{loadingsfig}}{the sum of the absolute values of the rows of netmatfig}
#'       \item{\code{netmat}}{a matrix can be used with image() function}
#'       \item{\code{loadingsummary}}{loadings summary for each row of netmatfig}
#' }
#' @export
#' @import ggplot2
#' @import grid
#' @import scales
#'
plotNetwork_change = function(component,title='',qmin=0.005, qmax=0.995, path = '~/Dropbox/JINGCA/Data/community_affiliation_mmpplus.csv',make.diag=NA) {
  # component:
  # vectorized network of length choose(n,2)


  # load communities for plotting:
  mmp_modules = read.csv(path)
  mmp_order = order(mmp_modules$Community_Vector)

  #check community labels:
  #table(mmp_modules$Community_Label)
  #table(mmp_modules$Community_Label,mmp_modules$Community_Vector)

  #labels = c('VI','SM','DS','VS','DM','CE','SC')
  #coords = c(0,70.5,124.5,148.5,197.5,293.5,360.5)


  zmin = quantile(component,qmin)
  zmax = quantile(component,qmax)

  netmat = vec2net(component,make.diag)

  meltsub = create.graph.long(netmat,mmp_order)
  #g2 = ggplot(meltsub, aes(X1,X2,fill=value))+ geom_tile()+ scale_fill_gradient2(low = "blue",  high = "red",limits=c(zmin,zmax),oob=squish)+labs(title = paste0("Component ",component), x = "Node 1", y = "Node 2")+coord_cartesian(clip='off',xlim=c(-0,390))

  g2 = ggplot(meltsub, aes(X1,X2,fill=value))+
    geom_tile()+
    scale_fill_gradient2(low = "blue",  high = "red",limits=c(zmin,zmax),oob=squish)+
    labs(title = title, x = "Node 1", y = "Node 2")+
    coord_cartesian(clip='off',xlim=c(-0,100))

  #for (i in 1:7) {
  #  if (i!=3) {
  #    g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+10), xmax = (coords[i]+10), ymin = -7, ymax = -7)
  #  } else{
  #    g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+1), xmax = (coords[i]+1), ymin = -7, ymax = -7)
  #  }
  #}
  # which nodes are prominent:
  loadingsummary = apply(abs(netmat),1,sum,na.rm=TRUE)
  loadingsum2 = loadingsummary[mmp_order]

  Community = factor(mmp_modules$Community_Label)[mmp_order]

  g3 = qplot(c(1:100),loadingsum2,col=Community,size=I(3))+xlab('MMP Index')+ylab('L1 Norm of the Rows')

  return(list(netmatfig = g2, loadingsfig = g3, netmat=netmat, loadingsummary = loadingsummary))
}


############
#‘ Function for plotting networks with ggplot
#' create graph dataset with netmat and mmp_order
#' a data.frame called with vectorization of reordered netmat by mmp_order.
#' @param gmatrix netmat
#' @param sort_indices mmp_order
#'
#' @return a data.frame with vectors:
#'  ## X1: vector of numerics.
#'  ## X2: vector of numerics.
#'  ## value: vectorization of reordered netmat by mmp_order.
#'
#' @export
#'
create.graph.long = function(gmatrix,sort_indices=NULL) {
  nnode = nrow(gmatrix)
  X1 = c(1:nnode)%x%rep(1,nnode) # row number comes from Kronecker product
  X2 =  rep(1,nnode)%x%c(1:nnode) # column number
  if (!is.null(sort_indices)) {
    gmatrix = gmatrix[sort_indices,sort_indices]
  }
  value = as.vector(as.matrix(gmatrix))
  data.frame(X1,X2,value)
}



###############################################################
# Irina's auxillary functions for updates
###############################################################

###############

makecifti <- function(table, filename, template = '~/Dropbox/MyHCP/Data/100307/MNINonLinear/tfMRI_MOTOR_LR_Atlas.dtseries.nii', toolboxloc1 = '~/Dropbox/Applications2/hcp_ciftimatlab', toolboxloc2 = '~/Dropbox/mfunctions/robertoostenveld-cifti-matlab-27383b8',wbloc = '~/Applications2/workbench/bin_macosx64/wb_command',scratchloc = '~',matlabloc = '/Applications/MATLAB_R2019b.app/bin/matlab'){

  # Adapted from https://mandymejia.wordpress.com/2016/10/28/r-function-to-write-cifti-files/
  # table = Vxp matrix
  # table can be a vector if only a single CIFTI is needed
  # V = number of locations in both or one hemisphere
  # p = number of time points in .dtseries.nii
  # filename = vector length p, containing file names to be written (no extension)
  # template = name of existing CIFTI file that can be used as a template
  # toolboxloc = location and name of folder containing cifti-matlab or fieldtrip toolbox
  # scratchloc = location where CSV files will be written then deleted

  # this function writes a MATLAB script and executes it
  # must have MATLAB installed

  # Write table to CSV in scratch location
  fname <- file.path(scratchloc, 'tmp.csv')
  write.table(table, file=fname, row.names=FALSE, col.names=FALSE, sep=',')

  # Write MATLAB script
  line1 <- paste0("addpath '", toolboxloc1, "'")
  line2 <- paste0("addpath '", toolboxloc2, "'")
  line3 <- paste0("cifti = ciftiopen('",template,"','",wbloc,"');")
  line4 <- paste0("data = csvread('",scratchloc,"/tmp.csv')",";")
  line5 <- "[nrow,ncol]=size(data);"

  # Zero out entires, allowing for files without subcortical voxels:
  line6 <- "cifti.cdata = NaN(91282,ncol);"
  line7 <- "cifti.cdata(1:nrow,:)=data;"
  line8 <- paste0("ciftisave(cifti,'",filename,"','",wbloc,"');")
  matlab_lines <- c(line1, line2, line3, line4, line5, line6, line7, line8)
  writeLines(matlab_lines, con=file.path(scratchloc, 'myscript.m'))

  system(paste0(matlabloc," -nodisplay -r \"run('~/myscript.m'); exit\""))

  file.remove(fname)
  file.remove(file.path(scratchloc, 'myscript.m'))
}





##############
makeciftimmp <- function(v360, filename, template = '~/Dropbox/MyHCP/Data/100307/MNINonLinear/tfMRI_MOTOR_LR_Atlas.dtseries.nii', toolboxloc1 = '~/Dropbox/Applications2/hcp_ciftimatlab', toolboxloc2 = '~/Dropbox/mfunctions/robertoostenveld-cifti-matlab-27383b8',wbloc = '~/Applications2/workbench/bin_macosx64/wb_command',scratchloc = '~',matlabloc = '/Applications/MATLAB_R2019b.app/bin/matlab',mmpatlasloc="~/Dropbox/HCP_rsfmri_mmp/Data/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii"){

  # This function modifies makecifti to create a cifti file from 360 values corresponding to the glasser mmp atlas
  # makecifti is adapted from https://mandymejia.wordpress.com/2016/10/28/r-function-to-write-cifti-files/

  # Write table to CSV in scratch location
  fname <- file.path(scratchloc, 'tmp.csv')
  write.table(v360, file=fname, row.names=FALSE, col.names=FALSE, sep=',')

  # Write MATLAB script
  line1 <- paste0("addpath '", toolboxloc1, "'")
  line2 <- paste0("addpath '", toolboxloc2, "'")
  line3 <- paste0("cifti = ciftiopen('",template,"','",wbloc,"');")
  line4 = paste0("mmpatlas = ciftiopen('",mmpatlasloc,"','",wbloc,"');")
  line5 = paste0("ncii = size(cifti.cdata,1);")
  line6 = paste0("ncortex = 59412;")
  line7 <- paste0("data = csvread('",scratchloc,"/tmp.csv')",";")
  line8 <- "[~,ncol]=size(data);"

  # Zero out entires, allowing for files without subcortical voxels:
  line9 <- "cifti.cdata = NaN(91282,ncol);"
  # TO DO: Expand to subcortical
  line10 <- "mmpatlas_plus = [mmpatlas.cdata;zeros(ncii-ncortex,1)];"
  line11 <- "for r = 1:360"
  line12 <- "indices = mmpatlas_plus==r;"
  line13 <- "cifti.cdata(indices,:) = repmat(data(r,:),[sum(indices),1]);"
  line14 <- "end"
  line15 <- paste0("ciftisave(cifti,'",filename,"','",wbloc,"');")
  matlab_lines <- c(line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15)
  writeLines(matlab_lines, con=file.path(scratchloc, 'myscript.m'))

  system(paste0(matlabloc," -nodisplay -r \"run('~/myscript.m'); exit\""))

  file.remove(fname)
  file.remove(file.path(scratchloc, 'myscript.m'))
}
