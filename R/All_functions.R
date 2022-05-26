### All functions used in the package

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


#' Greedy Match
#'
#'
#'\code{Greedy Match} matches a column of Mx and My by minimizing chordal distance between vectors,
#'removes the matched columns and then finds the next pair.
#'This equivalent to maximizing absolute correlation for data in which each column has mean equal to zero.
#'Returns permuted columns of Mx and My. This function does not do any scaling or sign flipping.
#'For this matching to coincide with angle matching, the columns must have zero mean.
#' @param Mx Subject Score for X with n x n.comp matrix
#' @param My Subject Score for Y with n x n.comp matrix
#' @param Ux Matrix with n.comp x n, Mx = Lx^-1 \%*\% t Ux, Lx is the whitener matrix of dX.
#' @param Uy Matrix with n.comp x n, My = Ly^-1 \%*\% t Uy, Ly is the whitener matrix of dY.
#'
#' @return a list of matrices:
#' \describe{
#'        \item{\code{Mx}}{Columns of original Mx reordered from highest to lowest correlation with matched component in My}
#'        \item{\code{My}}{Columns of original My reordered from highest to lowest correlation with matched component in Mx}
#'        \item{\code{Ux}}{Permuted rows of original Ux corresponds to MapX}
#'        \item{\code{Uy}}{Permuted rows of original Uy corresponds to MapY}
#'        \item{\code{correlations}}{a vector of correlations for each pair of columns in permuted Mx and M}
#'        \item{\code{mapX}}{the sequence of the columns in original Mx.}
#'        \item{\code{mapY}}{the sequence of the columns in original MY.}
#' }
#'
#' @export
#'
greedymatch = function(Mx,My,Ux,Uy) {
  # input:
  # Mx: n x n.comp
  # My: n x n.comp
  # For this matching to coincide with angle matching, the columns must have zero mean.
  # NOTE: returns permuted columns of Mx and My. This function does not do any scaling or sign flipping.
  ## Mx: Columns of original Mx reordered from highest to lowest correlation with matched component in My

  # check the means are approximately zero:
  checkX = all(colMeans(Mx)<1e-05)
  checkY = all(colMeans(My)<1e-05)
  if(!checkX) warning('Columns of Mx do not appear to be centered -- correlations and angle match give different results')
  if(!checkY) warning('Columns of My do not appear to be centered -- correlations and angle match give different results')

  allCorr = abs(cor(Mx,My))
  px = ncol(Mx)
  py = ncol(My)
  allCorr = abs(cor(Mx,My))
  minpxpy = min(px,py)
  mapX = numeric(minpxpy)
  mapY = numeric(minpxpy)
  selcorr = numeric(minpxpy)
  for (t in 1:minpxpy) {
    selcorr[t] = max(allCorr)
    tindices=which(allCorr==selcorr[t],arr.ind=TRUE)
    mapX[t]=tindices[1]
    mapY[t]=tindices[2]
    allCorr[tindices[1],] = 0
    allCorr[,tindices[2]] = 0
  }

  # reorder Mx and My
  # grab all columns (applies when px ne py), non-matched columns
  # can still impact the selection of the penalty:
  notinX = c(1:px)[!c(1:px)%in%mapX]
  mapX = c(mapX,notinX)
  notinY = c(1:py)[!c(1:py)%in%mapY]
  mapY = c(mapY,notinY)
  reorderMx = Mx[,mapX]
  reorderMy = My[,mapY]
  reorderUx = Ux[mapX,]
  reorderUy = Uy[mapY,]
  return(list('Mx' = reorderMx, 'My' = reorderMy,'Ux' = reorderUx, 'Uy' = reorderUy, 'correlations' = selcorr, 'mapX' = mapX, 'mapY' = mapY))
}

# BRisk: 30 March 2020
#' Permutation test for joint component by Greedymatch
#'
#' @param MatchedMx matrix with nsubject x n.comp, comes from greedymatch
#' @param MatchedMy matrix with nsubject2 x n.comp, comes from greedymatch
#' @param nperm default value = 1000
#' @param alpha default value = 0.01
#' @param multicore default value = 0
#'
#' @return a list of matrixes
#'  ## rj: joint component rank
#'  ## pvalues: pvalue for the components(columns) not matched
#'  ## fwer_alpha: quantile of corr permutation with 1- alpha
#' @export
#'
permTestJointRank = function(MatchedMx,MatchedMy,nperm=1000,alpha=0.01,multicore=0) {
  # Mx: nSubject x n.comp
  # multicore: if = 0 or 1, uses a single processor; multicore > 1 requires packages doParallel and doRNG
  nSubject = nrow(MatchedMx)
  nSubject2 = nrow(MatchedMy)
  px = ncol(MatchedMx)
  py = ncol(MatchedMy)
  if (nSubject!=nSubject2) {
    stop('Mx and My have different numbers of rows. The number of subjects must be equal in X and Y.')
  }
  if (multicore>1) {
    require(doParallel)
    require(doRNG)
    registerDoParallel(multicore)

    corrperm = foreach (b = 1:nperm, .combine=rbind) %dorng% {
      new.My = MatchedMy[sample(1:nSubject),]
      max(abs(cor(MatchedMx,new.My)))
    }
  } else {
    corrperm = numeric(nperm)
    for (k in 1:nperm) {
      new.My = MatchedMy[sample(1:nSubject),]
      corrperm[k] = max(abs(cor(MatchedMx,new.My)))
    }
  }

  pperm = NULL
  maxrj = min(px,py)
  corrmatched = numeric(maxrj)
  for (i in 1:maxrj) {
    corrmatched[i] = abs(cor(MatchedMx[,i],MatchedMy[,i]))
  }
  for (j in 1:maxrj) pperm[j] = mean(corrmatched[j]<corrperm)
  rj = sum(pperm<alpha)

  return(list(rj=rj,pvalues = pperm,fwer_alpha = quantile(corrperm,1-alpha)))
}

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
#' @param alpha default = 0.8
#'
#' @return the sum of JB score for each component, the row of Sj.
#' @export
#'
calculateJB <- function(U, X, alpha = 0.8){
  # Calculate u^{top}X for each u_l and Xj, this is UX
  UX = U %*% X # Sx matrix with rx x p, Sx = Ux * Lx * Xc, in which X = Lx * Xc.
  p = ncol(X)
  # Calculate all individual components
  gamma = rowMeans(UX^3) # length r
  kappa = rowMeans(UX^4 - 3) # length r

  # TU must be r by n
  JB = sum(alpha*gamma^2 + (1-alpha)*kappa^2)

  return(JB)
}



# Curvilinear algorithm with r0 joint components
#' Curvilinear algorithm with r0 joint components
#'
#' @param Ux Matrix with n.comp x n, initial value of Ux, comes from greedyMatch.
#' @param Uy Matrix with n.comp x n, initial value of Uy, comes from greedyMatch.
#' @param xData matrix with n x px, Xw = Lx \%*\% Xc.
#' @param yData matrix with n x py, Yw = Ly \%*\% Yc.
#' @param invLx Inverse matrix of Lx, matrix n x n.
#' @param invLy Inverse matrix of Ly, matrix n x n.
#' @param rho the weight parameter of matching relative to non-gaussianity.
#' @param tau initial step size, default value is 0.01
#' @param alpha default value is 0.8
#' @param maxiter default value is 1000
#' @param tol the threshold of change in Ux and Uy to stop the curvlinear function
#' @param r0 the joint rank, comes from greedyMatch.
#'
#' @return a list of matrices:
#' \describe{
#'       \item{\code{Ux}}{Optimized Ux with matrix n.comp x n.}
#'       \item{\code{Uy}}{Optimized Uy with matrix n.comp x n.}
#'       \item{\code{tau}}{step size}
#'       \item{\code{iter}}{number of iterations.}
#'       \item{\code{error}}{PMSE(Ux,Uxnew)+PMSE(Uy,Uynew)}
#'       \item{\code{obj}}{Objective Function value}
#' }
#'
#'
#' @export
#'
curvilinear <- function(Ux, Uy, xData, yData, invLx, invLy, rho, tau = 0.01, alpha = 0.8, maxiter = 1000, tol = 1e-6, r0){

  tau1 = tau
  # Form standardized UX
  normLX = tcrossprod(invLx, Ux[1:r0, , drop = F])
  normLX = diag(sqrt(1/diag(crossprod(normLX))))%*% t(normLX)

  # Form standardized UY
  normLY = tcrossprod(invLy, Uy[1:r0, , drop = F])
  normLY = diag(sqrt(1/diag(crossprod(normLY))))%*% t(normLY)

  # Calculate objective value at current point
  obj1 = objectiveJoint(Ux, Uy, xData, yData, normLX, normLY, rho, alpha = alpha) #-JBpartX - JBpartY- 2*rho*Innerproduct
  print(obj1)

  error = 100
  iter = 0
  obj <- rep(0, maxiter)
  r = nrow(Ux)
  n = ncol(Ux)
  while ((error > tol)&(iter < maxiter)){
    iter <- iter + 1
    # Calculate gradient-type function at current point
    GUx <- calculateG(Ux, xData, invLx, A = normLY, rho, alpha = alpha, r0 = r0) # this is n by r
    GUy <- calculateG(Uy, yData, invLy, A = normLX, rho, alpha = alpha, r0 = r0) # this is n by r
    # Calculate skew-symmetric W
    Wx = GUx%*%Ux - crossprod(Ux, t(GUx))
    Wy = GUy%*%Uy - crossprod(Uy, t(GUy))

    # Generate new point Y
    Vtaux = tcrossprod(Ux, solve(diag(n) + tau*Wx/2, diag(n) - tau*Wx/2))
    Vtauy = tcrossprod(Uy, solve(diag(n) + tau*Wy/2, diag(n) - tau*Wy/2))
    # Form standardized UX
    normLXtau = tcrossprod(invLx, Vtaux[1:r0, ])
    normLXtau = diag(sqrt(1/diag(crossprod(normLXtau)))) %*% t(normLXtau)

    # Form standardized UY
    normLYtau = tcrossprod(invLy, Vtauy[1:r0,])
    normLYtau = diag(sqrt(1/diag(crossprod(normLYtau))))%*% t(normLYtau)

    # Check the value of objective at new Ytau
    obj2 = objectiveJoint(Vtaux, Vtauy, xData, yData, normLXtau, normLYtau, rho, alpha = alpha)
    # Make sure going along descent direction
    while ((obj2 > obj1)&(tau > 1e-14)){
      # Overshoot
      tau = 0.8*tau
      # New Vtau
      Vtaux = tcrossprod(Ux, solve(diag(n) + tau*Wx/2, diag(n) - tau*Wx/2))
      Vtauy = tcrossprod(Uy, solve(diag(n) + tau*Wy/2, diag(n) - tau*Wy/2))

      # Form standardized UX
      normLXtau = tcrossprod(invLx, Vtaux[1:r0,])
      normLXtau = diag(sqrt(1/diag(crossprod(normLXtau)))) %*% t(normLXtau)

      # Form standardized UY
      normLYtau = tcrossprod(invLy, Vtauy[1:r0,])
      normLYtau = diag(sqrt(1/diag(crossprod(normLYtau)))) %*% t(normLYtau)

      # New objective
      obj2 = objectiveJoint(Vtaux, Vtauy, xData, yData, normLXtau, normLYtau, rho, alpha = alpha)
    }
    if (obj2 < obj1){
      obj1 = obj2
      obj[iter] = obj2
      print(paste("Objective function value", obj2))
      Unewx = Vtaux
      Unewy = Vtauy
      normLY = normLYtau
      normLX = normLXtau
    }else{
      obj[iter] = obj1
      print(paste("Objective function value", obj1))
      Unewx = Ux
      Unewy = Uy
    }
    # How large is the change in objective
    error <- subspace_dist(Ux, Unewx) + subspace_dist(Uy, Unewy)
    Ux <- Unewx
    Uy <- Unewy
  }
  return(list(Ux = Unewx, Uy=Unewy, tau = tau, error = error, iter = iter, obj = obj[1:iter]))
}



### curvilinear_c function with r0 joint components
#' Curvilinear algorithm based on C code with r0 joint components
#'
#' @param Ux Matrix with n.comp x n, initial value of Ux, comes from greedyMatch.
#' @param Uy Matrix with n.comp x n, initial value of Uy, comes from greedyMatch.
#' @param xData matrix with n x px, Xw = Lx \%*\% Xc.
#' @param yData matrix with n x py, Yw = Ly \%*\% Yc.
#' @param invLx Inverse matrix of Lx, matrix n x n.
#' @param invLy Inverse matrix of Ly, matrix n x n.
#' @param rho the weight parameter of matching relative to non-gaussianity.
#' @param tau initial step size, default value is 0.01
#' @param alpha default value is 0.8
#' @param maxiter default value is 1000
#' @param tol the threshold of change in Ux and Uy to stop the curvlinear function
#' @param r0 the joint rank, comes from greedyMatch.
#' @return a list of matrices:
#' \describe{
#'       \item{\code{Ux}}{Optimized Ux with matrix n.comp x n.}
#'       \item{\code{Uy}}{Optimized Uy with matrix n.comp x n.}
#'       \item{\code{tau}}{step size}
#'       \item{\code{iter}}{number of iterations.}
#'       \item{\code{error}}{PMSE(Ux,Uxnew)+PMSE(Uy,Uynew)}
#'       \item{\code{obj}}{Objective Function value}
#' }
#' @export
#' @import Rcpp
#' @useDynLib singR
curvilinear_c <- function(Ux, Uy, xData, yData, invLx, invLy, rho, tau = 0.01, alpha = 0.8, maxiter = 1000, tol = 1e-6, r0) {

  return(updateUboth_c(Ux=Ux, Uy=Uy, xData=xData, yData=yData, invLx=invLx, invLy=invLy, rho=rho, r0=r0, alpha = alpha, tau = tau, maxiter = maxiter, tol = tol))


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


