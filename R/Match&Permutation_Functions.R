### Functions for Match and permutation


# Use greedy algorithm:

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
#' Title
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
    #require(doParallel)
    #require(doRNG)
    #registerDoParallel(multicore)

    #corrperm = foreach (b = 1:nperm, .combine=rbind) %dorng% {
      #new.My = MatchedMy[sample(1:nSubject),]
      #max(abs(cor(MatchedMx,new.My)))
    stop('This package does not support multicore')
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


#-------------------------
#Match based on L2 distances and Hungarian
#' match ICA
#'
#' @param S loading variable matrix
#' @param template template for match
#' @param M subject score matrix
#'
#' @return
#'
#' @import clue
matchICA=function(S,template,M=NULL) {

  n.comp=ncol(S)
  n.obs=nrow(S)
  if(n.comp>n.obs) warning('Input should be n x d')
  if(!all(dim(template)==dim(S))) warning('Template should be n x d')
  S = t(S)
  template = t(template)
  l2.mat1=matrix(NA,nrow=n.comp,ncol=n.comp)
  l2.mat2=l2.mat1
  for (j in 1:n.comp) {
    for (i in 1:n.comp) {
      l2.mat1[i,j]=sum((template[i,]-S[j,])^2)/n.obs
      l2.mat2[i,j]=sum((template[i,]+S[j,])^2)/n.obs
    }
  }
  l2.mat1=sqrt(l2.mat1)
  l2.mat2=sqrt(l2.mat2)
  l2.mat=l2.mat1*(l2.mat1<=l2.mat2)+l2.mat2*(l2.mat2<l2.mat1)
  map=as.vector(solve_LSAP(l2.mat))
  l2.1=diag(l2.mat1[,map])
  l2.2=diag(l2.mat2[,map])
  sign.change=-1*(l2.2<l2.1)+1*(l2.1<=l2.2)
  perm=diag(n.comp)[,map]%*%diag(sign.change)

  s.perm=t(perm)%*%S
  if(!is.null(M)) {
    M.perm=t(M)%*%perm
    return(list(S=t(s.perm),M=t(M.perm)))
  }  else {
    t(s.perm)
  }
}

# IGAY: edited on June 4th, 2019 to return perm matrix and omangles ordering
# IGAY: edited on July 15th, 2019 to avoid exiting with error when rx>ry
# IGAY: edited on Aug 12th, 2019 to avoid dropping to vector type from matrix when original rank is 1
#' Match the colums of Mx and My
#'
#'\code{angleMatchICA} match the colums of Mx and My, using the n x p parameterization of the JIN decomposition assumes
#' @param Mx Subject score for X   matrix of n x n.comp
#' @param My Subject score for Y   matrix of n x n.comp
#' @param Sx Variable loadings for X  matrix of n.comp x px
#' @param Sy Variable loadings for Y  matrix of n.comp x py
#'
#' @return a list of matrixes:
#'  ## Mx:
#'  ## My:
#'  ## matchedangles:
#'  ## allangles:
#'  ## perm:
#'  ## omangles:
#' @export
#' @import clue
angleMatchICA=function(Mx,My,Sx=NULL,Sy=NULL) {
  # match the colums of Mx and My, using the
  # n x p parameterization of the JIN decomposition
  # assumes

  rx = ncol(Mx)
  ry = ncol(My)
  #if(rx>ry) stop('rx must be less than or equal to ry')
  n.comp = max(ry, rx) # IGAY: adjusted this on Jyly 19th, 2019 to be maximal of the two

  n.obs=nrow(Mx)
  if(n.comp>=n.obs) warning('Input should be n x r')
  if(n.obs!=nrow(My)) warning('Mx and My need to have the same number of rows')
  Mx = t(Mx)
  Mx = Mx / sqrt(apply(Mx^2,1,sum))


  My = t(My)
  My = My / sqrt(apply(My^2,1,sum))

  if(rx<ry) {
    Mx = rbind(Mx,matrix(0,ry-rx,n.obs))
  }
  if(ry<rx) {
    My = rbind(My,matrix(0,rx-ry,n.obs))
  }
  angle.mat1=acos(Mx%*%t(My))
  angle.mat2=acos(-1*Mx%*%t(My))
  angle.mat=angle.mat1*(angle.mat1<=angle.mat2)+angle.mat2*(angle.mat2<angle.mat1)
  map=as.vector(solve_LSAP(angle.mat))
  angle.mat1.perm = angle.mat1[,map]
  angle.mat2.perm = angle.mat2[,map]
  angle1=diag(angle.mat1.perm)
  angle2=diag(angle.mat2.perm)
  matchedangles = apply(cbind(angle1,angle2),1,min)
  allangles = angle.mat1.perm*(angle.mat1.perm<=angle.mat2.perm)+angle.mat2.perm*(angle.mat2.perm<angle.mat1.perm)
  sign.change=-1*(angle2<angle1)+1*(angle1<=angle2)
  perm=diag(n.comp)[,map]%*%diag(sign.change)

  My.perm=t(perm)%*%My

  # reorder components by their matched angles
  smatchedangles = sort(matchedangles)
  omangles = order(matchedangles)
  sallangles = allangles[omangles[1:rx],omangles[1:ry]]
  sMx = Mx[omangles[1:rx],, drop = F]
  sMy.perm = My.perm[omangles[1:ry],, drop = F]

  if(!is.null(Sy)) {
    Sy.perm=t(perm)%*%Sy
    sSy.perm = Sy.perm[omangles[1:ry],, drop = F]
    sSx = Sx[omangles[1:rx],, drop = F]
    return(list(Mx=t(sMx),My = t(sMy.perm),matchedangles = smatchedangles,allangles = sallangles,Sx = sSx, Sy = sSy.perm, perm = perm, omangles = omangles))
  }
  else {
    return(list(Mx=t(sMx),My = t(sMy.perm),matchedangles = smatchedangles,allangles = sallangles, perm = perm, omangles = omangles))
  }
}




# Wrapper functions for rank estimation using permutation approach
#######################################################################

permmatRank_sequential_JB = function(xdata,maxcompdata=ncol(xdata),ncomplist,nperms,ninitialization_data=10,ninitialization_perm=5) {
  #xdata: p x n subjects
  #maxcompdata: number of components to estimate from data. This should be the maximum number of possible non-Gaussian components. By default equal to n
  #ncomplist is a vector corresponding to the components to be tested for Gaussianity.
  # In the simplest case, this, will be 1:maxcompdata
  #for each ncomp in ncomplist, ncomp refers to a test of whether the ncomp$th$ component is Gaussian
  #the maximum of this should be less than or equal to maxcompdata
  #nperms: number of samples for each permutation test
  #ncores: number of cores to use via registerDoParallel

  # subfunctions used in this function:
  permmatRank = function(xdata, ncomp, nperms,
                         ninitialization_perm) {
    nsubjects = ncol(xdata)
    px = nrow(xdata)
    ngauss = nsubjects - ncomp + 1
    permmatJB = rep(0, nperms)
    for(k in 1:nperms){
      # sample n - r subjects
      tempX = xdata[ , sample(1:nsubjects, ngauss)]
      newX = matrix(0, px, ngauss)
      for (j in 1:ngauss) {
        # permute jth column
        pmat = sample(1:px)
        newX[,j] = tempX[pmat, j]
      }
      # estimate non-gaussianity of the component
      permmatJB[k] = mlcaFP(newX, n.comp = 1, restarts.pbyd=ninitialization_perm,distribution='JB')$nongaussianity
    }
    permmatJB
  }


  # Estimate model with components=maxcompdata :
  estX = mlcaFP(xdata,restarts.pbyd=round(ninitialization_data/2),restarts.dbyd=round(ninitialization_data/2),distribution='JB', n.comp=maxcompdata)

  # Construct p-values for reach component
  nc = length(ncomplist)
  permmatJB_bigmat = matrix(0, nrow = nperms, ncol = nc)
  pvalues = rep(NA, nc)

  for (i in 1:length(ncomplist)) {
    # Estimate residuals from the first r-1 components
    if (ncomplist[i] > 1){
      newxdata = lm(xdata~estX$S[,1:(ncomplist[i]-1)])$residuals
    }else{
      newxdata = xdata
    }
    # Find values of non-gaussianity from all permutations
    permmatJB_bigmat[,i] = permmatRank(newxdata, ncomplist[i], nperms, ninitialization_perm)
    # Calculate corresponding p-value
    pvalues[i] = mean(estX$nongaussianity[ncomplist[i]]<permmatJB_bigmat[,i])
  }
  colnames(permmatJB_bigmat) = ncomplist
  return(list(pvalues = pvalues,sample=permmatJB_bigmat))
}


#' Permutation test to get joint components ranks
#'
#' @param matchedResults results generated by angleMatchICA
#' @param nperms the number of permutation
#'
#' @return a list of matrixes
#'  ## pvalues: pvalues for the matched colunmns don't have correlation.
#'  ## corrperm: correlation value for original Mx with each random permutation of My.
#'  ## corrmatched: the correlation for each pair of matched columns.
#' @export
#'
permmatRank_joint = function(matchedResults, nperms = 100){

  # Calcualte correlations based on original Ms via angle-matching
  Mx = t(matchedResults$Mx)
  My = t(matchedResults$My)
  corrmatched = cos(matchedResults$matchedangles)
  n = ncol(Mx)

  # Calculate maximal correlations based on permutations
  corrperm = rep(NA, nperms)
  for(b in 1:nperms){
    # Resample the subjects in the 2nd mixing matrix
    new.My= My[,sample(1:n), drop = F]
    # Calculate maximal correlation
    corrperm[b] = max(abs(atanh(cor(t(Mx),t(new.My)))))
  }

  # Calculate p-values from permutations
  rx = nrow(Mx)
  ry = nrow(My)
  maxrj = min(rx,ry) #maximum possible number of joint components

  pperm = rep(NA, maxrj)
  for (j in 1:maxrj) {
    pperm[j] = mean(atanh(corrmatched[j]) < corrperm)
  }

  return(list(pvalues = pperm, corrperm = corrperm, corrmatched = corrmatched))
}



###############
