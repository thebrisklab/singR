
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




#-------------------------------------
# Match mixing matrices:
# This function does not require M to be square:
#' Match mixing matrices
#'
#' @param M1 Subject score 1 matrix r x n.
#' @param M2 Subject score 2 matrix r x n.
#' @param S1 Loading 1 with matrix p x r.
#' @param S2 Loading 2 with matrix p x r.
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
  covmat = stats::cov(x.center)
  evdcov = eigen(covmat,symmetric = TRUE)
  whitener = evdcov$vectors%*%diag(1/sqrt(evdcov$values))%*%t(evdcov$vectors)
  #RETURNS PARAMETERIZATION AS IN fastICA (i.e., X is n x d)
  #NOTE: For whitened X, re-whitening leads to different X
  #The output for square A is equivalent to solve(K)
  return(list(whitener=whitener,Z=x.center%*%whitener,mean=apply(X,2,mean)))
}



#--------------------------------------------

##############



whitener.evd = function(xData) {
  est.sigma = stats::cov(xData)  ## Use eigenvalue decomposition rather than SVD.
  evd.sigma = svd(est.sigma)
  whitener = evd.sigma$u%*%diag(evd.sigma$d^(-1/2))%*%t(evd.sigma$u)
  list(whitener=whitener,zData = xData%*%whitener)
}

# here subjects are by rows, columns correspond to components
#' Average Mj for Mx and My
#'
#' @param mjX n x rj
#' @param mjY n x rj
#'
#' @return a new Mj
#' @export
#'
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


