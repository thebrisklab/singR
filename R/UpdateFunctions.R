###############################################################
# Update functions
##############################################################
# Function that calculates joint objective value
objectiveJoint <- function(Ux, Uy, X, Y, normLX, normLY, rho, alpha = 0.8){
  JBpartX = calculateJB(Uy, Y, alpha = alpha)
  JBpartY = calculateJB(Ux, X, alpha = alpha)
  Innerproduct = sum(rowSums(normLX*normLY)^2)
  return(-JBpartX - JBpartY- 2*rho*Innerproduct)
}

# Function that calculates T(U), JB gradient with respect to U - DONE, but need to test more
calculateT <- function(U, X, alpha = 0.8, adjusted = F){
  # Calculate u^{top}X for each u_l and Xj, this is UX
  UX = U %*% X # r times p
  p = ncol(X)
  # Calculate all individual components
  gamma = rowMeans(UX^3) # length r
  kappa = rowMeans(UX^4 - 3) # length r
  prod1 = tcrossprod(UX^2, X)/p # r by n
  prod2 = tcrossprod(UX^3, X)/p # r by n
  # TU must be r by n
  TU = 6*alpha*diag(gamma)%*%prod1 + 8*(1-alpha)*diag(kappa)%*%prod2
  if (adjusted){
    # Use gradient adjustment as in Virta et al
    TU = TU - 24 * (1-alpha) * diag(kappa)%*%U
  }
  return(TU)
}

# Funciton that calculates full gradient with respect to current U - one function for X and Y
calculateG <- function(U, DataW, invL, A, rho, alpha = 0.8, r0 = nrow(U)){
  TU <- t(calculateT(U, X = DataW, alpha = alpha, adjusted = F)) # this makes it n by r
  GU <- -TU
  invLU <- tcrossprod(invL, U)
  normsLU2 <- colSums(invLU^2)
  for (i in 1:r0){
    uLa <- as.numeric(crossprod(invLU[ , i], A[i, ]))
    GU[ , i] <- GU[ , i] - 4 * rho * uLa * (invL %*% A[i, ]/normsLU2[i] - uLa * invL %*% invLU[ , i, drop = F]/(normsLU2[i]^2))
  }
  return(GU)
}


# Curvilinear algorithm with r0 joint components
#' Curvilinear algorithm with r0 joint components
#'
#' The curvilinear algorithm is modified from \href{https://www.semanticscholar.org/paper/A-feasible-method-for-optimization-with-constraints-Wen-Yin/d419879cdb80a87c3b7ab88e9f7478c1e70780ca}{Wen and Yin paper}.
#'
#' @param Ux Matrix with n.comp x n, initial value of Ux, comes from greedyMatch.
#' @param Uy Matrix with n.comp x n, initial value of Uy, comes from greedyMatch.
#' @param xData matrix with n x px, Xw = Lx \%*\% Xc.
#' @param yData matrix with n x py, Yw = Ly \%*\% Yc.
#' @param invLx Inverse matrix of Lx, matrix n x n.
#' @param invLy Inverse matrix of Ly, matrix n x n.
#' @param rho the weight parameter of matching relative to non-gaussianity.
#' @param tau initial step size, default value is 0.01
#' @param alpha controls weighting of skewness and kurtosis. Default value is 0.8, which corresponds to the Jarque-Bera test statistic with 0.8 weighting on squared skewness and 0.2 on squared kurtosis.
#' @param maxiter default value is 1000
#' @param tol the threshold of change in Ux and Uy to stop the curvlinear function
#' @param rj the joint rank, comes from greedyMatch.
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
curvilinear <- function(Ux, Uy, xData, yData, invLx, invLy, rho, tau = 0.01, alpha = 0.8, maxiter = 1000, tol = 1e-6, rj){
  r0=rj
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
#' #' The curvilinear algorithm is modified from \href{https://www.semanticscholar.org/paper/A-feasible-method-for-optimization-with-constraints-Wen-Yin/d419879cdb80a87c3b7ab88e9f7478c1e70780ca}{Wen and Yin paper}.
#'
#' @param Ux Matrix with n.comp x n, initial value of Ux, comes from greedyMatch.
#' @param Uy Matrix with n.comp x n, initial value of Uy, comes from greedyMatch.
#' @param xData matrix with n x px, Xw = Lx \%*\% Xc.
#' @param yData matrix with n x py, Yw = Ly \%*\% Yc.
#' @param invLx Inverse matrix of Lx, matrix n x n.
#' @param invLy Inverse matrix of Ly, matrix n x n.
#' @param rho the weight parameter of matching relative to non-gaussianity.
#' @param tau initial step size, default value is 0.01
#' @param alpha controls weighting of skewness and kurtosis. Default value is 0.8, which corresponds to the Jarque-Bera test statistic with 0.8 weighting on squared skewness and 0.2 on squared kurtosis.
#' @param maxiter default value is 1000
#' @param tol the threshold of change in Ux and Uy to stop the curvilinear function
#' @param rj the joint rank, comes from greedyMatch.
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
curvilinear_c <- function(Ux, Uy, xData, yData, invLx, invLy, rho, tau = 0.01, alpha = 0.8, maxiter = 1000, tol = 1e-6, rj) {

   return(updateUboth_c(Ux=Ux, Uy=Uy, xData=xData, yData=yData, invLx=invLx, invLy=invLy, rho=rho, r0=rj, alpha = alpha, tau = tau, maxiter = maxiter, tol = tol))


}



# Calculate subspace distance using two semi-orthogonal matrices of the same size
# Here U1, U2 are both r times n with orthonormal rows
subspace_dist <- function(U1, U2){
  # Use frobICA
  pmse(U1, U2)
}

