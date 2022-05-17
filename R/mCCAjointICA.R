

# Back-project on given M to find S subject to orthogonal constraints (this is Orthogonal Procrustes Problem)
# Sold - r times p previous matrix of components (for joint structure only)
# Mold - n times r previous mixing matrix (for joint structure only)
# Mave - n times r new estimated mixing matrix (based on averaging, norm 1)
# eps - itearations for Dx adjustment
est.S.backproject <- function(Sold, Mold, Mave, eps = 1e-1){
  r <- ncol(Mold)
  p <- ncol(Sold)
  J = Mold %*% Sold # previous joint structure, n by p
  error = 10
  Dold = diag(Sold %*% crossprod(J, Mave))/sqrt(p - 1)
  while(error > eps){
    # Re-estimate new S (orthogonal Procrustes problem)
    svdJM <- svd(crossprod(J, Mave %*% diag(Dold)))
    Snew <- tcrossprod(svdJM$v[ , 1:r], svdJM$u[, 1:r])
    # Re-estimate new scaling
    Dnew = diag(Snew%*% crossprod(J, Mave))
    # Evaluate error
    error = max(abs(Dold - Dnew))
    # Update Dold
    Dold = Dnew
  }
  return(list(S = Snew * sqrt(p - 1), D = Dnew))
}

# Just joint ICA
#' Find the joint component through joint_ICA method.
#'
#' @param dX the first dataset.
#' @param dY the second dataset.
#' @param r0 the number of joint components.
#'
#' @return a list with matrices
#' ## S: the concatenated matrix of SjX and SjY, [SjX',SjY'],at here (pX+pY) x rj
#' ## Mjoint: subject score estimated by joint ICA, at here rj x nsubject. for the example dataset, the rj is 2.
#'
#'
#' @export
jointICA <- function(dX, dY, r0 = 2){
  # Normalization step, divide by \|X\|_F^2 each part
  dXsA <- dX/sqrt(mean(dX^2))
  dYsA <- dY/sqrt(mean(dY^2))

  # Concatenate them together [X, Y] and perform SVD (PCA)
  dXY <- cbind(dXsA, dYsA) # [X, Y] ~ UDV'
  svdXY <- svd(dXY, nu = r0, nv = r0) # 2 is the true number of joint components

  # Scale V to have sd 1 (this is essentially whitening step)
  #sds <- apply(svdXY$v, 2, sd)
  #Vscaled <- svdXY$v %*% diag(1/sds)
  Vscaled <- svdXY$v

  # JB on V, PCA loadings from above
  estXY_JB = lngca(xData = Vscaled, n.comp = 2, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
  # Get the corresponding M by Least Squares on original data
  Mjoint = est.M.ols(sData = estXY_JB$S, xData = t(dXY))

  # What are the errors?
  # Return the output
  return(list(S = estXY_JB$S, Mjoint = Mjoint))
}

# Implementation for mCCA + joint ICA
# Mx - # of components in initial PCA of X
# My - # of components in initial PCA of Y
# M - # of joints components, has to be <= min(Mx, My)
mCCA_jointICA <- function(dX, dY, Mx, My, M = min(Mx, My)){
  # Whitening/centering of X
  n = nrow(dX)
  pX = ncol(dX)
  pY = ncol(dY)


  # do column centering
  dXcentered <- dX - matrix(colMeans(dX), n, pX, byrow = T)
  dYcentered <- dY - matrix(colMeans(dY), n, pY, byrow = T)

  # In Anylisis pipeline in Siu et al., 2011, on page 844, apply squared normalization step
  dXcentered <- dXcentered / sqrt(mean(dXcentered^2))
  dYcentered <- dYcentered / sqrt(mean(dYcentered^2))

  # Perform separate SVD on each only keeping the first Mk right singular vectors
  svdX = svd(dXcentered, nu = Mx, nv = Mx)
  svdY = svd(dYcentered, nu = My, nv = My)

  Ex = svdX$v
  Ey = svdY$v

  # Form the Y in the paper
  YdX = dXcentered %*% Ex # This is essentially UxDx
  YdY = dYcentered %*% Ey # This is essentially UyDy

  # Perform CCA to get Dx and Dy. Here I assume YdX and YdY are column-centered. If they are not, then the crossprod in the middle needs to be adjusted with centering. The scaling should not be adjusted to match the paper's scaling
  cca_svd = svd(diag(1/svdX$d[1:Mx]) %*% crossprod(YdX, YdY) %*% diag(1/svdY$d[1:My]), nu = min(Mx, My), nv = min(Mx, My))
  Dx = YdX %*% diag(1/svdX$d[1:Mx]) %*% cca_svd$u[ , 1:M] # n by M
  Dy = YdY %*% diag(1/svdY$d[1:My]) %*% cca_svd$v[ , 1:M] # n by M

  # Form Ck from the paper - these are not mean 0, nor are they orthogonal
  Cx = crossprod(Dx, dXcentered) # M by pX; X approx DxCx
  Cy = crossprod(Dy, dYcentered) # M by pY; Y approx DyCy

  # JB on CxS, CyS from above
  estCxy_JB = mlcaFP(xData = t(cbind(Cx, Cy)), n.comp = 2, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')

  # Get the corresponding M by Least Squares on Cx, Cy (this is W^{-1} in the paper's notation)
  MjointC = est.M.ols(sData = estCxy_JB$S, xData = t(cbind(Cx, Cy)))
  # Cx = M Sx; Cy = M Sy

  # Get Ax and Ay here, correlated almost the same way as Dx and Dy but not quite orthogonal (approximately so), seems consistent with paper
  Ax = Dx %*% MjointC # X approx Ax Sx = Dx Cx = Dx M Sx
  Ay = Dy %*% MjointC # Y approx Ay Sy = Dy Cy = Dy M Sy

  # Need to adjust back for scaling, otherwise this is wrong M
  Mx = est.M.ols(sData = estCxy_JB$S[1:pX, ], xData = t(dX))
  My = est.M.ols(sData = estCxy_JB$S[(pX+1):(pX+pY), ], xData = t(dY))
  Mave = t(aveM(t(Mx), t(My))) # average two mixing matrices

  # Return the output
  return(list(Ax = Ax, Ay = Ay, S = estCxy_JB$S, Mx = Mx, My = My, Mave = Mave))
}
