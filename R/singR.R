#' sing method for data integration
#'
#' @param dX original dataset for decomposition, matrix of n x px.
#' @param dY original dataset for decomposition, matrix of n x py.
#' @param n.comp the component number used in separate lngca process.
#' @param df default value=0 when use JB, df>0 when use tileted gaussian.
#' @param rho_extent small,medium or large.
#' @param Cplus whether to use C code in curvilinear search.
#' @param tol difference tolerance in curvilinear search.
#'
#' @return Function outputs a list including the following:
#' \describe{
#'       \item{\code{Sx}}{variable loadings for dataset X with matrix r x px.}
#'       \item{\code{Sy}}{variable loadings for dataset Y with matrix r x py.}
#'       \item{\code{Mxjoint}}{Mj of X data.}
#'       \item{\code{Myjoint}}{Mj of Y data.}
#'       \item{\code{est.Mj}}{estimated Mj for both datasets.}
#'       \item{\code{C_plus}}{whether to use C version of curvilinear search}
#'       \item{\code{rho_extent}}{the weight of rho in search}
#'       \item{\code{df}}{degree of freedom, = 0 when use JB,>0 when use tiltedgaussian}
#' }
#' @export
#'
singR <- function(dX,dY,n.comp=12,df=0,rho_extent=c('small','medium','large'),Cplus=T,tol = 1e-10) {
  # JB on X
  estX_JB = lngca(xData = t(dX), n.comp = n.comp, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB',stand = F,df=df) # what is the df at here.
  Uxfull <- estX_JB$U  ## Ax = Ux %*% Lx, where Lx is the whitened matrix from covariance matrix of dX.
  Mx_JB = est.M.ols(sData = t(estX_JB$S), xData = t(dX)) ## NOTE: for centered X, equivalent to xData %*% sData/(px-1)

  # JB on Y
  estY_JB = lngca(xData = t(dY), n.comp = n.comp, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB',stand = F,df=df)
  Uyfull <- estY_JB$U
  My_JB = est.M.ols(sData = t(estY_JB$S), xData = t(dY))

  matchMxMy = greedymatch(t(Mx_JB), t(My_JB), Ux = Uxfull, Uy = Uyfull)
  permJoint <- permTestJointRank(matchMxMy$Mx,matchMxMy$My) # alpha = 0.01, nperm=1000
  pval_joint = permJoint$pvalues
  joint_rank = permJoint$rj


  # Center X and Y
  n = nrow(dX)
  pX = ncol(dX)
  pY = ncol(dY)
  dXcentered <- standard(dX)
  dYcentered <- standard(dY)

  # For X
  # Scale rowwise
  est.sigmaXA = tcrossprod(dXcentered)/(pX-1)
  whitenerXA = est.sigmaXA%^%(-0.5)
  xDataA = whitenerXA %*% dXcentered
  invLx = est.sigmaXA%^%(0.5)
  Sx=matchMxMy$Ux[1:joint_rank] %*% xDataA

  # For Y
  # Scale rowwise
  est.sigmaYA = tcrossprod(dYcentered)/(pY-1)  ## since already centered, can just take tcrossprod
  whitenerYA = est.sigmaYA%^%(-0.5)
  yDataA = whitenerYA %*% dYcentered
  invLy = est.sigmaYA%^%(0.5)
  Sy=matchMxMy$Uy[1:joint_rank] %*% yDataA

  JBall = calculateJB(Sx) + calculateJB(Sy)
  rho_extent=match.arg(rho_extent)
  if(rho_extent=="small") {
      rho=JBall/10
  } else if(rho_extent=="medium") {
      rho=JBall
  } else {
      rho=JBall*10
    }

  if(Cplus){
    out_indiv <- curvilinear_c(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, rho = rho, tol = tol, maxiter = 1500, r0 = joint_rank)
  }else{
    out_indiv <- curvilinear(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, rho = rho, tol = tol, maxiter = 1500, r0 = joint_rank)
  }

  Sx = t(out_indiv$Ux[1:joint_rank, ] %*% xDataA)
  Sy = t(out_indiv$Uy[1:joint_rank, ] %*% yDataA)

  Mxjoint = tcrossprod(invLx, out_indiv$Ux[1:joint_rank, ])
  Myjoint = tcrossprod(invLy, out_indiv$Uy[1:joint_rank, ])

  est.Mj = aveM(Mxjoint,Myjoint)

  return(list(Sx=Sx,Sy=Sy,Mxjoint=Mxjoint,Myjoint=Myjoiny,est.Mj=est.Mj,Cplus=Cplus,rho_extent=rho_extent,df=df))

}
