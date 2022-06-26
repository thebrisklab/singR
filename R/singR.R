#' sing method for data integration
#'
#' @param dX original dataset for decomposition, matrix of n x px.
#' @param dY original dataset for decomposition, matrix of n x py.
#' @param n.comp the component number used in separate lngca process.
#' @param df default value=0 when use JB, df>0 when use tileted gaussian.
#' @param rho_extent small,medium or large.
#' @param Cplus whether to use C code in curvilinear search.
#' @param tol difference tolerance in curvilinear search.
#' @param stand whether to use standardization, if it was true, it will make the column and row means to 0 and columns sd to 1. If false, it will only make the row means to 0.
#' @param distribution "JB" or "tiltedgaussian"
#' @return Function outputs a list including the following:
#' \describe{
#'       \item{\code{Sjx}}{variable loadings for dataset X with matrix rj x px.}
#'       \item{\code{Sjy}}{variable loadings for dataset Y with matrix rj x py.}
#'       \item{\code{Mxjoint}}{Mj of X data with matrix n x rj.}
#'       \item{\code{Myjoint}}{Mj of Y data with matrix n x rj.}
#'       \item{\code{est.Mj}}{estimated Mj for both datasets with matrix n x rj.}
#'       \item{\code{C_plus}}{whether to use C version of curvilinear search.}
#'       \item{\code{rho_extent}}{the weight of rho in search}
#'       \item{\code{df}}{degree of freedom, = 0 when use JB,>0 when use tiltedgaussian.}
#' }
#' @export
#' @examples
#' \dontrun{
#' #get simulation data
#' data(exampledata)
#' data=exampledata
#'
#' # use JB stat to compute with singR
#' output_JB=singR(dX=data$dX,dY=data$dY,n.comp=12,df=0,rho_extent="small",distribution="JB")
#' output_tilted=singR(dX=data$dX,dY=data$dY,n.comp=12,df=3,rho_extent="small",distribution="tiltedgaussian")
#'
#' # use frobICA to show the difference of JB and tiltedgaussian
#' frobICA_JB=frobICA(M1 = t(output_JB$est.Mj),M2 = t(data$mj),standardize = T) #0.0071682
#' frobICA_tilted=frobICA(M1 = t(output_tilted$est.Mj),M2 = t(data$mj),standardize = T) #0.0071295
#'
#' }

singR <- function(dX,dY,n.comp=12,df=0,rho_extent=c('small','medium','large'),Cplus=T,tol = 1e-10,stand=F,distribution="JB") {
  # Center X and Y
  if (stand) {
    n = nrow(dX)
    pX = ncol(dX)
    pY = ncol(dY)
    dXcentered <- standard(dX)
    dYcentered <- standard(dY)
  }else{
    n = nrow(dX)
    pX = ncol(dX)
    pY = ncol(dY)
    dXcentered <- dX - matrix(rowMeans(dX), n, pX, byrow = F)
    dYcentered <- dY - matrix(rowMeans(dY), n, pY, byrow = F)
  }


# JB on X
  estX_JB = lngca(xData = t(dXcentered), n.comp = n.comp, whiten = 'sqrtprec', restarts.pbyd = 20, distribution=distribution,stand = F,df=df) # what is the df at here.
  Uxfull <- estX_JB$U  ## Ax = Ux %*% Lx, where Lx is the whitened matrix from covariance matrix of dX.
  Mx_JB = est.M.ols(sData = t(estX_JB$S), xData = t(dXcentered)) ## NOTE: for centered X, equivalent to xData %*% sData/(px-1)

  # JB on Y
  estY_JB = lngca(xData = t(dYcentered), n.comp = n.comp, whiten = 'sqrtprec', restarts.pbyd = 20, distribution=distribution,stand = F,df=df)
  Uyfull <- estY_JB$U
  My_JB = est.M.ols(sData = t(estY_JB$S), xData = t(dYcentered))

  matchMxMy = suppressWarnings(greedymatch(t(Mx_JB), t(My_JB), Ux = Uxfull, Uy = Uyfull)) # ignore the warnings of greedymatch, that the column is not scaled.
  permJoint <- permTestJointRank(matchMxMy$Mx,matchMxMy$My) # alpha = 0.01, nperm=1000
  pval_joint = permJoint$pvalues
  joint_rank = permJoint$rj



  # For X
  # Scale rowwise
  est.sigmaXA = tcrossprod(dXcentered)/(pX-1)
  whitenerXA = est.sigmaXA%^%(-0.5)
  xDataA = whitenerXA %*% dXcentered
  invLx = est.sigmaXA%^%(0.5)
  Sx=matchMxMy$Ux[1:joint_rank,] %*% xDataA

  # For Y
  # Scale rowwise
  est.sigmaYA = tcrossprod(dYcentered)/(pY-1)  ## since already centered, can just take tcrossprod
  whitenerYA = est.sigmaYA%^%(-0.5)
  yDataA = whitenerYA %*% dYcentered
  invLy = est.sigmaYA%^%(0.5)
  Sy=matchMxMy$Uy[1:joint_rank,] %*% yDataA

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

  Sjx = out_indiv$Ux[1:joint_rank, ] %*% xDataA
  Sjy = out_indiv$Uy[1:joint_rank, ] %*% yDataA

  Mxjoint = tcrossprod(invLx, out_indiv$Ux[1:joint_rank, ])
  Myjoint = tcrossprod(invLy, out_indiv$Uy[1:joint_rank, ])

  est.Mj = aveM(Mxjoint,Myjoint)

  return(list(Sjx=Sjx,Sjy=Sjy,Mxjoint=Mxjoint,Myjoint=Myjoint,est.Mj=est.Mj,Cplus=Cplus,rho_extent=rho_extent,df=df))

}
