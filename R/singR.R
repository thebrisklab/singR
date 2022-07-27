#' SImultaneous Non-Gaussian Component analysis for data integration.
#'
#' This function combines all steps from the \href{https://projecteuclid.org/journals/annals-of-applied-statistics/volume-15/issue-3/Simultaneous-non-Gaussian-component-analysis-SING-for-data-integration-in/10.1214/21-AOAS1466.full}{SING paper}
#' @param dX original dataset for decomposition, matrix of n x px.
#' @param dY original dataset for decomposition, matrix of n x py.
#' @param n.comp.X the number of non-Gaussian components in dataset X. If null, will estimate the number using ICtest::FOBIasymp.
#' @param n.comp.Y the number of non-Gaussian components in dataset Y. If null, will estimate the number using ICtest::FOBIasymp.
#' @param df default value=0 when use JB, if df>0, estimates a density for the loadings using a tilted Gaussian (non-parametric density estimate).
#' @param rho_extent Controls similarity of the scores in the two datasets. Numerical value and three options in character are acceptable. small, medium or large is defined from the JB statistic. Try "small" and see if the loadings are equal, then try others if needed. If numeric input, it will multiply the input by JBall to get the rho.
#' @param Cplus whether to use C code (faster) in curvilinear search.
#' @param tol difference tolerance in curvilinear search.
#' @param stand whether to use standardization, if true, it will make the column and row means to 0 and columns sd to 1. If false, it will only make the row means to 0.
#' @param distribution "JB" or "tiltedgaussian"; "JB" is much faster. In SING, this refers to the "density" formed from the vector of loadings. "tiltedgaussian" with large df can potentially model more complicated patterns.
#' @param maxiter the max iteration number for the curvilinear search.
#' @param individual whether to return the individual non-Gaussian components, default value = F.
#' @param whiten whitening method used in lngca. Defaults to "svd" which uses the n left eigenvectors divided by sqrt(px-1). Optionally uses the square root of the n x n "precision" matrix.
#' @return Function outputs a list including the following:
#' \describe{
#'       \item{\code{Sjx}}{variable loadings for joint NG components in dataset X with matrix rj x px.}
#'       \item{\code{Sjy}}{variable loadings for joint NG components in dataset Y with matrix rj x py.}
#'       \item{\code{Six}}{variable loadings for individual NG components in dataset X with matrix riX x px.}
#'       \item{\code{Siy}}{variable loadings for individual NG components in dataset Y with matrix riX x py.}
#'       \item{\code{Mix}}{scores of individual NG components in X with matrix n x riX.}
#'       \item{\code{Miy}}{scores of individual NG components in Y with matrix n x riY.}
#'       \item{\code{est.Mjx}}{Estimated subject scores for joint components in dataset X with matrix n x rj.}
#'       \item{\code{est.Mjy}}{Estimated subject scores for joint components in dataset Y with matrix n x rj.}
#'       \item{\code{est.Mj}}{Average of est.Mjx and est.Mjy as the subject scores for joint components in both datasets with matrix n x rj.}
#'       \item{\code{C_plus}}{whether to use C version of curvilinear search.}
#'       \item{\code{rho_extent}}{the weight of rho in search}
#'       \item{\code{df}}{degree of freedom, = 0 when use JB, >0 when use tiltedgaussian.}
#' }
#' @export
#' @examples
#' \dontrun{
#' #get simulation data
#' library(singR)
#' data(exampledata)
#'
#' # use JB stat to compute with singR
#' output_JB=singR(dX=exampledata$dX,dY=exampledata$dY,
#' df=0,rho_extent="small",distribution="JB",individual=TRUE)
#'
#' # use tiltedgaussian distribution to compute with singR.
#' # tiltedgaussian may be more accurate but is considerably slower,
#' # and is not recommended for large datasets.
#' output_tilted=singR(dX=exampledata$dX,dY=exampledata$dY,
#' df=5,rho_extent="small",distribution="tiltedgaussian",individual=TRUE)
#'
#' # use pmse to measure difference from the truth
#' pmse(M1 = t(output_JB$est.Mj),M2 = t(exampledata$mj),standardize = T)
#'
#' pmse(M1 = t(output_tilted$est.Mj),M2 = t(exampledata$mj),standardize = T)
#'
#' }

singR <- function(dX,dY,n.comp.X=NULL,n.comp.Y=NULL,df=0,rho_extent=c('small','medium','large'),Cplus=T,tol = 1e-10,stand=F,distribution="JB",maxiter=1500,individual=F,whiten = c('sqrtprec','eigenvec','none')) {


  #match.arg(c('small','medium','large'))
  #match.arg(rho_extent)
  whiten=match.arg(arg = NULL,choices = whiten)
  if ((nrow(dX)>ncol(dX)|nrow(dY)>ncol(dY)) & whiten != "sqrtprec") stop("with n > p, use whiten = 'sqrtprec'")


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

  if(is.null(n.comp.X)) {
    n.comp.X=NG_number(dXcentered)
  }
  if(is.null(n.comp.Y)) {
    n.comp.Y=NG_number(dYcentered)
  }

# JB on X
  estX_JB = lngca(xData = dXcentered, n.comp = n.comp.X, whiten = whiten, restarts.pbyd = 20, distribution=distribution,stand = F,df=df) # what is the df at here.
  Uxfull = estX_JB$U  ## Ax = Ux %*% Lx, where Lx is the whitened matrix from covariance matrix of dX.
  Mx_JB = est.M.ols(sData = estX_JB$S, xData = dXcentered) ## NOTE: for centered X, equivalent to xData %*% sData/(px-1)

  # JB on Y
  estY_JB = lngca(xData = dYcentered, n.comp = n.comp.Y, whiten = whiten, restarts.pbyd = 20, distribution=distribution,stand = F,df=df)
  Uyfull = estY_JB$U
  My_JB = est.M.ols(sData = estY_JB$S, xData = dYcentered)

  matchMxMy = greedymatch(scale(Mx_JB,scale=F), scale(My_JB,scale = F), Ux = Uxfull, Uy = Uyfull) # ignore the warnings of greedymatch, that the column is not scaled.
  permJoint <- permTestJointRank(matchMxMy$Mx,matchMxMy$My) # alpha = 0.01, nperm=1000
  # pval_joint = permJoint$pvalues
  joint_rank = permJoint$rj

  if(joint_rank==0) {
    stop("the joint_rank=0, to find individual non-Gaussian components can use lngca function separately.")
  }

  # For X
  # Scale rowwise
  est.sigmaXA = tcrossprod(dXcentered)/(pX-1)
  whitenerXA = est.sigmaXA%^%(-0.5)
  xDataA = whitenerXA %*% dXcentered
  invLx = est.sigmaXA%^%(0.5)



  # For Y
  # Scale rowwise
  est.sigmaYA = tcrossprod(dYcentered)/(pY-1)  ## since already centered, can just take tcrossprod
  whitenerYA = est.sigmaYA%^%(-0.5)
  yDataA = whitenerYA %*% dYcentered
  invLy = est.sigmaYA%^%(0.5)

  if(joint_rank==1){
    Sx=matrix(matchMxMy$Ux[1:joint_rank,],nrow=1) %*% xDataA
    Sy=matrix(matchMxMy$Uy[1:joint_rank,],nrow=1) %*% yDataA
  }
  if(joint_rank>1){
    Sx=matchMxMy$Ux[1:joint_rank,] %*% xDataA
    Sy=matchMxMy$Uy[1:joint_rank,] %*% yDataA
  }

  JBall = calculateJB(Sx) + calculateJB(Sy)

  rho_extent=rho_extent[1]

  if(is.character(rho_extent)) {

    if(rho_extent=="small") {
        rho=JBall/10
    } else if(rho_extent=="medium") {
        rho=JBall
    } else if(rho_extent=="large") {
        rho=JBall*10
      }
  } else if(is.numeric(rho_extent)){
    rho=rho_extent * JBall
  }


  if(Cplus){
    out_indiv <- curvilinear_c(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, rho = rho, tol = tol, maxiter = maxiter, rj = joint_rank)
  }else{
    out_indiv <- curvilinear(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, rho = rho, tol = tol, maxiter = maxiter, rj = joint_rank)
  }

 # calculate the joint components and the Mjoint
  if(joint_rank==1) {
    Sjx = matrix(out_indiv$Ux[1:joint_rank, ],nrow=1) %*% xDataA
    Sjy = matrix(out_indiv$Uy[1:joint_rank, ],nrow=1) %*% yDataA
    Mxjoint = tcrossprod(invLx, matrix(out_indiv$Ux[1:joint_rank, ],nrow=1))
    Myjoint = tcrossprod(invLy, matrix(out_indiv$Uy[1:joint_rank, ],nrow=1))
  }

  if(joint_rank>1) {
    Sjx = out_indiv$Ux[1:joint_rank, ] %*% xDataA
    Sjy = out_indiv$Uy[1:joint_rank, ] %*% yDataA
    Mxjoint = tcrossprod(invLx, out_indiv$Ux[1:joint_rank, ])
    Myjoint = tcrossprod(invLy, out_indiv$Uy[1:joint_rank, ])
  }

  sign_jx = signchange(Sjx,Mxjoint) # the input of signchange is r x px for S, and n x r for M, which is optional.
  Sjx = sign_jx$S
  Mxjoint = sign_jx$M

  sign_jy = signchange(Sjy,Mxjoint)
  Sjy = sign_jy$S
  Myjoint = sign_jy$M

  est.Mj = aveM(Mxjoint,Myjoint)

# get the individual NG components
  if(individual) {
    Six=NULL
    Siy=NULL
    Mx_I=NULL
    My_I=NULL

      if(n.comp.X>joint_rank){
        if(joint_rank+1==n.comp.X){
          Six = matrix(out_indiv$Ux[(joint_rank+1):n.comp.X, ],nrow = 1) %*% xDataA
          Mx_I = tcrossprod(invLx, matrix(out_indiv$Ux[(joint_rank+1):n.comp.X, ],nrow = 1))
        }else{
          Six = out_indiv$Ux[(joint_rank+1):n.comp.X, ] %*% xDataA
          Mx_I = tcrossprod(invLx, out_indiv$Ux[(joint_rank+1):n.comp.X, ])
        }
        sign_ix=signchange(Six,Mx_I)
        Six = sign_ix$S
        Mx_I = sign_ix$M

      }
      if(n.comp.Y>joint_rank){
        if(joint_rank+1==n.comp.Y){
          Siy = matrix(out_indiv$Uy[(joint_rank+1):n.comp.Y, ],nrow = 1) %*% yDataA
          My_I = tcrossprod(invLy, matrix(out_indiv$Uy[(joint_rank+1):n.comp.Y, ],nrow = 1))
        }else{
          Siy = out_indiv$Uy[(joint_rank+1):n.comp.Y, ] %*% yDataA
          My_I = tcrossprod(invLy, out_indiv$Uy[(joint_rank+1):n.comp.Y, ])
        }
        sign_iy=signchange(Siy,My_I)
        Siy = sign_iy$S
        My_I = sign_iy$M

      }
    return(list(Sjx=Sjx,Sjy=Sjy,Six=Six,Siy=Siy,Mix=Mx_I,Miy=My_I,est.Mj=est.Mj,est.Mjx=Mxjoint,est.Mjy=Myjoint,Cplus=Cplus,rho_extent=rho_extent,df=df))
  }
  return(list(Sjx=Sjx,Sjy=Sjy,est.Mj=est.Mj,est.Mjx=Mxjoint,est.Mjy=Myjoint,Cplus=Cplus,rho_extent=rho_extent,df=df))

}
