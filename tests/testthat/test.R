library(singR)

data("exampledata")
data=exampledata



test_that("standardization",{
  dXcentered <- standard(data$dX)
  expect_equal(dim(dXcentered),c(48,1089))
})


#test singR
test_that("test for singR",{
  output <- singR(dX = data$dX,dY = data$dY,n.comp = 12,df = 0,rho_extent = "small",Cplus = TRUE,stand=FALSE)
  expect_equal(dim(output$Sjx),c(2,1089))
  expect_equal(dim(output$Sjy),c(2,4950))
  expect_equal(dim(output$Mxjoint),c(48,2))
})


# JB on X
estX_JB = lngca(xData = t(data$dX), n.comp = 12, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB',stand = F,df=0) # what is the df at here.
Uxfull <- estX_JB$U  ## Ax = Ux %*% Lx, where Lx is the whitened matrix from covariance matrix of dX.
Mx_JB = est.M.ols(sData = t(estX_JB$S), xData = t(data$dX)) ## NOTE: for centered X, equivalent to xData %*% sData/(px-1)

# JB on Y
estY_JB = lngca(xData = t(data$dY), n.comp = 12, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB',stand = F)
Uyfull <- estY_JB$U
My_JB = est.M.ols(sData = t(estY_JB$S), xData = t(data$dY))

test_that("linear non-Gaussian component analysis", {

  expect_equal(dim(estX_JB$U),c(12,48))
  expect_equal(dim(estX_JB$M),c(48,12))
  expect_equal(dim(Mx_JB),c(12,48))

})


matchMxMy = suppressWarnings(greedymatch(t(Mx_JB), t(My_JB), Ux = Uxfull, Uy = Uyfull))
permJoint <- permTestJointRank(matchMxMy$Mx,matchMxMy$My) # alpha = 0.01, nperm=1000
pval_joint = permJoint$pvalues
joint_rank = permJoint$rj

test_that("greedy match and permTest", {

  expect_equal(dim(matchMxMy$Mx),c(48,12))
  expect_equal(dim(matchMxMy$Ux),c(12,48))
  pval_joint = permJoint$pvalues
  expect_equal(length(pval_joint),12)
  expect_equal(permJoint$rj,2)

})


n = nrow(data$dX)
pX = ncol(data$dX)
pY = ncol(data$dY)
dXcentered <- data$dX - matrix(rowMeans(data$dX), n, pX, byrow = F)
dYcentered <- data$dY - matrix(rowMeans(data$dY), n, pY, byrow = F)



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

JBall = calculateJB(matchMxMy$Ux[1:2, ], X = xDataA) + calculateJB(matchMxMy$Uy[1:2, ], X = yDataA)



# JB and tolerance parameters
#alpha = 0.8
tol = 1e-10

# curvilinear with small rho
rho = JBall/10
out_indiv_small <- curvilinear_c(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, rho = rho, tol = tol, maxiter = 1500, r0 = joint_rank)

test_that("Curvilinear search", {

  expect_equal(dim(out_indiv_small$Ux),c(12,48))
  expect_equal(dim(out_indiv_small$Uy),c(12,48))
})






#' spmwm = 3*matrix(rnorm(100000),nrow=100)+1
#' dim(spmwm)
#' apply(spmwm,1,mean) # we want these to be 0
#' apply(spmwm,2,mean) # we want these to be 0
#' apply(spmwm,2,sd) # we want each of these variances to be 1
#'
#' spmwm_cp=standard(spmwm)
#' max(abs(apply(spmwm_cp,1,mean)))
#' max(abs(apply(spmwm_cp,2,mean)))
#' max(abs(apply(spmwm_cp,2,sd)-1))
stand=standard(data$dX,max.iter = 100)
mean(abs(apply(stand,2,sd)-1))

max(abs(apply(stand,2,mean)))

true_Sx=rbind(data$sjX,data$siX)
est_Sx=estX_JB$S
frobICA(S1 = t(true_Sx),S2 = t(est_Sx),standardize = T)

