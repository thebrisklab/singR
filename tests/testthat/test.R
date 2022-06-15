library(singR)

data = data("exampledata")
# JB on X
estX_JB = lngca(xData = t(data$dX), n.comp = 12, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
Uxfull <- estX_JB$Ws  ## Ax = Ux %*% Lx, where Lx is the whitened matrix from covariance matrix of dX.
Mx_JB = est.M.ols(sData = estX_JB$S, xData = t(data$dX)) ## NOTE: for centered X, equivalent to xData %*% sData/(px-1)

# JB on Y
estY_JB = lngca(xData = t(data$dY), n.comp = 12, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
Uyfull <- estY_JB$Ws
My_JB = est.M.ols(sData = estY_JB$S, xData = t(data$dY))

test_that("linear non-Gaussian component analysis", {

  expect_equal(dim(estX_JB$Ws),c(48,12))

  expect_equal(dim(Mx_JB),c(12,48))

})


matchMxMy = greedymatch(t(Mx_JB), t(My_JB), Ux = t(Uxfull), Uy = t(Uyfull))
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
alpha = 0.8
tol = 1e-10

# curvilinear with small rho
rho = JBall/10
out_indiv_small <- curvilinear_c(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, rho = rho, tol = tol, alpha = alpha, maxiter = 1500, r0 = joint_rank)

test_that("Curvilinear search", {

  expect_equal(dim(out_indiv_small$Ux),c(12,48))
  expect_equal(dim(out_indiv_small$Uy),c(12,48))
})

