library(singR)

data <- generateData_v3(nsubject = 48, snr = c(1, 1), vars = c(0.005, 0.005))
n = nrow(data$dX)
pX = ncol(data$dX)
pY = ncol(data$dY)
dXcentered <- data$dX - matrix(rowMeans(data$dX), n, pX, byrow = F)
dYcentered <- data$dY - matrix(rowMeans(data$dY), n, pY, byrow = F)

test_that("generateData_v3", {
  expect_equal(dim(data$dX), c(48,1089))
  expect_equal(dim(data$dY), c(48,4950))
})

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
