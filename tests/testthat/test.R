library(singR)

data("exampledata")
data=exampledata


# test standard
test_that("standardization",{
  dXcentered <- standard(data$dX)
  expect_equal(dim(dXcentered),c(48,1089))
})

# test NG_number
test_that("NG_number",{
  expect_equal(NG_number(data$dX),4)
})



#test singR
test_that("test for singR",{
  output <- singR(dX = data$dX,dY = data$dY,df = 0,rho_extent = "small",Cplus = TRUE,stand=FALSE,n.comp.X = 4,n.comp.Y = 4,individual = T)
  expect_equal(dim(output$Sjx),c(2,1089))
  expect_equal(dim(output$Sjy),c(2,4950))
  expect_equal(dim(output$est.Mj),c(48,2))
})

n.comp=4

# JB on X
estX_JB = lngca(xData = data$dX, n.comp = n.comp, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB',stand = F,df=0) # what is the df at here.
Uxfull <- estX_JB$U  ## Ax = Ux %*% Lx, where Lx is the whitened matrix from covariance matrix of dX.
Mx_JB = est.M.ols(sData = estX_JB$S, xData = data$dX) ## NOTE: for centered X, equivalent to xData %*% sData/(px-1)

# JB on Y
estY_JB = lngca(xData = data$dY, n.comp = n.comp, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB',stand = F)
Uyfull <- estY_JB$U
My_JB = est.M.ols(sData = estY_JB$S, xData = data$dY)

test_that("linear non-Gaussian component analysis", {

  expect_equal(dim(estX_JB$U),c(n.comp,48))
  expect_equal(dim(estX_JB$M),c(48,n.comp))
  expect_equal(dim(Mx_JB),c(48,n.comp))

})


matchMxMy = suppressWarnings(greedymatch(Mx_JB, My_JB, Ux = Uxfull, Uy = Uyfull))
permJoint <- permTestJointRank(matchMxMy$Mx,matchMxMy$My) # alpha = 0.01, nperm=1000
#pval_joint = permJoint$pvalues
joint_rank = permJoint$rj

test_that("greedy match and permTest", {

  expect_equal(dim(matchMxMy$Mx),c(48,n.comp))
  expect_equal(dim(matchMxMy$Ux),c(n.comp,48))
  pval_joint = permJoint$pvalues
  expect_equal(length(pval_joint),n.comp)
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
out_indiv_small <- curvilinear_c(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, rho = rho, tol = tol, maxiter = 1500, rj = joint_rank)

test_that("Curvilinear search", {

  expect_equal(dim(out_indiv_small$Ux),c(n.comp,48))
  expect_equal(dim(out_indiv_small$Uy),c(n.comp,48))
})

