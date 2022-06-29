library(ICtest)

n <- 1500
S <- cbind(runif(n), rchisq(n, 2), rexp(n), rnorm(n), rnorm(n), rnorm(n))
A <- matrix(rnorm(36), ncol = 6)
X <- S %*% t(A)


FOBIasymp(X, k = 2)
FOBIasymp(X, k = 3)
FOBIasymp(X, k = 4)

library(singR)
data("exampledata")
data=exampledata


NG_number(standard(data$dX))
dX=data$dX
n = nrow(dX)
pX = ncol(dX)

dXcentered <- dX - matrix(rowMeans(dX), n, pX, byrow = F)


NG_number(dXcentered)

NG_number(data$dX)
NG_number(standard(data$dX))

data=standard(data$dX)
FOBIasymp(t(data),k=9)

