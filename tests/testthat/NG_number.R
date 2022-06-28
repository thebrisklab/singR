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

FOBIasymp(t(data$dX),k=2)
a=FOBIasymp(t(data$dX),k=3)
a=FOBIasymp(t(data$dX),k=4)
FOBIasymp(t(data$dX),k=5)

NG_number <- function(data){ #data nxp
  k_max=nrow(data)
  data=t(data)
  k=0
  FOBI=FOBIasymp(data,k=k)
  while (FOBI$p.value < 0.05) {
    k=k+1
    FOBI=FOBIasymp(data,k=k)
  }
  return(k)
}


NG_number(standard(data$dX))
dX=data$dX
n = nrow(dX)
pX = ncol(dX)
pY = ncol(dY)
dXcentered <- dX - matrix(rowMeans(dX), n, pX, byrow = F)
dYcentered <- dY - matrix(rowMeans(dY), n, pY, byrow = F)

NG_number(dXcentered)














