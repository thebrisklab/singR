set.seed(123)
library(singR)
library(fastICA)

# set simulation data
# A n x n, 100 x 100
# S n x p, 100 x 50, 3 non-gaussian,97 gaussian
# Y n x p, 100 x 50
p=50
S <- rbind(runif(p),rexp(p),rexp(p))
Snorm <- matrix(rnorm(97*50),97,50)
S <- rbind(S,Snorm)
A <- matrix(seq(0.1,100,0.1),100,100)

Y <- A %*% S # in lngca, each row of S is an independent component.
Y <- t(S) %*% t(A) # in fast ICA, each column of S is an independent component. The two Y is the same.

Y.est=lngca(Y,n.comp = 3)
Y.ica=fastICA(Y,n.comp = 3)


pmse(S1 = S[1:3,],S2 = Y.est$S)

pmse(S1 = t(Y.ica$S),S2 = S[1:3,]) # the input of S need to be r x p, so transpose for the fastICA output.

pmse(S1=t(Y.ica$S),S2 = Y.est$S)


##### set another dataset X, with the same A and different S
p.2=80
S.2 <- rbind(rexp(p.2),runif(p.2),rexp(p.2))
Snorm.2 <- matrix(rnorm(97*80),97,80)
S.2 <- rbind(S.2,Snorm.2)
A.2 <- A
X <- A.2 %*% S.2

output.2 <- singR(dX = X,dY = Y,n.comp.X = 3,n.comp.Y = 3)
