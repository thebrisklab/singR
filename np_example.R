set.seed(123)
library(singR)
library(fastICA)

# set simulation data
# A n x n, 100 x 100
# S n x p, 100 x 50, 3 non-gaussian,97 gaussian
# Y n x p, 100 x 50
p=50
S <- rbind(runif(p),rexp(p),rexp(p))

S.jb=apply(X = t(S),MARGIN = 2,FUN = jb.stat)
sum(S.jb[[1]]$gps)
sum(S.jb[[2]]$gps)
sum(S.jb[[3]]$gps)

Snorm <- matrix(rnorm(97*50),97,50)
S <- rbind(S,Snorm)
A <- matrix(seq(0.1,100,0.1),100,100)

Y <- A %*% S # in lngca, each row of S is an independent component.
Y.t <- t(S) %*% t(A) # in fast ICA, each column of S is an independent component. The two Y is the same.

Y.est=lngca(Y,n.comp = 10)
Y.ica=fastICA(Y.t,n.comp = 3)


pmse(S1 = t(S[1:3,]),S2 = t(Y.est$S))

pmse(S1 = Y.ica$S,S2 = t(S[1:3,]))

pmse(S1= Y.ica$S,S2 = t(Y.est$S)) # the input of S need to be p x r, so transpose for the lngca and S output.


##### set another dataset X, with the same A and different S
# A   n x n, 100 x 100
# S.2 n x p, 100 x 80
# X   n x p, 100 x 80
p.2=80
S.2 <- rbind(rexp(p.2),runif(p.2),rexp(p.2))
Snorm.2 <- matrix(rnorm(97*80),97,80)
S.2 <- rbind(S.2,Snorm.2)
A.2 <- A
X <- A.2 %*% S.2

output.2 <- singR(dX = X,dY = Y,n.comp.X = 3,n.comp.Y = 3)


pmse(M1  = t(output.2$est.Mjx),M2 = t(output.2$est.Mjy))  # the M input is r x n.
pmse(S1 = t(S[1:3,]),S2 = t(output.2$Sjy))
pmse(S1 = t(S.2[1:3,]),S2 = t(output.2$Sjx))
