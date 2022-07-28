set.seed(123)
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




library(ICtest)
FOBIasymp(t(Y),k=1)

data("exampledata")
NG_number(exampledata$dX)
FOBIasymp(t(exampledata$dX),k=5)
FOBIasymp.2(Y,k=100,type = 'S2')
FOBIasymp(Y,k=5)

FOBIboot(Y,k=5)
FOBI(t(Y))

FOBIasymp.2(Y,k=60)
NG_number(t(Y))

FOBIasymp.2(Y,k=99,type='S1')
NG_number.2(t(Y),type = 'S3')

NG_number.2(exampledata$dX,type = 'S1')
