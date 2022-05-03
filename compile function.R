library(devtools)
document()
?library(Rcpp)
library(RcppArmadillo)
use_rcpp_armadillo()
use_package("parallel")

install_github("thebrisklab/singR")

