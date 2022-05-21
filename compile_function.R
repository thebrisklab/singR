library(devtools)
document()
?library(Rcpp)
library(RcppArmadillo)
use_rcpp_armadillo()
use_package("parallel",type = "NULL")

install_github("thebrisklab/singR")

install.packages("testthat")
library(testthat)
testthat::skip_on_cran()










usethis::use_test("name")

devtools::test()
usethis::use_testthat()

devtools::check()
devtools::test()
