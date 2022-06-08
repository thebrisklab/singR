library(devtools)
document()
?library(Rcpp)
library(RcppArmadillo)
use_rcpp_armadillo()

install_github("thebrisklab/singR")

install.packages("testthat")
library(testthat)
testthat::skip_on_cran()



usethis::use_data(simdata, overwrite = TRUE)

data(simdata)
usethis::use_data_raw()


class(mj)

usethis::use_test("name")

devtools::test()
usethis::use_testthat()
use_test("NameOfTest")


devtools::check()
devtools::test()

use_citation()
spell_check()
