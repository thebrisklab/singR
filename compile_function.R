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



usethis::use_data(simdata, overwrite = TRUE)
load("c:/Software/Data/Small_Simulated_data.Rdata")
data("Small_Simulated_data.Rdata")
data(simdata)

class(mj)

usethis::use_test("name")

devtools::test()
usethis::use_testthat()

devtools::check()
devtools::test()

simdata <- list(dX=dX,dY=dY,mj=mj,sIx=new_sIx,sIy=new_sIy,sjx=new_sjx,sjy=new_sjy)
