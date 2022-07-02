library(devtools)
document()
?library(Rcpp)
library(RcppArmadillo)
use_rcpp_armadillo()

install_github("thebrisklab/singR",force=TRUE)

install.packages("testthat")
library(testthat)
testthat::skip_on_cran()
use_package("ICtest")


usethis::use_data(simdata, overwrite = TRUE,compress = TRUE)
usethis::use_data(exampledata, overwrite = TRUE,compress = TRUE)

data(simdata)
data(exampledata)
usethis::use_data_raw()
exampledata = data


class(mj)

usethis::use_test("name")

devtools::test()
usethis::use_testthat()
use_test("NameOfTest")


devtools::check(env_vars = c(NOT_CRAN = "FALSE")) # cran check
devtools::test() # test that

use_citation()
library(spelling)
devtools::spell_check()
spelling::updatae_wordlist()

usethis::use_github_action_check_standard()
usethis::use_github_actions_badge(name = "R-CMD-check",
                                  repo_spec = NULL)
usethis::use_coverage()
usethis::use_github_action("test-coverage")

usethis::use_vignette("my-vignette","singR")

library(testthat)

tools::checkRdaFiles()
