# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

# To install;
# devtools::install_github("Kuan-Liu-Lab/Bayesmsm")
library(testthat)
library(bayesmsm)

test_check("bayesmsm")

usethis::use_test()

devtools::check()
