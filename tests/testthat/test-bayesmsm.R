library(testthat)

test_that("bayesmsm works with no errors", {
  system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm")
  testdata <- read.csv(system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm"))
  expect_silent(
    model <- bayesmsm(ymodel = y ~ a_1+a_2,
                      nvisit = 2,
                      reference = c(rep(0,2)),
                      comparator = c(rep(1,2)),
                      family = "gaussian",
                      data = testdata,
                      wmean = rep(1, 1000),
                      nboot = 1000,
                      optim_method = "BFGS",
                      estimand = 'RD',
                      parallel = FALSE,
                      ncore = 6)
  )
})


# > test_check("bayesmsm")
# [ FAIL 1 | WARN 0 | SKIP 0 | PASS 4 ]
#
# ══ Failed tests ════════════════════════════════════════════════════════════════
# ── Failure ('test-bayesmsm.R:6:3'): bayesmsm works with no errors ──────────────
# `... <- NULL` produced messages.
#
# [ FAIL 1 | WARN 0 | SKIP 0 | PASS 4 ]
# Error: Test failures
# Execution halted
#
# 1 error ✖ | 1 warning ✖ | 1 note ✖
