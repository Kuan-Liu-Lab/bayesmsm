library(testthat)

test_that("plot_APO works with no errors", {
  system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm")
  testdata <- read.csv(system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm"))
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

  expect_silent(plot_APO(model, "effect_comparator"))

  expect_silent(plot_APO(model$bootdata, "effect_reference"))
})
