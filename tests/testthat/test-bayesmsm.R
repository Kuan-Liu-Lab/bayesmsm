test_that("bayesmsm works with no errors", {
  library(MCMCpack)
  library(doParallel)
  library(foreach)
  # system.file("extdata", "sim_causal.csv", package = "bayesmsm")
  testdata <- read.csv(system.file("extdata",
                                   "sim_causal.csv",
                                   package = "bayesmsm"))
  model <- bayesmsm(ymodel = Y ~ A1 + A2 + A3,
                    nvisit = 3,
                    reference = c(rep(0,3)),
                    comparator = c(rep(1,3)),
                    treatment_effect_type = "sq",
                    family = "binomial",
                    data = testdata,
                    wmean = rep(1,500),
                    nboot = 10,
                    optim_method = "BFGS",
                    seed = 890123,
                    parallel = FALSE)

  # Check if 'model' is a list
  expect_true(is.list(model))

  # Check if 'model' has 12 elements
  expect_length(model, 12)

  # Check if 'model' has the correct names
  expected_names <- c('RD_mean', 'RR_mean', 'OR_mean',
                      'RD_sd', 'RR_sd', 'OR_sd',
                      'RD_quantile', 'RR_quantile', 'OR_quantile',
                      'bootdata', 'reference', 'comparator')
  expect_named(model, expected_names)
})
