test_that("plot_est_box works with no errors", {
  # system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm")
  testdata <- read.csv(system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm"))
  # testdata <- data.frame(y = rnorm(1000), a_1 = rbinom(1000, 1, 0.45), a_2 = rbinom(1000, 1, 0.55))
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
                    parallel = TRUE,
                    ncore = 2)

  expect_silent(plot_est_box(model$bootdata))

  expect_silent(plot_est_box(model))
})
