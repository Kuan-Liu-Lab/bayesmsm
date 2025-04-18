test_that("plot_ATE works with no errors", {
  # system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm")
  testdata <- read.csv(system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm"))
  # testdata <- data.frame(y = rnorm(1000), a_1 = rbinom(1000, 1, 0.45), a_2 = rbinom(1000, 1, 0.55))
  model <- bayesmsm(ymodel = y ~ a_1+a_2,
                    nvisit = 2,
                    reference = c(rep(0,2)),
                    comparator = c(rep(1,2)),
                    treatment_effect_type = "sq",
                    family = "gaussian",
                    data = testdata,
                    wmean = rep(1, 1000),
                    nboot = 10,
                    optim_method = "BFGS",
                    seed = 890123,
                    parallel = FALSE)

  expect_silent(plot_ATE(model, ATE = "RD"))


})
