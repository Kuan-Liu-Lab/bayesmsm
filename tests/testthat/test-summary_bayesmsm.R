test_that("summary_bayesmsm works with no errors", {
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

  expect_true(all(is.numeric(unlist(summary_bayesmsm(model)))))
  expect_false(anyNA(summary_bayesmsm(model)))
})
