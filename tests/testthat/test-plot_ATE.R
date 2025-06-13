test_that("plot_ATE works with no errors", {

  testdata <- simData(
    n                = 50,
    n_visits         = 2,
    covariate_counts = c(2, 2),
    # inline treatment models
    amodel = list(
      c("(Intercept)" =  0, "L1_1" =  0.5, "L2_1" = -0.5),
      c("(Intercept)" =  0, "L1_2" =  0.5, "L2_2" = -0.5, "A_prev" = 0.3)
    ),
    # inline outcome model
    ymodel = c("(Intercept)" = 0,
               "A1"         = 0.2,
               "A2"         = 0.3,
               "L1_2"       = 0.1,
               "L2_2"       = -0.1),
    y_type       = "continuous",
    right_censor = FALSE,
    seed         = 101
  )
  model <- bayesmsm(ymodel = Y ~ A1 + A2,
                    nvisit = 2,
                    reference = c(rep(0,2)),
                    comparator = c(rep(1,2)),
                    treatment_effect_type = "sq",
                    family = "binomial",
                    data = testdata,
                    wmean = rep(1,50),
                    nboot = 10,
                    optim_method = "BFGS",
                    seed = 890123,
                    parallel = FALSE)

  expect_silent(plot_ATE(model, ATE = "RD"))
})
