test_that("bayesweight works with no errors", {
  library(R2jags)
  library(coda)

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
  weights <- bayesweight(trtmodel.list = list(
                         A1 ~ L1_1 + L2_1,
                         A2 ~ L1_2 + L2_2 + A1),
                         data = testdata,
                         n.chains = 1,
                         n.iter = 20,
                         n.burnin = 10,
                         n.thin = 1,
                         seed = 890123,
                         parallel = FALSE)

  # Check that the weights object has the expected length (equal to the number of observations in testdata)
  expect_equal(length(weights$weights), nrow(testdata))

  # Check that no NA values are present in the resulting weights
  expect_false(any(is.na(weights$weights)))

  # Check if the weights object is numeric
  expect_true(is.numeric(weights$weights))

  # Check that weights are non-negative
  expect_true(all(weights$weights >= 0))
})
