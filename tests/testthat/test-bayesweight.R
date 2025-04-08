test_that("bayesweight works with no errors", {
  library(R2jags)
  library(coda)

  testdata <- read.csv(system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm"))
  weights <- bayesweight(trtmodel.list = list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                              a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                         data = testdata,
                         n.chains = 1,
                         n.iter = 100,
                         n.burnin = 50,
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
