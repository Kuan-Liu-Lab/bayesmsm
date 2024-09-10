test_that("bayesweight works with no errors", {
  library(R2jags)
  library(coda)
  library(doParallel)
  library(foreach)

  testdata <- read.csv(system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm"))
  weights <- bayesweight(trtmodel.list = list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                              a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                         data = testdata,
                         n.iter = 2500,
                         n.burnin = 1500,
                         n.thin = 5,
                         n.chains = 2,
                         seed = 890123,
                         parallel = TRUE)

  # Check that the weights object is created and is of the expected class
  expect_s3_class(weights, "mcmc.list")

  # Check that the weights object has the expected length (equal to the number of observations in testdata)
  expect_equal(length(weights), nrow(testdata))

  # Check that no NA values are present in the resulting weights
  expect_false(any(is.na(weights)))

  # Check if the weights object is numeric
  expect_true(is.numeric(weights))

  # Check that weights are non-negative
  expect_true(all(weights >= 0))

})
