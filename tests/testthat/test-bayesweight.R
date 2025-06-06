test_that("bayesweight works with no errors", {
  library(R2jags)
  library(coda)

  simdat <- read.csv(system.file("extdata",
                       "sim_causal.csv",
                       package = "bayesmsm"))
  weights <- bayesweight(trtmodel.list = list(
                         A1 ~ L11 + L21,
                         A2 ~ L11 + L21 + L12 + L22 + A1,
                         A3 ~ L11 + L21 + L12 + L22 + A1 + L13 + L23 + A2),
                         data = simdat,
                         n.chains = 1,
                         n.iter = 20,
                         n.burnin = 10,
                         n.thin = 1,
                         seed = 890123,
                         parallel = FALSE)

  # Check that the weights object has the expected length (equal to the number of observations in simdat)
  expect_equal(length(weights$weights), nrow(simdat))

  # Check that no NA values are present in the resulting weights
  expect_false(any(is.na(weights$weights)))

  # Check if the weights object is numeric
  expect_true(is.numeric(weights$weights))

  # Check that weights are non-negative
  expect_true(all(weights$weights >= 0))
})
