test_that("bayesweight_cen works with no errors", {
  library(R2jags)
  library(doParallel)

  simdat_cen <- read.csv(system.file("extdata", "sim_causal.csv", package = "bayesmsm"))
  weights_cen <- bayesweight_cen(trtmodel.list = list(A1 ~ L11 + L21,
                                                      A2 ~ L11 + L21 + L12 + L22 + A1,
                                                      A3 ~ L11 + L21 + L12 + L22 + A1 + L13 + L23 + A2),
                                 cenmodel.list = list(C1 ~ L11 + L21,
                                                      C2 ~ L11 + L21 + A1,
                                                      C3 ~ L11 + L21 + A1 + L12 + L22 + A2),
                                 data = simdat_cen,
                                 n.iter = 1500,
                                 n.burnin = 500,
                                 n.thin = 1,
                                 parallel = TRUE,
                                 n.chains = 2,
                                 seed = 890123)

  # Check that the weights_cen object has the expected dimensions (length or rows matching the data size)
  expect_equal(length(weights_cen), nrow(simdat_cen))

  # Check if the weights_cen object is numeric
  expect_true(is.numeric(weights_cen))

  # Check that weights are non-negative
  expect_true(all(na.omit(weights_cen) >= 0))

})
