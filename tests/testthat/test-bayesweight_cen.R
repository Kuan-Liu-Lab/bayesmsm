test_that("bayesweight_cen works with no errors", {
  library(R2jags)

  # 1) Define coefficient lists for 2 visits
  amodel <- list(
    # Visit 1: logit P(A1=1) = -0.3 + 0.4*L1_1 - 0.2*L2_1
    c("(Intercept)" = -0.3, "L1_1" = 0.4, "L2_1" = -0.2),
    # Visit 2: logit P(A2=1) = -0.1 + 0.3*L1_2 - 0.1*L2_2 + 0.5*A_prev
    c("(Intercept)" = -0.1, "L1_2" = 0.3, "L2_2" = -0.1, "A_prev" = 0.5)
  )

  # 2) Define outcome model: logistic on treatments and last covariates
  ymodel <- c(
    "(Intercept)" = -0.8,
    "A1"          = 0.2,
    "A2"          = 0.4,
    "L1_2"        = 0.3,
    "L2_2"        = -0.3
  )

  # 3) Define right-censoring models at each visit
  cmodel <- list(
    # Censor at visit 1 based on baseline covariates and A1
    c("(Intercept)" = -1.5, "L1_1" = 0.2, "L2_1" = -0.2, "A" = 0.2),
    # Censor at visit 2 based on visit-2 covariates and A2
    c("(Intercept)" = -1.5, "L1_2" = 0.1, "L2_2" = -0.1, "A" = 0.3)
  )

  # 4) Load package and simulate data
  testdata <- simData(
    n                = 50,
    n_visits         = 2,
    covariate_counts = c(2, 2),
    amodel           = amodel,
    ymodel           = ymodel,
    y_type           = "binary",
    right_censor     = TRUE,
    cmodel           = cmodel,
    seed             = 123
  )


  weights_cen <- bayesweight_cen(
                  trtmodel.list = list(
                  A1 ~ L1_1 + L2_1,
                  A2 ~ L1_2 + L2_2 + A1),
                  cenmodel.list = list(
                    C1 ~ L1_1 + L2_1 + A1,
                    C2 ~ L1_2 + L2_2 + A2),
                  data = testdata,
                  n.chains = 1,
                  n.iter = 20,
                  n.burnin = 10,
                  n.thin = 1,
                  seed = 890123,
                  parallel = FALSE)

  # Check that the weights_cen object has the expected dimensions (length or rows matching the data size)
  expect_equal(length(weights_cen$weights), nrow(testdata))

  # Check if the weights_cen object is numeric
  expect_true(is.numeric(weights_cen$weights))

  # Check that weights are non-negative
  expect_true(all(na.omit(weights_cen$weights) >= 0))
})
