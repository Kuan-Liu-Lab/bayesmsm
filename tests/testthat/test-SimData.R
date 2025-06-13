# Define simple models for testing
amodel <- list(
  # Visit 1: intercept-only
  c("(Intercept)" = 0, "L1_1" = 0, "L2_1" = 0),
  # Visit 2: intercept-only, no dependence on previous A
  c("(Intercept)" = 0, "L1_2" = 0, "L2_2" = 0, "A_prev" = 0)
)
ymodel <- c("(Intercept)" = 0)  # constant baseline for Y

# 1) Test without right_censor: output structure and no NAs
test_that("simData returns complete data.frame without censoring", {
  dat <- simData(
    n = 50,
    n_visits = 2,
    covariate_counts = c(2, 2),
    amodel = amodel,
    ymodel = ymodel,
    y_type = "binary",
    right_censor = FALSE,
    seed = 101
  )
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 50)
  # Expect columns L1_1,L2_1,A1,L1_2,L2_2,A2,Y
  expected_cols <- c("L1_1","L2_1","A1","L1_2","L2_2","A2","Y")
  expect_true(all(expected_cols %in% names(dat)))
  expect_false(any(is.na(dat)))
})

# 2) Test reproducibility with seed
test_that("simData produces identical outputs when seeded", {
  dat1 <- simData(
    n = 30, n_visits = 2, covariate_counts = c(2,2),
    amodel = amodel, ymodel = ymodel, y_type = "binary",
    right_censor = FALSE, seed = 202
  )
  dat2 <- simData(
    n = 30, n_visits = 2, covariate_counts = c(2,2),
    amodel = amodel, ymodel = ymodel, y_type = "binary",
    right_censor = FALSE, seed = 202
  )
  expect_equal(dat1, dat2)
})

# 3) Error if right_censor=TRUE without cmodel
test_that("simData errors when right_censor=TRUE without cmodel", {
  expect_error(
    simData(
      n = 10, n_visits = 2, covariate_counts = c(2,2),
      amodel = amodel, ymodel = ymodel, y_type = "binary",
      right_censor = TRUE, seed = 303
    ),
    "Provide cmodel list"
  )
})

# 4) Test right_censor: some C2 observed and Y NAs post-censor
test_that("simData generates correct censor indicators and missingness", {
  cmodel <- list(
    # baseline censor always 0.5
    c("(Intercept)" = 0, "L1_1" = 0, "L2_1" = 0, "A" = 0),
    # second visit censor always 0
    c("(Intercept)" = -10, "L1_2" = 0, "L2_2" = 0, "A" = 0)
  )
  dat_cens <- simData(
    n = 50, n_visits = 2, covariate_counts = c(2,2),
    amodel = amodel, ymodel = ymodel, y_type = "binary",
    right_censor = TRUE, cmodel = cmodel, seed = 404
  )
  # C1 should be present and some ones
  expect_true("C1" %in% names(dat_cens))
  expect_true(any(dat_cens$C1 %in% c(0,1)))
  # C2 should be present; since cmodel2 intercept=-10, C2 should be all 0 or NA
  expect_true("C2" %in% names(dat_cens))
  # if any C1==1 then L2,A2,Y==NA; otherwise not NA
  cens1 <- which(dat_cens$C1 == 1)
  if (length(cens1)>0) {
    expect_true(all(is.na(dat_cens$L1_2[cens1])))
    expect_true(all(is.na(dat_cens$A2[cens1])))
    expect_true(all(is.na(dat_cens$Y[cens1])))
  }
  # rows with C1==0 should have C2==0
  uncens1 <- which(dat_cens$C1 == 0)
  expect_true(all(dat_cens$C2[uncens1] == 0))
})
