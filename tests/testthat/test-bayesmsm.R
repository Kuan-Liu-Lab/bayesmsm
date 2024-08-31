test_that("bayesmsm works with no errors", {
  library(MCMCpack)
  library(doParallel)
  library(foreach)
  # system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm")
  testdata <- read.csv(system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm"))
  # testdata <- data.frame(y = rnorm(1000), a_1 = rbinom(1000, 1, 0.45), a_2 = rbinom(1000, 1, 0.55))
  model <- bayesmsm(ymodel = y ~ a_1+a_2,
                      nvisit = 2,
                      reference = c(rep(0,2)),
                      comparator = c(rep(1,2)),
                      family = "gaussian",
                      data = testdata,
                      wmean = rep(1, 1000),
                      nboot = 100,
                      optim_method = "BFGS",
<<<<<<< HEAD
=======
                      # estimand = 'RD',
>>>>>>> c02a61c92d3be0eeca4a22bb0196a74db590a6ec
                      seed = 890123,
                      parallel = TRUE,
                      ncore = 2)

  # Check if 'model' is a list
  expect_true(is.list(model))

<<<<<<< HEAD
=======
  # Check if 'model' has 6 elements
  expect_length(model, 6)

  # Check if 'model' has the correct names
  # Kuan Aug30: names has to be updated to 'RD_mean', 'RD_sd', 'RD_quantile', 'bootdata', 'reference', 'comparator';
  expected_names <- c('RD_mean', 'RD_sd', 'RD_quantile', 'bootdata', 'reference', 'comparator')
  expect_named(model, expected_names)

  # Check the types of the elements within 'model'
  # expect_type(model$mean, "double")
  # expect_type(model$sd, "double")
  # expect_type(model$quantile, "double")
  # expect_true(is.data.frame(model$bootdata))
  # expect_type(model$reference, "double")
  # expect_type(model$comparator, "double")


>>>>>>> c02a61c92d3be0eeca4a22bb0196a74db590a6ec
  testdata2 <- read.csv(system.file("extdata", "binary_outcome_data.csv", package = "bayesmsm"))
  model2 <- bayesmsm(ymodel = y ~ a_1+a_2,
                    nvisit = 2,
                    reference = c(rep(0,2)),
                    comparator = c(rep(1,2)),
                    family = "binomial",
                    data = testdata2,
                    wmean = rep(1, 1000),
                    nboot = 100,
                    optim_method = "BFGS",
<<<<<<< HEAD
=======
                    # estimand = 'OR',
>>>>>>> c02a61c92d3be0eeca4a22bb0196a74db590a6ec
                    seed = 890123,
                    parallel = TRUE,
                    ncore = 2)

  # Check if 'model' is a list
  expect_true(is.list(model2))

<<<<<<< HEAD
=======
  # Check if 'model' has 12 elements
  expect_length(model2, 12)

  # Check if 'model' has the correct names
  # Kuan Aug30, update this! you added the RD_thing;
  # expected_names <- c("mean", "sd", "quantile", "bootdata", "reference", "comparator")
  # expect_named(model2, expected_names)

  # Check the types of the elements within 'model'
  # double checks didn't pass in this test file (kuan aug 30)
  # expect_type(model2$mean, "double")
  # expect_type(model2$sd, "double")
  # expect_type(model2$quantile, "double")
  # expect_true(is.data.frame(model2$bootdata))
  # expect_type(model2$reference, "double")
  # expect_type(model2$comparator, "double")
>>>>>>> c02a61c92d3be0eeca4a22bb0196a74db590a6ec

})
