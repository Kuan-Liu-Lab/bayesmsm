#' Bayesian Treatment Effect Weight Estimation for Censored Data
#'
#' This function estimates Bayesian importance sampling weights for treatment models and censoring models across multiple time points via JAGS
#'
#' @param trtmodel.list A list of formulas corresponding to each time point with the time-specific treatment variable on the left-hand side and pre-treatment covariates to be balanced on the right-hand side. The formulas must be in temporal order, and must contain all covariates to be balanced at that time point. Interactions and functions of covariates are allowed.
#' @param cenmodel.list A list of formulas for the censored data at each time point, with censoring indicators on the left-hand side and covariates on the right-hand side. The formulas must be in temporal order, and must contain all covariates to be balanced at that time point.
#' @param data A data set in the form of a data frame containing the variables in "trtmodel.list" and "cenmodel.list". This must be a wide data set with exactly one row per unit.
#' @param n.chains Integer specifying the number of MCMC chains to run. Set to 1 for non-parallel computation. For parallel computation, it is required to use at least 2 chains. The default is 2.
#' @param n.iter Integer specifying the total number of iterations for each chain (including burn-in). The default is 25000.
#' @param n.burnin Integer specifying the number of burn-in iterations for each chain. The default is 15000.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampler. The default is 5.
#' @param seed Starting seed for the JAGS model. The default is NULL.
#' @param parallel Logical scalar indicating whether to run the MCMC chains in parallel. The default is TRUE.
#'
#' @return A list of the calculated weights and the JAGS model where "weights" is a vector of posterior mean weights, computed by taking the average of the weights across all MCMC iterations and `model_string` is a character of the JAGS model based on the input of "trtmodel.list".
#'
#' @importFrom R2jags jags
#' @importFrom coda mcmc as.mcmc geweke.diag
#' @import parallel
#' @import doParallel
#' @importFrom foreach "%dopar%"
#' @importFrom stats as.formula terms var
#'
#' @export
#'
#' @examples
#' simdat_cen <- read.csv(system.file("extdata",
#'                        "sim_causal.csv",
#'                        package = "bayesmsm"))
#' weights_cen <- bayesweight_cen(
#'                 trtmodel.list = list(
#'                 A1 ~ L11 + L21,
#'                 A2 ~ L11 + L21 + L12 + L22 + A1,
#'                 A3 ~ L11 + L21 + L12 + L22 + A1 + L13 + L23 + A2),
#'                 cenmodel.list = list(
#'                 C1 ~ L11 + L21,
#'                 C2 ~ L11 + L21 + A1,
#'                 C3 ~ L11 + L21 + A1 + L12 + L22 + A2),
#'                 data = simdat_cen,
#'                 n.chains = 1,
#'                 n.iter = 20,
#'                 n.burnin = 10,
#'                 n.thin = 1,
#'                 seed = 890123,
#'                 parallel = FALSE)
#' summary(weights_cen)
#' simdat_cen <- read.csv("censoring_data.csv", header = TRUE, fileEncoding="UTF-8-BOM")
#' weights_cen <- bayesweight_cen(
#'                 trtmodel.list = list(a_1 ~ w1 + w2 + L1_1 + L2_1,
#'                 a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
#'                 cenmodel.list = list(c ~ w1 + w2 + L1_1 + L2_1 + a_1),
#'                 data = simdat_cen,
#'                 n.chains = 1,
#'                 n.iter = 20,
#'                 n.burnin = 10,
#'                 n.thin = 1,
#'                 seed = 890123,
#'                 parallel = FALSE)
#' summary(weights_cen)
bayesweight_cen <- function(trtmodel.list,
                            cenmodel.list,
                            data,
                            n.chains = 2,
                            n.iter = 25000,
                            n.burnin = 15000,
                            n.thin = 5,
                            seed = NULL,
                            parallel = TRUE) {

  create_marginal_treatment_models <- function(trtmodel.list) {

    # Use lapply to iterate over indices
    trtmodel.list_s <- lapply(seq_along(trtmodel.list), function(index) {
      # Extract the response variable (treatment variable) from each model
      response_var <- all.vars(trtmodel.list[[index]])[1]  # assuming the response is the first variable on the LHS

      # Create the marginal model formula
      if (index == 1) {
        # The first treatment model does not depend on any previous treatments
        formula_s <- as.formula(paste(response_var, "~ 1"))
      } else {
        # Subsequent treatment models depend on all previous treatments
        previous_treatments <- sapply(seq_len(index - 1), function(j) {
          all.vars(trtmodel.list[[j]])[1]
        })
        formula_s <- as.formula(paste(response_var, "~", paste(previous_treatments, collapse = " + ")))
      }

      return(formula_s)
    })

    return(trtmodel.list_s)
  }

  # Generate trtmodel.list_s
  trtmodel.list_s <- create_marginal_treatment_models(trtmodel.list)
  cenmodel.list_s <- create_marginal_treatment_models(cenmodel.list)

  # Extract variables from treatment and censoring models
  extract_variables_list <- function(input) {
    extract_from_formula <- function(formula) {
      formula_terms <- terms(formula)
      response_variable <- attr(formula_terms, "response")
      response_name <- if (response_variable > 0) {
        all_vars <- all.vars(formula)
        all_vars[response_variable]
      } else {NA}
      predictor_names <- attr(formula_terms, "term.labels")
      list(response = response_name, predictors = predictor_names)
    }
    if (is.list(input)) {
      lapply(input, extract_from_formula)
    } else {
      extract_from_formula(input)
    }
  }

  trtmodel <- extract_variables_list(trtmodel.list)
  trtmodel_s <- extract_variables_list(trtmodel.list_s)
  cenmodel <- extract_variables_list(cenmodel.list)
  cenmodel_s <- extract_variables_list(cenmodel.list_s)

  # Define JAGS model for treatment and censoring
  write_jags_model <- function(trtmodel.list, cenmodel.list) {

    n_trt <- length(trtmodel.list)
    n_cen <- length(cenmodel.list)
    shift <- n_trt - n_cen      # how many early visits have no cen-model

    var_info_trt <- lapply(trtmodel.list, extract_variables_list)
    var_info_cen <- lapply(cenmodel.list, extract_variables_list)

    model_string <- "model{\n"

    for (v in seq_along(var_info_trt)) {

      cen_idx <- v - shift
      has_cen <- cen_idx >= 1 && cen_idx <= n_cen

      visit_trt <- var_info_trt[[v]]
      response_trt <- visit_trt$response
      predictors_trt <- visit_trt$predictors

      # visit_cen <- var_info_cen[[v]]
      # response_cen <- visit_cen$response
      # predictors_cen <- visit_cen$predictors
      if (has_cen) {
        visit_cen    <- var_info_cen[[cen_idx]]
        response_cen <- visit_cen$response
        predictors_cen <- visit_cen$predictors
      }

      model_string <- paste0(model_string, "\nfor (i in 1:N", v, ") {\n")

      # Conditional treatment assignment model
      model_string <- paste0(model_string,
                             "\n# conditional model;\n",
                             response_trt, "[i] ~ dbern(p", v, "[i])\n",
                             "logit(p", v, "[i]) <- b", v, "0")
      for (p in seq_along(predictors_trt)) {
        model_string <- paste0(model_string, " + b", v, p, "*", predictors_trt[p], "[i]")
      }
      model_string <- paste0(model_string, "\n")

      # Censoring model
      # model_string <- paste0(model_string,
      #                        response_cen, "[i] ~ dbern(cp", v, "[i])\n",
      #                        "logit(cp", v, "[i]) <- s", v, "0")
      # for (p in seq_along(predictors_cen)) {
      #   model_string <- paste0(model_string, " + s", v, p, "*", predictors_cen[p], "[i]")
      # }
      # model_string <- paste0(model_string, "\n")
      if (has_cen) {
        model_string <- paste0(model_string,
                               response_cen,"[i] ~ dbern(cp",v,"[i])\n",
                               "logit(cp",v,"[i]) <- s",v,"0")
        for (p in seq_along(predictors_cen))
          model_string <- paste0(model_string," + s",v,p,"*",predictors_cen[p],"[i]")
        model_string <- paste0(model_string,"\n")
      }

      # Marginal treatment assignment model
      model_string <- paste0(model_string,
                             "\n# marginal model;\n",
                             response_trt, "s[i] ~ dbern(p", v, "s[i])\n",
                             "logit(p", v, "s[i]) <- bs", v, "0")
      if (v > 1) {
        for (j in 1:(v - 1)) {
          prev_response_trt <- var_info_trt[[j]]$response
          model_string <- paste0(model_string, " + bs", v, j, "*", prev_response_trt, "s[i]")
        }
      }
      model_string <- paste0(model_string, "\n")

      # Marginal censoring model
      # model_string <- paste0(model_string,
      #                        response_cen, "s[i] ~ dbern(cp", v, "s[i])\n",
      #                        "logit(cp", v, "s[i]) <- ts", v, "0")
      # if (v > 1) {
      #   for (j in 1:(v - 1)) {
      #     prev_response_trt <- var_info_trt[[j]]$response
      #     model_string <- paste0(model_string, " + ts", v, j, "*", prev_response_trt, "s[i]")
      #   }
      # }
      # model_string <- paste0(model_string, "\n}\n")
      if (has_cen) {
        model_string <- paste0(model_string,
                               response_cen,"s[i] ~ dbern(cp",v,"s[i])\n",
                               "logit(cp",v,"s[i]) <- ts",v,"0")
        if (v > 1) {
          for (j in 1:(v-1)) {
            prev_response_trt <- var_info_trt[[j]]$response
            model_string <- paste0(model_string,
                                   " + ts",v,j,"*",prev_response_trt,"s[i]")
          }
        }
        model_string <- paste0(model_string,"\n")
      }
      model_string <- paste0(model_string,"}\n")

    }

    # Priors section
    model_string <- paste0(model_string, "\n# Priors\n")
    for (v in seq_along(var_info_trt)) {

      cen_idx <- v - shift
      has_cen <- cen_idx >= 1 && cen_idx <= n_cen

      num_preds_trt <- length(var_info_trt[[v]]$predictors)
      if (has_cen) {
        num_preds_cen <- length(var_info_cen[[cen_idx]]$predictors)
      }

      # Treatment priors
      for (p in 0:num_preds_trt) {
        model_string <- paste0(model_string, "b", v, p, " ~ dunif(-10, 10)\n")
      }

      # Censoring priors
      # for (p in 0:num_preds_cen) {
      #   model_string <- paste0(model_string, "s", v, p, " ~ dunif(-10, 10)\n")
      # }
      if (has_cen) {
        for (p in 0:num_preds_cen)
          model_string <- paste0(model_string,"s", v, p," ~ dunif(-10, 10)\n")
      }

      # Marginal treatment priors
      model_string <- paste0(model_string, "bs", v, "0 ~ dunif(-10, 10)\n")
      if (v > 1) {
        for (j in 1:(v - 1)) {
          model_string <- paste0(model_string, "bs", v, j, " ~ dunif(-10, 10)\n")
        }
      }

      # Marginal censoring priors
      # model_string <- paste0(model_string, "ts", v, "0 ~ dunif(-10, 10)\n")
      # if (v > 1) {
      #   for (j in 1:(v - 1)) {
      #     model_string <- paste0(model_string, "ts", v, j, " ~ dunif(-10, 10)\n")
      #   }
      # }
      if (has_cen) {
        model_string <- paste0(model_string,"ts", v, "0 ~ dunif(-10, 10)\n")
        if (v > 1) {
          for (j in 1:(v - 1)) {
            model_string <- paste0(model_string,"ts", v, j, " ~ dunif(-10, 10)\n")
          }
        }
      }
    }

    # Add the closing brace for the model block
    model_string <- paste0(model_string, "}\n")

    return(model_string)
  }

  # extracting parameters list;
  jags_model_parameter <- function(trtmodel.list,
                                   cenmodel.list){

    n_trt <- length(trtmodel.list)
    n_cen <- length(cenmodel.list)
    shift <- n_trt - n_cen

    all_parameters <- c()
    var_info_trt <- lapply(trtmodel.list, extract_variables_list)
    var_info_cen <- lapply(cenmodel.list, extract_variables_list)

    for (v in seq_along(var_info_trt)) {

      cen_idx <- v - shift
      has_cen <- cen_idx >= 1 && cen_idx <= n_cen

      visit_trt <- var_info_trt[[v]]
      response_trt <- visit_trt$response
      predictors_trt <- visit_trt$predictors

      if (has_cen) {
        visit_cen <- var_info_cen[[cen_idx]]
        response_cen <- visit_cen$response
        predictors_cen <- visit_cen$predictors
      }

      # Treatment model
      for (p in seq_along(predictors_trt)) {
        all_parameters <- c(all_parameters, sprintf("b%d%d", v, p))
      }

      # Censoring model
      if (has_cen) {
        for (p in seq_along(predictors_cen)) {
          all_parameters <- c(all_parameters, sprintf("s%d%d", v, p))
        }
      }

      all_parameters <- c(all_parameters, sprintf("bs%d0", v))

      if (v > 1) {
        for (j in 1:(v - 1)) {
          prev_response_trt <- var_info_trt[[j]]$response
          all_parameters <- c(all_parameters, sprintf("bs%d%d", v, j))
        }
      }

      if (has_cen) {
        all_parameters <- c(all_parameters, sprintf("ts%d0", v))
        if (v > 1) {
          for (j in 1:(v - 1)) {
            prev_response_trt <- var_info_trt[[j]]$response
            all_parameters <- c(all_parameters, sprintf("ts%d%d", v, j))
          }
        }
      }
    }

    # Priors section
    for (v in seq_along(var_info_trt)) {

      cen_idx <- v - shift
      has_cen <- cen_idx >= 1 && cen_idx <= n_cen

      num_preds_trt <- length(var_info_trt[[v]]$predictors)
      if (has_cen) {
        num_preds_cen <- length(var_info_cen[[cen_idx]]$predictors)
      }

      # Treatment priors
      for (p in 0:num_preds_trt) {
        all_parameters <- c(all_parameters, sprintf("b%d%d", v, p))
      }

      # Censoring priors
      if (has_cen) {
        for (p in 0:num_preds_cen) {
          all_parameters <- c(all_parameters, sprintf("s%d%d", v, p))
        }
      }

      # Marginal treatment priors
      all_parameters <- c(all_parameters, sprintf("bs%d0", v))
      if (v > 1) {
        for (j in 1:(v - 1)) {
          all_parameters <- c(all_parameters, sprintf("bs%d%d", v, j))
        }
      }

      # Marginal censoring priors
      if (has_cen) {
        all_parameters <- c(all_parameters, sprintf("ts%d0", v))
        if (v > 1) {
          for (j in 1:(v - 1)) {
            all_parameters <- c(all_parameters, sprintf("ts%d%d", v, j))
          }
        }
      }
    }

    return(unique(all_parameters))
  }

  # Prepare data for JAGS
  prepare_jags_data <- function(data, trtmodel.list, cenmodel.list) {

    n_trt  <- length(trtmodel.list)
    n_cen  <- length(cenmodel.list)
    shift  <- n_trt - n_cen

    variable_info_trt <- lapply(trtmodel.list, extract_variables_list)
    variable_info_cen <- lapply(cenmodel.list, extract_variables_list)

    # Collect all variables from trtmodel and cenmodel
    all_vars <- unique(unlist(lapply(c(variable_info_trt, variable_info_cen), function(info) {
      c(info$response, info$predictors)
    })))
    if (n_cen > 0) {
      cen_vars <- unique(unlist(lapply(variable_info_cen, \(info)
                                       c(info$response, info$predictors))))
      all_vars <- unique(c(all_vars, cen_vars))
    }

    # Initialize the list to pass to JAGS
    jags.data <- list()

    # Add each column of data to the jags.data list, with NA values removed
    for (col in all_vars) {
      if (!is.na(col)) {
        jags.data[[col]] <- data[[col]][!is.na(data[[col]])]
      }
    }

    # Determine N for each visit based on the length of non-NA values of the response variable
    for (v in seq_along(variable_info_trt)) {
      response_trt <- variable_info_trt[[v]]$response
      jags.data[[paste0("N", v)]] <- sum(!is.na(data[[response_trt]]))

      # Marginal Treatment History
      jags.data[[paste0(response_trt, "s")]] <- data[[response_trt]][!is.na(data[[response_trt]])] # Remove the NAs

      # Marginal Censoring History
      # jags.data[[paste0("C", v, "s")]] <- data[[paste0("C", v)]][!is.na(data[[paste0("C", v)]])]
      cen_idx <- v - shift
      has_cen <- cen_idx >= 1 && cen_idx <= n_cen
      if (has_cen) {
        response_cen <- variable_info_cen[[cen_idx]]$response
        jags.data[[paste0(response_cen, "s")]] <-
          data[[response_cen]][!is.na(data[[response_cen]])]
      }
    }

    return(jags.data)
  }

  # Generate the model string
  model_str <- write_jags_model(trtmodel.list, cenmodel.list)

  # Prepare JAGS data
  jags.data <- prepare_jags_data(data, trtmodel.list, cenmodel.list)
  jags.params <- jags_model_parameter(trtmodel.list, cenmodel.list)

  # Define seed for each chain
  if(!is.null(seed)){
    set.seed(seed)
  }
  new.seed.by.chain <- sample.int(2^30, size=n.chains)

  # Run JAGS model
  if (parallel == TRUE) {
    if (n.chains == 1) {
      stop("Parallel MCMC requires at least 2 chains. Computing is running on 1 core per chain.")
    }
    available_cores <- detectCores(logical = FALSE)
    if (n.chains >= available_cores) {
      stop(paste("Parallel MCMC requires 1 core per chain. You have", available_cores, "cores. We recommend using", available_cores - 2, "cores."))
    }
    # Run JAGS model in parallel
    cl <- parallel::makeCluster(n.chains)
    doParallel::registerDoParallel(cl)

    # jags.model.wd <- paste(getwd(), '/censoring_model.txt',sep='')

    posterior <- foreach::foreach(chain_idx=1:n.chains, .packages=c('R2jags'),
                                  .combine='rbind') %dopar%{

                                    con <- textConnection(model_str)
                                    jagsfit <- jags(data = jags.data,
                                                    parameters.to.save = jags.params,
                                                    model.file = con,
                                                    n.chains = 1,
                                                    n.iter = n.iter,
                                                    n.burnin = n.burnin,
                                                    n.thin = n.thin,
                                                    jags.seed = new.seed.by.chain[chain_idx])
                                    close(con)

                                    # Combine MCMC output from multiple chains
                                    out.mcmc <- as.mcmc(jagsfit)
                                    return(do.call(rbind, lapply(out.mcmc, as.matrix)))

                                  }

    parallel::stopCluster(cl)

  } else if (parallel == FALSE) {

    if (n.chains != 1) {
      stop("Non-parallel MCMC requires exactly 1 chain.")
    }

    con <- textConnection(model_str)
    # Run JAGS model without parallel computing
    jagsfit <- jags(data = jags.data,
                    parameters.to.save = jags.params,
                    model.file = con,
                    n.chains = 1,
                    n.iter = n.iter,
                    n.burnin = n.burnin,
                    n.thin = n.thin,
                    jags.seed = new.seed.by.chain[1])
    close(con)

    # Extract MCMC output
    out.mcmc <- as.mcmc(jagsfit)
    diagnostics <- geweke.diag(out.mcmc)

    # Check diagnostics for convergence issues
    significant_indices <- which(abs(diagnostics[[1]]$z) > 1.96)
    if (length(significant_indices) > 0) {
      warning("Some parameters have not converged with Geweke index > 1.96. More iterations may be needed.")
    }

    posterior <- as.matrix(out.mcmc[[1]])

  }

  expit <- function(x){exp(x) / (1+exp(x))}


  # Initialize arrays for storing probabilities
  n_visits <- length(trtmodel)
  n_posterior <- dim(posterior)[1]
  # Calculate the number of observations for the last visit (i.e. complete data)
  n_obs <- dim(data)[1]

  # psc <- array(dim = c(n_visits, n_posterior, n_obs))
  # psm <- array(dim = c(n_visits, n_posterior, n_obs))
  # csc <- array(dim = c(n_visits, n_posterior, n_obs))
  # csm <- array(dim = c(n_visits, n_posterior, n_obs))
  psc <- array(NA_real_, dim = c(n_visits, n_posterior, n_obs))
  psm <- array(NA_real_, dim = c(n_visits, n_posterior, n_obs))
  csc <- array(1,        dim = c(n_visits, n_posterior, n_obs))  # ← 1 by default
  csm <- array(1,        dim = c(n_visits, n_posterior, n_obs))  # ← 1 by default

  # Calculate probabilities for each visit
  parameter_map <- colnames(posterior)

  for (nvisit in 1:n_visits) {

    ## ─── which censoring model (if any) belongs to this visit ─────
    n_trt <- length(trtmodel.list)
    n_cen <- length(cenmodel.list)
    shift <- n_trt - n_cen
    cen_idx <- nvisit - shift                   # mapped index into cen-list
    has_cen <- cen_idx >= 1 && cen_idx <= n_cen # TRUE if censoring exists

    predictors_c <- trtmodel[[nvisit]]$predictors
    predictors_s <- trtmodel_s[[nvisit]]$predictors
    # predictors_cen_c <- cenmodel[[nvisit]]$predictors
    # predictors_cen_s <- cenmodel_s[[nvisit]]$predictors

    design_matrix_c <- cbind(1, data[predictors_c])
    design_matrix_s <- cbind(1, data[predictors_s])
    # design_matrix_cen_c <- cbind(1, data[predictors_cen_c])
    # design_matrix_cen_s <- cbind(1, data[predictors_cen_s])

    beta_indices_c <- match(c(sprintf("b%d0", nvisit),
                              sapply(1:length(predictors_c), function(p) sprintf("b%d%d", nvisit, p))),
                            parameter_map)

    beta_indices_s <- match(c(sprintf("bs%d0", nvisit),
                              if (nvisit > 1) sapply(1:length(predictors_s), function(p) sprintf("bs%d%d", nvisit, p)) else NULL),
                            parameter_map)

    # beta_indices_cen_c <- match(c(sprintf("s%d0", nvisit),
    #                               sapply(1:length(predictors_cen_c), function(p) sprintf("s%d%d", nvisit, p))),
    #                             parameter_map)
    #
    # beta_indices_cen_s <- match(c(sprintf("ts%d0", nvisit),
    #                               if (nvisit > 1) sapply(1:length(predictors_cen_s), function(p) sprintf("ts%d%d", nvisit, p)) else NULL),
    #                             parameter_map)

    if (has_cen) {
      predictors_cen_c <- cenmodel[[cen_idx]]$predictors
      predictors_cen_s <- cenmodel_s[[cen_idx]]$predictors

      design_matrix_cen_c <- cbind(1, data[predictors_cen_c])
      design_matrix_cen_s <- cbind(1, data[predictors_cen_s])

      # beta_indices_cen_c <- match(c(sprintf("s%d0", nvisit),
      #                               sapply(1:length(predictors_cen_c), function(p) sprintf("s%d%d", nvisit, p))),
      #                             parameter_map)
      beta_indices_cen_c <- match(c(sprintf("s%d0", nvisit),
                        if (length(predictors_cen_c))
                          sprintf("s%d%d",  nvisit,
                                  seq_along(predictors_cen_c))),
                      parameter_map)

      # beta_indices_cen_s <- match(c(sprintf("ts%d0", nvisit),
      #                               if (nvisit > 1) sapply(1:length(predictors_cen_s), function(p) sprintf("ts%d%d", nvisit, p)) else NULL),
      #                             parameter_map)
      beta_indices_cen_s <- match(c(sprintf("ts%d0", nvisit),
                        if (nvisit > 1 && length(predictors_cen_s))
                          sprintf("ts%d%d", nvisit,
                                  seq_along(predictors_cen_s))),
                      parameter_map)
    }

    for (j in 1:n_posterior) {
      psc[nvisit, j, ] <- expit(posterior[j, beta_indices_c] %*% t(design_matrix_c))
      psm[nvisit, j, ] <- expit(posterior[j, beta_indices_s] %*% t(design_matrix_s))

      if (has_cen) {
        csc[nvisit, j, ] <- expit(posterior[j, beta_indices_cen_c] %*% t(design_matrix_cen_c))
        csm[nvisit, j, ] <- expit(posterior[j, beta_indices_cen_s] %*% t(design_matrix_cen_s))
      }
    }
  }

  numerator_trt <- apply(psm, c(2, 3), prod)
  denominator_trt <- apply(psc, c(2, 3), prod)
  numerator_cen <- apply(csm, c(2, 3), prod)
  denominator_cen <- apply(csc, c(2, 3), prod)

  weights <- (numerator_trt / denominator_trt) * (numerator_cen / denominator_cen)
  wmean <- colMeans(weights)

  # Return the weights and the model string
  return(list(weights = wmean, model_string = model_str))

}
