#' Summary function to generate result table from bayesmsm
#'
#' This function generates a ready to use result table that contents the estimated APO and ATE and their 95\% credible intervals
#'
#' @param model A model object from bayesmsm
#'
#' @return A summary table of the results from bayesmsm.
#' @importFrom stats sd quantile
#' @export
#'
#' @examples
#' # Continuous outcome
#' testdata <- read.csv(system.file("extdata",
#'                      "continuous_outcome_data.csv",
#'                      package = "bayesmsm"))
#' model <- bayesmsm(ymodel = y ~ a_1+a_2,
#'                   nvisit = 2,
#'                   reference = c(rep(0,2)),
#'                   comparator = c(rep(1,2)),
#'                   treatment_effect_type = "sq",
#'                   family = "gaussian",
#'                   data = testdata,
#'                   wmean = rep(1, 1000),
#'                   nboot = 10,
#'                   optim_method = "BFGS",
#'                   seed = 890123,
#'                   parallel = FALSE)
#' summary_bayesmsm(model)
summary_bayesmsm <- function(model) {
  # Extract bootstrapped data
  bootdata <- model$bootdata

  # Calculate summary statistics for each metric
  summary_stats <- function(metric) {
    mean_val <- mean(bootdata[[metric]])
    sd_val <- sd(bootdata[[metric]])
    quantiles <- quantile(bootdata[[metric]], probs = c(0.025, 0.975))

    return(c(mean = mean_val, sd = sd_val, `2.5%` = quantiles[1], `97.5%` = quantiles[2]))
  }

  # Initialize a list to store summary statistics
  results_list <- list()

  # Add summary statistics for reference and comparator
  results_list$Reference <- summary_stats("effect_reference")
  results_list$Comparator <- summary_stats("effect_comparator")

  # Add summary statistics for RD
  results_list$RD <- summary_stats("RD")

  # Check if RR and OR exist in the bootdata and add them if they do
  if ("RR" %in% colnames(bootdata)) {
    results_list$RR <- summary_stats("RR")
  }

  if ("OR" %in% colnames(bootdata)) {
    results_list$OR <- summary_stats("OR")
  }

  # Convert the list to a data frame for better presentation
  summary_table <- do.call(rbind, results_list)

  # Set the column names explicitly
  colnames(summary_table) <- c("mean", "sd", "2.5%", "97.5%")

  # Return the summary table
  return(summary_table)
}
