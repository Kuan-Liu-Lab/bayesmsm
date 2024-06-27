#' Plot Average Treatment Effect (ATE) Density from Bootstrap Results
#'
#' @param input A model object, data frame or vector containing the bootstrap estimates of ATE.
#' @param col_density Color for the density plot (default is "blue").
#' @param fill_density Fill color for the density plot (default is "lightblue").
#' @param main Title of the plot (default is "Density of ATE Estimates").
#' @param xlab X-axis label (default is "ATE").
#' @param ylab Y-axis label (default is "Density").
#' @param xlim Limits for the x-axis (default is NULL).
#' @param ylim Limits for the y-axis (default is NULL).
#' @param ... Additional graphical parameters passed to the plot function.
#'
#' @export
#'
#' @examples
#' testdata <- read.csv(system.file("extdata", "continuous_outcome_data.csv", package = "bayesmsm"))
#' model <- bayesmsm(ymodel = y ~ a_1+a_2,
#'                            nvisit = 2,
#'                            reference = c(rep(0,2)),
#'                            comparator = c(rep(1,2)),
#'                            family = "gaussian",
#'                            data = testdata,
#'                            wmean = rep(1, 1000),
#'                            nboot = 1000,
#'                            optim_method = "BFGS",
#'                            parallel = TRUE,
#'                            ncore = 2)
#' plot_ATE(model)
#'
plot_ATE <- function(model,
                     estimand = c("RD", "RR", "OR"),
                     col_density = "blue",
                     fill_density = "lightblue",
                     main = "Posterior Predictive Distribution of Average Treatment Effect (ATE)",
                     xlab = "ATE", ylab = "Posterior Predictive Distribution",
                     xlim = NULL, ylim = NULL, ...) {

  estimand <- match.arg(estimand)

  # Check if the input model contains the necessary columns
  if (!is.list(model) || !"bootdata" %in% names(model)) {
    stop("input must be a model object containing a 'bootdata' data frame.")
  }

  if (!estimand %in% names(model$bootdata)) {
    stop(paste("bootdata within the model object must have a '", estimand, "' column.", sep = ""))
  }

  ate_values <- model$bootdata[[estimand]]

  # Calculate the density of ATE estimates
  ate_density <- density(ate_values)
  ci <- quantile(ate_values, probs = c(0.025, 0.975))
  density_ci <- density(ate_values, from = ci[1], to = ci[2])

  # Set layout to allocate space for legend
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))

  # Plot the density
  par(mar = c(5, 4, 4, 0)) # Adjust margins for the plot
  plot(ate_density, col = col_density, main = paste(main, " (", estimand, ")", sep = ""), xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
  polygon(c(density_ci$x, rev(density_ci$x)), c(rep(min(ate_density$y), length(density_ci$x)), rev(density_ci$y)), col = rgb(0, 0, 1, alpha = 0.3))
  abline(v = mean(ate_values), col = "purple", lwd = 2, lty = 3)
  abline(v = ci[1], col = "darkgreen", lty = 2)
  abline(v = ci[2], col = "darkgreen", lty = 2)

  # Plot the legend in a separate plot
  par(mar = c(5, 0, 4, 0)) # Adjust margins for the legend
  plot.new()
  legend("center", legend = c("ATE Density",
                              paste("Mean:", round(mean(ate_values), 3)),
                              paste("95% CI: [", round(ci[1], 3), ",", round(ci[2], 3), "]")),
         col = c(col_density, "purple", "darkgreen"),
         lwd = 2, lty = c(1, 3, 2),
         bty = "n")
}
