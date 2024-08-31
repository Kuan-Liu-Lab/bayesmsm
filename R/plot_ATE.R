#' Plot Average Treatment Effect Density from Bootstrap Results
#'
#' @param input A model object, data frame or vector containing the bootstrap estimates of ATE.
#' @param ATE define causal estimand of interest from RD, OR, RR.
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
#'                            nboot = 100,
#'                            optim_method = "BFGS",
#'                            parallel = TRUE,
#'                            ncore = 2)
#' plot_ATE(model)
#'
plot_ATE <- function(input,
                     ATE = "RD",
                     col_density = "blue",
                     fill_density = "lightblue",
                     main = "Posterior Predictive Distribution of Average Treatment Effect",
                     xlab = "ATE", ylab = "Posterior Predictive Distribution",
                     xlim = NULL, ylim = NULL, ...) {
  # Check if input is either a data frame or part of a model object
  if (is.list(input) && "bootdata" %in% names(input)) {
    # If input is a list and has bootdata, check for ATE column within bootdata
    if (ATE %in% names(input$bootdata)) {
      ate_values <- unlist(input$bootdata[ATE])
    } else {
      stop("bootdata within the model object must have an 'ATE' column.")
    }
  }

  ate_density <- density(ate_values)
  ci <- quantile(ate_values, probs = c(0.025, 0.975))
  density_ci <- density(ate_values, from = ci[1], to = ci[2])

  plot(ate_density, col = col_density, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
  polygon(c(density_ci$x, rev(density_ci$x)), c(rep(min(ate_density$y), length(density_ci$x)), rev(density_ci$y)), col = rgb(0, 0, 1, alpha = 0.3))
  abline(v = mean(ate_values), col = "purple", lwd = 2, lty = 3)
  abline(v = ci[1], col = "darkgreen", lty = 2)
  abline(v = ci[2], col = "darkgreen", lty = 2)

  legend_text <- c(paste(ATE, "Density",sep = " "),
                   paste("Mean:", round(mean(ate_values), 3)),
                   paste("95% CI: [", round(ci[1], 3), ",", round(ci[2], 3), "]"))

  legend("topright", legend = legend_text,
         col = c(col_density, "purple", "darkgreen"),
         lwd = 2, lty = c(1, 3, 2))
}
