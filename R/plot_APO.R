#' Plot Average Potential Outcomes (APO)
#'
#' This function plots the density of APO for a specified effect type from bootstrap simulation results.
#'
#' @param input A data frame or model object containing 'bootdata', which include 'effect_comparator' and 'effect_reference' columns.
#' @param effect_type A character string specifying which effect to plot: 'effect_comparator' or 'effect_reference'.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A density plot showing the distribution of the specified average potential outcome (reference or comparison).
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
#'                            parallel = FALSE,
#'                            ncore = 2)
#' plot_APO(model$bootdata, effect_type = "effect_comparator")
#' plot_APO(model, effect_type = "effect_reference")
#'
plot_APO <- function(input, effect_type, ...) {
  # Validate input
  if ("bootdata" %in% names(input)) {
    bootdata <- input$bootdata
  } else if (is.data.frame(input)) {
    bootdata <- input
  } else {
    stop("Input must be a data frame or a model object containing a 'bootdata' data frame.")
  }
  if (!is.data.frame(bootdata) || !("effect_comparator" %in% names(bootdata)) || !("effect_reference" %in% names(bootdata))) {
    stop("bootdata must be a data frame containing 'effect_comparator' and 'effect_reference' columns.")
  }
  if (!is.character(effect_type) || length(effect_type) != 1) {
    stop("effect_type must be a single character string specifying the effect to plot.")
  }
  if (!effect_type %in% c("effect_comparator", "effect_reference")) {
    stop("effect_type must be either 'effect_comparator' or 'effect_reference'.")
  }

  # Extract the relevant column
  effect <- bootdata[, effect_type, drop = FALSE]

  # Calculate density
  density_effect <- density(effect[[1]])

  # Define titles and colors based on effect_type
  titles <- c(effect_comparator = "Comparator Level", effect_reference = "Reference Level")
  colors <- c(effect_comparator = "blue", effect_reference = "red")

  # Calculate mean
  mean_effect <- mean(effect[[1]])

  # Calculate CI
  ci <- quantile(effect[[1]], probs = c(0.025, 0.975))
  density_ci <- density(effect[[1]], from = ci[1], to = ci[2])

  # Plotting
  plot(density_effect, main = paste("Average Potential Outcome (APO) of", titles[effect_type]), xlab = "Effect", ylab = "Density", col = colors[effect_type], lwd = 2, ...)

  # Shade the area under the curve within the 95% CI
  polygon(c(density_ci$x, rev(density_ci$x)), c(rep(min(density_effect$y), length(density_ci$x)), rev(density_ci$y)), col = rgb(0, 0, 1, alpha = 0.3))

  # Add vertical lines for the mean and 95% CI bounds
  abline(v = mean_effect, col = "darkgrey", lty = 3)
  abline(v = ci[1], col = "darkgreen", lty = 2)
  abline(v = ci[2], col = "darkgreen", lty = 2)

  # Legend with mean and 95% CI bounds
  legend_text <- c(titles[effect_type],
                   paste("Mean:", round(mean_effect, 3)),
                   paste("95% CI: [", round(ci[1], 3), ",", round(ci[2], 3), "]"))

  legend("topright", legend = legend_text,
         col = c(colors[effect_type], "purple", "darkgreen"),
         lwd = 2, lty = c(1, 3, 12))
}
