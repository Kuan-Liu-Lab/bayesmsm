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
plot_APO <- function(model, effect_type,
                     col_density = "blue",
                     fill_density = "lightblue",
                     main = "Average Potential Outcome (APO)",
                     xlab = "Effect", ylab = "Density",
                     xlim = NULL, ylim = NULL, ...) {

  # Validate input
  if ("bootdata" %in% names(model)) {
    bootdata <- model$bootdata
  } else if (is.data.frame(model)) {
    bootdata <- model
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
  effect <- bootdata[[effect_type]]

  # Calculate density
  density_effect <- density(effect)

  # Define titles and colors based on effect_type
  titles <- c(effect_comparator = "Comparator Level", effect_reference = "Reference Level")
  colors <- c(effect_comparator = "blue", effect_reference = "red")

  # Calculate mean and CI
  mean_effect <- mean(effect)
  ci <- quantile(effect, probs = c(0.025, 0.975))
  density_ci <- density(effect, from = ci[1], to = ci[2])

  # Set layout to allocate space for legend
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))

  # Plot the density
  par(mar = c(5, 4, 4, 0)) # Adjust margins for the plot
  plot(density_effect, col = col_density, main = paste(main, titles[effect_type], sep = ": "), xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
  polygon(c(density_ci$x, rev(density_ci$x)), c(rep(min(density_effect$y), length(density_ci$x)), rev(density_ci$y)), col = rgb(0, 0, 1, alpha = 0.3))
  abline(v = mean_effect, col = "purple", lwd = 2, lty = 3)
  abline(v = ci[1], col = "darkgreen", lty = 2)
  abline(v = ci[2], col = "darkgreen", lty = 2)

  # Plot the legend in a separate plot
  par(mar = c(5, 0, 4, 0)) # Adjust margins for the legend
  plot.new()
  legend("center", legend = c(titles[effect_type],
                              paste("Mean:", round(mean_effect, 3)),
                              paste("95% CI: [", round(ci[1], 3), ",", round(ci[2], 3), "]")),
         col = c(col_density, "purple", "darkgreen"),
         lwd = 2, lty = c(1, 3, 2),
         bty = "n")
}
