#' Error bar plots for treatment effects
#'
#' @param input A data frame or model object containing 'bootdata', which include 'effect_comparator', 'effect_reference', and 'RD' columns.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return An error bar plot of the mean effects and their 95\% confidence intervals for comparator level, reference level, and ATE.
#' @export
#'
#' @examples
#' testdata <- read.csv(system.file("extdata",
#'                                  "continuous_outcome_data.csv",
#'                                  package = "bayesmsm"))
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
#' plot_est_box(model$bootdata) # without reference & comparator information
#'                              # below labels
#' plot_est_box(model) # with reference & comparator information below labels
#'
plot_est_box <- function(input, ...) {
  # Extract bootdata from the model or use the data frame directly
  bootdata <- if (is.data.frame(input)) {
    input
  } else if ("bootdata" %in% names(input)) {
    input$bootdata
  } else {
    stop("Input must be a data frame or a model object containing 'bootdata'.")
  }

  # Identify if the family is binomial to include RR and OR
  is_binomial <- "RR" %in% names(bootdata) && "OR" %in% names(bootdata)

  # Validate bootdata
  required_columns <- c("effect_comparator", "effect_reference", "RD")
  if (is_binomial) {
    required_columns <- c(required_columns, "RR", "OR")
  }
  if (!all(required_columns %in% names(bootdata))) {
    stop("bootdata must contain 'effect_comparator', 'effect_reference', and 'RD' columns.")
  }

  # Calculate means and standard deviations
  means <- sapply(bootdata[required_columns], mean)
  lowerbd <- sapply(bootdata[required_columns], function(x) quantile(x, probs = 0.025))
  upperbd <- sapply(bootdata[required_columns], function(x) quantile(x, probs = 0.975))

  # Define the position for each point
  position <- 1:length(means)

  # Define some offsets for text placement
  text_offset <- (max(upperbd) - min(lowerbd)) * 0.05

  old_par <- par(no.readonly = TRUE) # Save current graphical parameters
  on.exit(par(old_par)) # Restore graphical parameters on exit
  par(mar = c(5, 4, 4, 3) + 0.1) # Adjust margins

  # Plotting
  plot(position, means, ylim = range(lowerbd - text_offset, upperbd + text_offset), pch = 19, xaxt = "n", # round down vs round up;
       xlab = "Treatment Level", ylab = "Effect", main = "Treatment Effect Estimates", ...)
  axis(1, at = position, labels = if (is_binomial) c("Comparator Level", "Reference Level", "RD", "RR", "OR") else c("Comparator Level", "Reference Level", "RD"))

  # Error bars
  arrows(position, lowerbd, position, upperbd, angle = 90, code = 3, length = 0.1, ...)

  # Adding text for means and CIs
  for (i in seq_along(means)) {
    text(position[i], upperbd[i] + text_offset, labels = paste("Mean:", round(means[i], 2)), cex = 0.8, pos = 3)
    text(position[i], upperbd[i] + 2 * text_offset, labels = paste("95% CI: [", round(lowerbd[i], 2), ", ", round(upperbd[i], 2), "]"), cex = 0.8, pos = 3)
  }

  # Check if the input is a model and extract treatment sequences if they exist
  has_treatment_info <- "reference" %in% names(input) && "comparator" %in% names(input)

  # Conditional treatment sequence information below x-axis labels
  if (has_treatment_info) {
    mtext(paste("(", paste(input$reference, collapse = ", "), ")", sep = ""), side = 1, at = position[2], line = 2)
    mtext(paste("(", paste(input$comparator, collapse = ", "), ")", sep = ""), side = 1, at = position[1], line = 2)
  }
}
