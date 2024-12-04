#' A function to calculate the effect of an intervention given the parameter estimates and intervention levels
#'
#' @param intervention_levels A numeric vector indicating the levels of intervention for each predictor variable.
#' @param variables A list of the names of the response variable and predictor variables extracted from the model.
#' @param param_estimates A vector of parameter estimates from the model.
#' @param treatment_effect_type Character string specifying the type of treatment effect to estimate. Options are "sq" for sequential treatment effects, which estimates effects for specific treatment sequences across visits, and "cum" for cumulative treatment effects, which assumes a single cumulative treatment variable representing the total exposure. The default is "sq".
#'
#' @return A numeric value representing the calculated effect of the specified intervention.
#'
calculate_effect <- function(intervention_levels, variables, param_estimates, treatment_effect_type) {
  if (treatment_effect_type == "cum") {
    # For cumulative treatment, only consider b1
    effect <- param_estimates[1] * intervention_levels[1]
  } else {
    # Start with the intercept term
    effect<-effect_intercept<-param_estimates[1]

    # Go through each predictor and add its contribution
    for (i in 1:length(variables$predictors)) {
      term <- variables$predictors[i]
      term_variables <- unlist(strsplit(term, ":"))
      term_index <- which(names(param_estimates) == term)

      # Calculate the product of intervention levels for the interaction term
      term_contribution <- param_estimates[term_index]
      for (term_variable in term_variables) {
        var_index <- which(variables$predictors == term_variable)
        term_contribution <- term_contribution * intervention_levels[var_index]
      }

      # Add the term contribution to the effect
      effect <- effect + term_contribution
    }
  }
  return(effect)
}
