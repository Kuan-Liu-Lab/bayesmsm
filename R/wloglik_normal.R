#' Weighted Log-Likelihood for Simple Linear Regression
#'
#' This function calculates the weighted log-likelihood for a simple linear
#' regression model where the outcome y is continuous, and we assume the residual
#' follows a normal distribution.
#'
#' @param param Numeric vector containing parameters (coefficients for treatment and standard deviation).
#' @param Y Numeric vector of dependent variable (outcome).
#' @param A Matrix of independent variables (treatment variables).
#' @param weight Numeric vector of treatment weights for each observation.
#' @return Numeric value of the weighted log-likelihood.
#' @export
#'
wloglik_normal<-function(param,
                         Y,
                         A,
                         weight){
  #number of observations;
  n <- length(Y)
  theta <- param[1:dim(A)[2]] #causal parameters on the mean
  #number of parameter is determined by number of treatment variables, plus intercept;
  sigma <- param[(dim(A)[2]+1)] # the remaining the parameter represent the standard deviation;
  mmat <- as.matrix(A) #design matrix of the causal outcome model, e.g., A = cbind(1, a_1, a_2);
  logl<- -0.5*log(sigma^2) - 0.5*((Y - mmat%*%theta)^2)/(sigma^2)
  wlogl<-sum(weight*logl)

  return(wlogl)
}
