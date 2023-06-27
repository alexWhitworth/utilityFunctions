
#' @title Unequal Probability Sampling. Estimate mean
#' @description Mean estimators for unequal probability samples. See details. 
#' @param y A numeric vector of outcomes for all units
#' @param pi A numeric vector  of inclusion probabilities for each unit.
#' @param trim Logical. Should extreme weights be trimmed?
#' @param trim_max Numeric scalar for the max relative weight to trim. Defaults to \code{100}
#' @param method One of \code{c('Horvitz-Thomson', 'Generalized-unequal-prob', 'Hansen-Hurwitz', 'weighted-mean')}.
#' See Details.
#' @section Details:
#' Multiple methods are implemented. Hansen-Hurwitz is the only provided sampling with replacement 
#' estimator. Weighted mean is also provided but produces a biased estimator. In general, 
#' 'Hansen-Hurwitz' and 'Horvitz-Thompson' are preferred if \code{y} is approximately proportional
#' to \code{pi}; 'Generalized-unequal-prob' is preferred if \code{y} and \code{pi} are not well
#' related. No estimator will be uniformly better than any other
#' @reference Thompson, S.K. Sampling, (3rd Ed). Chapter 6. (2012)
mean_est <- function(y, pi, na.rm= TRUE, trim= TRUE, trim_max= 100, 
                           method= c('Horvitz-Thomson', 'Generalized-unequal-prob', 'Hansen-Hurwitz'
                                     , 'weighted-mean')) {
  method <- match.arg(method)
  pi <- pi / sum(pi)  # ensure unit norm
  if (trim) {
    pi <- ifelse(pi < 1 / trim_max, 1 / trim_max, pi)
  }
  N <- sum(1 / pi, na.rm= na.rm)
  
  mu <- switch(method
          , 'Hansen-Hurwitz' = 1 / N * (1 / len(y)) * sum(y / pi, na.rm= na.rm) # sample w/ replacement
          , 'Horvitz-Thompson' = 1 / N * sum(y / pi, na.rm= na.rm) # sample w/o replacement
          , 'Generalized-unequal-prob' = sum(y/pi, na.rm= na.rm)
          , 'weighted-mean' = weighetd.mean(y, w= 1/pi, na.rm= na.rm)
  )
  return(mu)
}

#' @title Unequal Probability Sampling. Estimate variance
#' @description Mean estimators for unequal probability samples. See details. 
#' @param y A numeric vector of outcomes for all units
#' @param pi A numeric vector  of inclusion probabilities for each unit.
#' @param mu A numeric scalar for the mean estimate. See \code{\link{mean_est}}.
#' @param trim Logical. Should extreme weights be trimmed?
#' @param trim_max Numeric scalar for the max relative weight to trim. Defaults to \code{100}
#' @param method One of \code{c('Horvitz-Thomson', 'Generalized-unequal-prob', 'Hansen-Hurwitz', 'weighted-mean')}.
#' See Details.
#' @section Details:
#' Multiple methods are implemented. Hansen-Hurwitz is the only provided sampling with replacement 
#' estimator. Weighted mean is also provided but produces a biased estimator. In general, 
#' 'Hansen-Hurwitz' and 'Horvitz-Thompson' are preferred if \code{y} is approximately proportional
#' to \code{pi}; 'Generalized-unequal-prob' is preferred if \code{y} and \code{pi} are not well
#' related. No estimator will be uniformly better than any other
#' 
#' Note: We assume independence for all units \eqn{i, j, \in \{1, \ldots, N \}} if \eqn{i \neq j}.
#' @reference Thompson, S.K. Sampling, (3rd Ed). Chapter 6. (2012)
var_est <- function(y, pi, mu, na.rm= TRUE, trim= TRUE, trim_max= 100, 
                    method= c('Horvitz-Thomson', 'Generalized-unequal-prob', 'Hansen-Hurwitz'
                              , 'weighted-mean')) {
  
  method <- match.arg(method)
  pi <- pi / sum(pi)  # ensure unit norm
  if (trim) {
    pi <- ifelse(pi < 1 / trim_max, 1 / trim_max, pi)
  }
  N <- sum(1 / pi, na.rm= na.rm)
  
  # population total
  tau <- switch(method
          , 'Hansen-Hurwitz' = (1 / len(y)) * sum(y / pi, na.rm= na.rm)
          , 'Horvitz-Thompson' = sum(y / pi, na.rm= na.rm)
          , 'Generalized-unequal-prob' = mu * N
  )
  
  # assumed independence of units: i and j, i != j
  v <- switch(method
          , 'Hansen-Hurwitz' = 1 / (N * (N - 1)) * sum( (y / pi - tau)^2 )
          , 'Horvitz-Thompson' = sum( (1 - pi) / pi *y^2)
          , 'Generalized-unequal-prob' = sum( (1 - pi) / pi * (y - mu)^2 / pi, na.rm= na.rm)
          , 'weighted-mean' = Hmisc::wtd.var(y, w= 1/pi, na.rm= na.rm) 
  )
  return(v / N^2) # variance of mean
}
  

#' @title Unequal Probability Sampling Estimators
#' @description A function calculating unequal probability sampling estimators (mean and variance).
#' @param y A numeric vector of outcomes for all units
#' @param pi A numeric vector  of inclusion probabilities for each unit.
#' @param trim Logical. Should extreme weights be trimmed?
#' @param trim_max Numeric scalar for the max relative weight to trim. Defaults to \code{100}
#' @param method One of \code{c('Horvitz-Thomson', 'Generalized-unequal-prob', 'Hansen-Hurwitz', 'weighted-mean')}.
#' See Details.
#' @section Details:
#' Multiple methods are implemented. Hansen-Hurwitz is the only provided sampling with replacement 
#' estimator. Weighted mean is also provided but produces a biased estimator. In general, 
#' 'Hansen-Hurwitz' and 'Horvitz-Thompson' are preferred if \code{y} is approximately proportional
#' to \code{pi}; 'Generalized-unequal-prob' is preferred if \code{y} and \code{pi} are not well
#' related. No estimator will be uniformly better than any other
#' @reference Thompson, S.K. Sampling, (3rd Ed). Chapter 6. (2012)
prevalence_est <- function(y, pi, na.rm= TRUE, trim= TRUE, trim_max= 100, 
                           method= c('Horvitz-Thomson', 'Generalized-unequal-prob', 'Hansen-Hurwitz'
                                     , 'weighted-mean')) {
  method <- match.arg(method)
  hat_mu <- mean_est(y, pi, trim= trim, trim_max= trim_max, method= method)
  hat_var <- var_est(y, pi, mu= hat_mu, trim= trim, trim_max= trim_max, method= method)
  return(data.table::data.table(hat_mu= hat_mu, hat_v= har_var))
}
