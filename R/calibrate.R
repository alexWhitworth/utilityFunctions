
# calibrate model predictions 2-dimensionally

#' @title Calibrate model under/over fit
#' @description Calibrate model predictions conditional on a two-dimensional grid of predictors. 
#' Specifically, check if the model is over (or under) predicting vs true rate. 
#' \eqn{Cal_{j,k} = frac{\sum_{i=1}^N \hat Y_i}{\sum_{i=1}^N Y_i} | X_j = x_j, X_k = x_k}
#' @param DT A \code{data.table} object with your data
#' @param var1 A character scalar for the name of \eqn{X_1} in \code{DT} to condition on.
#' @param var2 A character scalar for the name of \eqn{X_2} in \code{DT}  to condition on.
#' @param y A character scalar for the name of \eqn{Y} in \code{DT}
#' @param yhat A character scalar for the name of \eqn{\hat Y} in \code{DT}
#' @return A class \code{data.table} object 
calibrate <- function(DT, var1, var2, y, yhat= 'hat_y') {
  v1 <- get(var1, as.environment(DT))
  v2 <- get(var2, as.environment(DT))
  
  # if v1, v2 are numeric/continuous, bin into deciles
  if (is.numeric(v1)) {
    v1 <- gtools::quantcut(v1, q= 10)
    levels(v1) <- paste0("p", seq(0, 90, 10), "_p", seq(10, 100, 10))
  }
  if (is.numeric(v2)) {
    v2 <- gtools::quantcut(v2, q= 10)
    levels(v2) <- paste0("p", seq(0, 90, 10), "_p", seq(10, 100, 10))
  }
  # form tmp DT with necessary inputs for calculation
  DT <- data.table::data.table(
    var1= var1, lvls1= v1
    , var2= var2, lvls2= v2
    , y= get(y, as.environment(DT))
    , yhat= get(yhat, as.environment(DT))
  )
  return(DT[, .(.N, cal= sum(yhat) / sum(y)), .(var1,lvls1,var2, lvls2)][order(lvls1, lvls2)])
}