
# the PDP plot in dbarts sucks. rewrite it in ggplot2

#' @title Partial Dependency Plot
#' @description Create a partial dependence plot using \code{library(ggplot2)} instead of base R's
#' graphics. 
#' @param pdp A class \code{dbarts::pdbart} object. The result of \code{\link[dbarts]{pdbart}}.
#' @param xind Integer scalar indicating which variable to be plotted.
#' @param plquants Numeric vector of length two. In pdp plots, beliefs about \eqn{f(x)} are 
#' indicated by plotting the posterior median, lower, and upper quantiles. \code{plquants} indicates
#' the upper and lower quantiles to be displayed.
#' @param rng_y Numeric vector of length two. Used to limit the range of the Y variable in the plot.
#' @factor Logical. Is the \eqn{X} variable a factor or continuous variable? Defaults to \code{FALSE}.
#' @return A ggplot plot object
#' @seealso \code{\link[dbarts]{pdbart}}.
plot.pdp <- function(pdp, xind, plquants= c(0.05, 0.95)
                     , rng_y= c(0.7, 1), factor= FALSE) {
  
  tsum <- pnorm(apply(pdp$fd[[xind]], 2, quantile, probs= c(plquants[1], 0.5, plquants[2])))
  dt <- data.frame(x= pdp$levs[[xind]], pmin= tsum[1,], p50= tsum[2,], pmax= tsum[3,])
  
  if (is.null(rng_y)) {
    rng_y <- pnorm(range(pdp$fd))
  } 
  if (factor == FALSE) {
    return(
      ggplot(dt, aes(x= x, y= p50)) + geom_point() + geom_line() +
        geom_ribbon(aes(ymin= pmin, ymax= pmax), fill= 'grey70', alpha= 0.5) +
        labs(x= pdp$xlbs[xind], y= 'partial dependence'
             , title= paste('Partial Dependence:', pdp$xlbs[xind])) +
        scale_y_continuous(breaks= round(seq(rng_y[1], rng_y[2], length.out= 10), 3)
                           , limits= rng_y
                           , labels= scales::percent) +
        theme(plot.title= element_text(size= 14, face= 'bold')
              , axis.title= element_text(size= 12, face= 'bold')
              , axis.text= element_text(size= 11))
    )  
  } else {
    return(
      ggplot(dt, aes(x= factor(x), y= p50)) + 
        geom_crossbar(aes(ymin= pmin, ymax= pmax), fill= 'grey70', alpha= 0.5) +
        labs(x= pdp$xlbs[xind], y= 'partial dependence'
             , title= paste('Variable Importance:', pdp$xlbs[xind])) +
        scale_y_continuous(breaks= round(seq(rng_y[1], rng_y[2], length.out= 10), 3)
                           , limits= rng_y
                           , labels= scales::percent) +
        theme(plot.title= element_text(size= 14, face= 'bold')
              , axis.title= element_text(size= 12, face= 'bold')
              , axis.text= element_text(size= 11))
    )
  }
  
}
