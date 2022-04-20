# @description calculate absolute mean standardized differences between vectors X and Y
# @param x \code{numeric} vector, assumed to be "treated" units
# @param y \code{numeric} vector, assumed to be "control" units
# @param v_xy Expected \code{list(var_x, var_y)}, for pre-computed variances (ie. scaling factors)
# to be used in AMSD calculation. The scaling factors should be the same pre and post matching
# for comparability. 
amsd.numeric <- function(x, y, v_xy= NULL, na.rm= TRUE) {
  x_bar <- mean(x, na.rm= na.rm)
  y_bar <- mean(y, na.rm=na.rm)
  if (is.null(v_xy)) {
    v_x <- var(x, na.rm= na.rm)
    v_y <- var(y, na.rm= na.rm)
  } else {
    v_x <- v_xy[[1]]
    v_y <- v_xy[[2]]
  }
  
  return(abs(x_bar - y_bar) / sqrt((v_x + v_y) / 2))
}

# @description calculate absolute mean standardized differences between vectors X and Y.
# with \code{calc_var= TRUE} can be used to calcualte variances vs return AMSD. Useful for 
# pre-calculating the scaling factor \code{v_xy}
# @param DT a class \code{'data.table'} object
# @param v character scalar with the name of the variable to be calculated
# @param t character scalar with the name of the treatment variable
# @param v_xy Expected \code{list(var_x, var_y)}, for pre-computed variances (ie. scaling factors)
# to be used in AMSD calculation. The scaling factors should be the same pre and post matching
# for comparability. 
# @param calc_var 
amsd.multinomial <- function(DT, v, t, v_xy= NULL, na.rm= TRUE, calc_var= FALSE) {
  # 01. grab variables, calculate summary statistics
  tmpDT <- data.table(
    x= get(v, DT)
    , t= get(t, DT)
  )
  sumstats <- tmpDT[, .(
    .N
  ), .(x, t)][ tmpDT[, .N, t], on= "t"][, pct := N / i.N]
  
  # 02. transform to mean and variance
  sumstats <- dcast.data.table(sumstats, formula= x ~ t, value.var= 'pct')
  setnames(sumstats, names(sumstats)[2:3], c('t_0', 't_1'))
  sumstats[, sP := sqrt(1/2 * (t_0 * (1-t_0) + t_1 * (1-t_1)))]
  
  # 03. return amsd
  if (calc_var) {
    out <- sumstats[, .(x, sP)]
    setnames(out, names(out), c(v, 'sPooled'))
    return(out)
  } else{
    if (is.null(v_xy)) {
      sumstats[, amsd := abs(t_1 - t_0) / sP]
    } else {
      setnames(v_xy, v, "x")
      sumstats <- sumstats[v_xy, on= "x"]
      sumstats[, amsd := abs(t_1 - t_0) / sPooled]
      setnames(v_xy, "x", v)
    }
    out <- sumstats[, .(x= paste(v, x, sep= ":"), amsd)]  
    setnames(out, "x", "var")
    return(out)
  }
}
