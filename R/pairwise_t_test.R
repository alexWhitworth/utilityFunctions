#' @title Pairwise T Test
#' @description Pairwise T-test for independent two-sample tests  with assumed approximately
#' equal variances. Also allows Welch-Satterthwaite DF. Can also be used for single T-tests
#' vs multiple pairwise comparisons
#' @param y_bar A \code{numeric} vector of sample means
#' @param sigma A \code{numeric} vector of associated sample standard deviations
#' @param y_bar A \code{numeric} vector of sample sizes
#' @param p.adjust.method A \code{character} string specifying the correcton method. See
#' \code{\link[stats]{p.adjust}}.
#' @param welch Logical. Do you want to calculate a Welch T-test with Welch-Satterthwaite 
#' degrees of freedom?
#' @return an object of class 'pairwise.htest'
#' @examples {
#'     ## single t-test
#'     pairwise_t_test(y_bar= rnorm(2), sigma= rnorm(2, mean=1)
#'                     , n= sample.int(100, size= 2), p.adjust.method= 'none')
#'     ## pairwise tests
#'     Y <- rnorm(10); sig <- rnorm(10, mean= 1); n <- sample.int(100, size= 10)
#'     pairwise_t_test(y_bar= Y, sigma= sig, n= n, p.adjust.method= 'holm')
#' }
pairwise_t_test <- function(y_bar, sigma, n, p.adjust.method= c("holm", "hochberg", "hommel", 
                                                                "bonferroni", "BH", "BY",
                                                                "fdr", "none")
                            , welch= FALSE) {
  if (length(y_bar) != length(sigma)) {
    stop('y_bar, sigma, and n must have equal length.')
  }
  if (length(y_bar) != length(n)) {
    stop('y_bar, sigma, and n must have equal length.')
  }
  
  # always two sided
  p.adjust.method <- match.arg(p.adjust.method)
  N <- length(y_bar)
  t_mat <- p_mat <- rel_diff <- matrix(NA, nrow= N-1, ncol= N-1)
  # 01. calculate T-stat + P-stat
  for (i in 1:(N-1)) {
    for (j in (i+1):(N)) {
      if (welch) {
        sp <- sqrt(sigma[i]^2 / n[i] + sigma[j]^2 / n[j])
        tstat <- (y_bar[i] - y_bar[j]) / sp
        df <- df_welch_sat(sd= c(sigma[i], sigma[j]), n= c(n[i], n[j]))
      } else {
        sp <- sqrt(
          ((n[i] - 1) * sigma[i]^2 + (n[j] - 1) * sigma[j]^2) / (n[i] + n[j] - 2)
        )
        tstat <- (y_bar[i] - y_bar[j]) / (sp * sqrt(1 / n[i] + 1 / n[j]))
        df <- n[i] + n[j] - 2
      }
      t_mat[i, (j-1)] <- tstat
      rel_diff[i, (j-1)] <- (y_bar[i] - y_bar[j]) / y_bar[j]
      if (tstat >= 0) {
        p_mat[i, (j-1)] <- (1 - pt(tstat, df= df)) + pt(-1 * tstat, df= df)
      } else {
        p_mat[i, (j-1)] <- pt(tstat, df= df) + (1 - pt(-1 * tstat, df= df))
      }
    }
  }
  # 02. control family-wise error rate
  dd <- dim(p_mat)
  p_mat <- stats::p.adjust(p= p_mat, method= p.adjust.method, n= sum(!is.na(p_mat)))
  dim(p_mat) <- dd
  # 03. return
  ret <- list(rel_diff= rel_diff, method= NULL, data.name= NULL
              , p.value= p_mat, p.adjust.method= p.adjust.method)
  class(ret) <- 'pairwise.htest'
  return(ret)
}

# degrees of freedom via Welch-Satterthwaite equation
df_welch_sat <- function(sd, n) {
  return(
    (sd[1]^2 / n[1] + sd[2]^2 / n[2])^2 / 
      ( (sd[1]^2 / n[1])^2 / (n[1] - 1) + (sd[2]^2 / n[2])^2 / (n[2] - 1) )
  )
}



print.pairwise.htest <- function (x, digits = max(1L, getOption("digits") - 4L), ...) 
{
  cat("\n\tPairwise comparisons using", x$method, "\n\n")
  cat("relative difference: ", "\n")
  rd <- format(x$rel_diff, digits= digits)
  print(rd, quote= FALSE, ...)
  cat("data: ", x$data.name, "\n\n")
  pp <- format.pval(x$p.value, digits = digits, na.form = '-')
  attributes(pp) <- attributes(x$p.value) 
  print(pp, quote = FALSE, ...)
  cat("\nP value adjustment method:", x$p.adjust.method, "\n")
  invisible(x)
}
