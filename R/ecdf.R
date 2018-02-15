#' Empirical CDF estimation for heterogeneity in panel data without bias-correction
#'
#' \code{neecdf} implements the naive estimation of the empirical CDF
#' for the heterogeneous mean, autocovariance, and autocorrelation.
#' The procedure is proposed in Okui and Yanagi (2017).
#' See the package vignette via \code{vignette("panelhetero")} for details.
#'
#' @param data matrix of panel data in which each row is individual time series
#' @param acov_order non-negative integer for the order of the autocovariance
#' @param acor_order positive integer for the order of the autocorrelation
#'
#' @import ggplot2
#' @importFrom stats na.omit
#'
#' @return list that contains the following elements.
#' \item{mean}{graph of the naive empirical CDF estimation for the mean}
#' \item{acov}{graph of the naive empirical CDF estimation for the autocovariance}
#' \item{acor}{graph of the naive empirical CDF estimation for the autocorrelation}
#' \item{mean_func}{function that returns naive empirical CDF estimates for the mean}
#' \item{acov_func}{function that returns naive empirical CDF estimates for the autocovariance}
#' \item{acor_func}{function that returns naive empirical CDF estimates for the autocorrelation}
#' \item{quantity}{matrix that contains the estimated quantities}
#' \item{acov_order}{the same as the argument}
#' \item{acor_order}{the same as the argument}
#' \item{N}{the number of cross-sectional units}
#' \item{S}{the length of time series}
#'
#' @export
#'
neecdf <- function(data, acov_order = 0, acor_order = 1) {

  # initialization
  x <- NULL

  # handling error
  error2(data = data, acov_order = acov_order, acor_order = acor_order)

  # omitting NA
  data <- na.omit(data)

  # sample sizes
  N <- nrow(data)
  S <- ncol(data)

  # estimated means, autocovariances, autocorrelations
  mean_est <- rowMeans(data)
  acov_est <- apply(data, MARGIN = 1, acov, acov_order = acov_order)
  acor_est <- apply(data, MARGIN = 1, acor, acor_order = acor_order)

  # limits for figures by ggplot2
  mean_lim <- c(min(mean_est), max(mean_est))
  acov_lim <- c(min(acov_est), max(acov_est))
  acor_lim <- c(min(acor_est), max(acor_est))

  # figures by ggplot2
  mean_plot <- ggplot(data = data.frame(x = mean_lim), aes(x = x))
  mean_plot <- mean_plot + stat_function(fun = ecdfest, args = list(X = mean_est), geom = "step")
  mean_plot <- mean_plot + labs(x = "x", y = "") + ylim(0, 1)

  acov_plot <- ggplot(data = data.frame(x = acov_lim), aes(x = x))
  acov_plot <- acov_plot + stat_function(fun = ecdfest, args = list(X = acov_est), geom = "step")
  acov_plot <- acov_plot + labs(x = "x", y = "") + ylim(0, 1)

  acor_plot <- ggplot(data = data.frame(x = acor_lim), aes(x = x))
  acor_plot <- acor_plot + stat_function(fun = ecdfest, args = list(X = acor_est), geom = "step")
  acor_plot <- acor_plot + labs(x = "x", y = "") + ylim(0, 1)

  # functions
  mean_func <- function(x) {
    ecdfest(x = x, X = mean_est)
  }

  acov_func <- function(x) {
    ecdfest(x = x, X = acov_est)
  }

  acor_func <- function(x) {
    ecdfest(x = x, X = acor_est)
  }

  # results
  quantity <- cbind(mean_est, acov_est, acor_est)
  colnames(quantity) <- c("mean", "autocovariance", "autocorrelation")
  result <- list(mean = mean_plot, acov = acov_plot, acor = acor_plot,
                 mean_func = mean_func, acov_func = acov_func, acor_func = acor_func,
                 quantity = quantity, acov_order = acov_order,
                 acor_order = acor_order, N = N, S = S)

  return(result)

}


#' computing empirical CDF estimate
#'
#' @param x point at which the CDF is estimated
#' @param X vector of cross-sectional data
#'
ecdfest <- Vectorize(FUN = function(x, X) {

  est <- mean(ifelse(X <= x, 1, 0))

  return(est)

}, vectorize.args = "x")
