#' computing the autocovariance of time series data
#'
#' @param y vector of time series data
#' @param acov_order non-negative integer for the order of the autocovariance
#'
acov <- function(y, acov_order) {

  # the order
  k <- acov_order

  # the length of time series
  S <- length(y)

  # the mean
  mean_est <- mean(y)

  # the autocovariance
  y1 <- y[1:(S - k)] - mean_est
  y2 <- y[(k + 1):S] - mean_est
  acov_est <- sum(y1 * y2) / (S - k)

  return(acov_est)

}

#' computing the autocorrelation of time series data
#'
#' @param y vector of time series data
#' @param acor_order positive integer for the order of the autocorrelation
#'
acor <- function(y, acor_order) {

  # the autocorrelation
  acor_est <- acov(y = y, acov_order = acor_order) / acov(y = y, acov_order = 0)

  return(acor_est)

}
