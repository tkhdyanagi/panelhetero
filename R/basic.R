#' Compute the autocovariance of time series
#'
#' @param y A vector of time series
#' @param acov_order A non-negative integer of the order of autocovariance
#'
#' @returns A scalar of the autocovariance
#'
acov <- function(y, acov_order) {

  # Variable definitions
  k <- acov_order
  S <- length(y)
  mean_est <- mean(y)

  # Autocovariance
  y1 <- y[1:(S - k)] - mean_est
  y2 <- y[(k + 1):S] - mean_est
  acov_est <- sum(y1 * y2) / (S - k)

  # Return
  return(acov_est)

}

#' Compute the autocorrelation of time series
#'
#' @param y A vector of time series
#' @param acor_order A positive integer of the order of autocorrelation
#'
#' @returns A scalar of the autocorrelation
#'
acor <- function(y, acor_order) {
  acov(y = y, acov_order = acor_order) / acov(y = y, acov_order = 0)
}
