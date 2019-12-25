#' function 1 to handle errors
#'
#' @param data matrix of panel data in which each row is individual time series
#' @param acov_order non-negative integer for the order of the autocovariance
#' @param acor_order positive integer for the order of the autocorrelation
#' @param R positive integer for the number of bootstrap replications
#'
error1 <- function(data, acov_order, acor_order, R) {

  # error for data
  if (!is.matrix(data)) {
    stop(message = "data must be a matrix.")
  }

  # error for acov_order
  if ((!is.numeric(acov_order)) || (length(acov_order) > 1) || (acov_order < 0)) {
    stop(message = "acov_order must be a non-negative integer.")
  }

  # error for acor_order
  if ((!is.numeric(acor_order)) || (length(acor_order) > 1) || (acor_order < 1)) {
    stop(message = "acor_order must be a positive integer.")
  }

  # error for R
  if ((!is.numeric(R)) || (length(R) > 1) || (R < 1)) {
    stop(message = "R must be a positive integer.")
  }

}


#' function 2 to handle errors
#'
#' @param data matrix of panel data in which each row is individual time series
#' @param acov_order non-negative integer for the order of the autocovariance
#' @param acor_order positive integer for the order of the autocorrelation
#'
error2 <- function(data, acov_order, acor_order) {

  # error for data
  if (!is.matrix(data)) {
    stop(message = "data must be a matrix.")
  }

  # error for acov_order
  if ((!is.numeric(acov_order)) || (length(acov_order) > 1) || (acov_order < 0)) {
    stop(message = "acov_order must be a non-negative integer.")
  }

  # error for acor_order
  if ((!is.numeric(acor_order)) || (length(acor_order) > 1) || (acor_order < 1)) {
    stop(message = "acor_order must be a positive integer.")
  }

}


#' function 3 to handle errors
#'
#' @param data matrix of panel data in which each row is individual time series
#' @param acov_order non-negative integer for the order of the autocovariance
#' @param acor_order positive integer for the order of the autocorrelation
#' @param mean_bw bandwidth for the mean
#' @param acov_bw bandwidth for the autocovariance
#' @param acor_bw bandwidth for the autocorrelation
#'
error3 <- function(data, acov_order, acor_order, mean_bw, acov_bw, acor_bw) {

  # error for data
  if (!is.matrix(data)) {
    stop(message = "data must be a matrix.")
  }

  # error for acov_order
  if ((!is.numeric(acov_order)) || (length(acov_order) > 1) || (acov_order < 0)) {
    stop(message = "acov_order must be a non-negative integer.")
  }

  # error for acor_order
  if ((!is.numeric(acor_order)) || (length(acor_order) > 1) || (acor_order < 1)) {
    stop(message = "acor_order must be a positive integer.")
  }

  # error for mean_bw
  if ((!is.null(mean_bw)) && ((!is.numeric(mean_bw)) || (length(mean_bw) > 1) || (mean_bw <= 0))) {
    stop(message = "mean_bw must be a positive bandwidth.")
  }

  # error for acov_bw
  if ((!is.null(acov_bw)) && ((!is.numeric(acov_bw)) || (length(acov_bw) > 1) || (acov_bw <= 0))) {
    stop(message = "acov_bw must be a positive bandwidth.")
  }

  # error for acor_bw
  if ((!is.null(acor_bw)) && ((!is.numeric(acor_bw)) || (length(acor_bw) > 1) || (acor_bw <= 0))) {
    stop(message = "acor_bw must be a positive bandwidth.")
  }

}