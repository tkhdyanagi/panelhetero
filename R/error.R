#' Error handling 1
#'
#' @param data A matrix of panel data
#' @param acov_order A non-negative integer of the order of autocovariance
#' @param acor_order A positive integer of the order of autocorrelation
#' @param R A positive integer of the number of bootstrap repetitions
#'
#' @returns NULL
#'
error1 <- function(data, acov_order, acor_order, R) {

  # Error for data
  if (!is.matrix(data)) {
    stop(message = "data must be a matrix.")
  }

  # Error for acov_order
  if (!is.numeric(acov_order) | length(acov_order) > 1 | acov_order < 0) {
    stop(message = "acov_order must be a non-negative integer.")
  }

  # Error for acor_order
  if (!is.numeric(acor_order) | length(acor_order) > 1 | acor_order < 1) {
    stop(message = "acor_order must be a positive integer.")
  }

  # Error for R
  if (!is.numeric(R) | length(R) > 1 | R < 1) {
    stop(message = "R must be a positive integer.")
  }

}

#' Error handling 2
#'
#' @param data A matrix of panel data
#' @param acov_order A non-negative integer of the order of autocovariance
#' @param acor_order A positive integer of the order of autocorrelation
#'
#' @returns NULL
#'
error2 <- function(data, acov_order, acor_order) {

  # Error for data
  if (!is.matrix(data)) {
    stop(message = "data must be a matrix.")
  }

  # Error for acov_order
  if (!is.numeric(acov_order) | length(acov_order) > 1 | acov_order < 0) {
    stop(message = "acov_order must be a non-negative integer.")
  }

  # Error for acor_order
  if (!is.numeric(acor_order) | length(acor_order) > 1 | acor_order < 1) {
    stop(message = "acor_order must be a positive integer.")
  }

}

#' Error handling 3
#'
#' @param data A matrix of panel data
#' @param acov_order A non-negative integer of the order of autocovariance
#' @param acor_order A positive integer of the order of autocorrelation
#' @param mean_bw A scalar of bandwidth used for the estimation of
#' the denisty of mean
#' @param acov_bw A scalar of bandwidth used for the estimation of
#' the denisty of autocovariance
#' @param acor_bw A scalar of bandwidth used for the estimation of
#' the denisty of autocorrelation
#'
#' @returns NULL
#'
error3 <- function(data, acov_order, acor_order, mean_bw, acov_bw, acor_bw) {

  # Error for data
  if (!is.matrix(data)) {
    stop(message = "data must be a matrix.")
  }

  # Error for acov_order
  if (!is.numeric(acov_order) | length(acov_order) > 1 | acov_order < 0) {
    stop(message = "acov_order must be a non-negative integer.")
  }

  # error for acor_order
  if (!is.numeric(acor_order) | length(acor_order) > 1 | acor_order < 1) {
    stop(message = "acor_order must be a positive integer.")
  }

  # Error for mean_bw
  if (!is.null(mean_bw)) {
    if (!is.numeric(mean_bw) | length(mean_bw) > 1 | mean_bw <= 0) {
      stop(message = "mean_bw must be a positive bandwidth.")
    }
  }

  # error for acov_bw
  if (!is.null(acov_bw)) {
    if (!is.numeric(acov_bw) | length(acov_bw) > 1 | acov_bw <= 0) {
      stop(message = "acov_bw must be a positive bandwidth.")
    }
  }

  # error for acor_bw
  if (!is.null(acor_bw)) {
    if (!is.numeric(acor_bw) | length(acor_bw) > 1 | acor_bw <= 0) {
      stop(message = "acor_bw must be a positive bandwidth.")
    }
  }

}