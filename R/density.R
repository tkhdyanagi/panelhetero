#' Kernel density estimation for heterogeneity in panel data without bias-correction
#'
#' \code{nekd} implements the naive estimation of the kernel density
#' for the heterogeneous mean, autocovariance, and autocorrelation.
#' The procedure is proposed in Okui and Yanagi (2018).
#' See the package vignette via \code{vignette("panelhetero")} for details.
#'
#' @param data matrix of panel data in which each row is individual time series
#' @param acov_order non-negative integer for the order of the autocovariance
#' @param acor_order positive integer for the order of the autocorrelation
#' @param mean_bw bandwidth for the mean
#' @param acov_bw bandwidth for the autocovariance
#' @param acor_bw bandwidth for the autocorrelation
#'
#' @import ggplot2
#' @importFrom KernSmooth dpik
#' @importFrom stats na.omit
#'
#' @return list that contains the following elements.
#' \item{mean}{graph of the naive kernel density estimation for the mean}
#' \item{acov}{graph of the naive kernel density estimation for the autocovariance}
#' \item{acor}{graph of the naive kernel density estimation for the autocorrelation}
#' \item{mean_func}{function that returns naive kernel density estimates for the mean}
#' \item{acov_func}{function that returns naive kernel density estimates for the autocovariance}
#' \item{acor_func}{function that returns naive kernel density estimates for the autocorrelation}
#' \item{bandwidth}{vector of the selected bandwidths}
#' \item{quantity}{matrix that contains the estimated quantities}
#' \item{acov_order}{the same as the argument}
#' \item{acor_order}{the same as the argument}
#' \item{N}{the number of cross-sectional units}
#' \item{S}{the length of time series}
#'
#' @export
#'
nekd <- function(data, acov_order = 0, acor_order = 1, mean_bw = NULL, acov_bw = NULL, acor_bw = NULL) {

  # initialization
  x <- NULL

  # handling errors
  error3(data = data, acov_order = acov_order, acor_order = acor_order, mean_bw = mean_bw, acov_bw = acov_bw, acor_bw = acor_bw)

  # omitting NA
  data <- na.omit(data)

  # sample sizes
  N <- nrow(data)
  S <- ncol(data)

  # estimated means, autocovariances, autocorrelations
  mean_est <- rowMeans(data)
  acov_est <- apply(data, MARGIN = 1, acov, acov_order = acov_order)
  acor_est <- apply(data, MARGIN = 1, acor, acor_order = acor_order)

  # plug-in bandwidths
  if (is.null(mean_bw)) {
    mean_bw <- dpik(x = mean_est, scalest = "minim", kernel = "normal")
  }
  if (is.null(acov_bw)) {
    acov_bw <- dpik(x = acov_est, scalest = "minim", kernel = "normal")
  }
  if (is.null(acor_bw)) {
    acor_bw <- dpik(x = acor_est, scalest = "minim", kernel = "normal")
  }

  # limit for figures by ggplot2
  mean_lim <- c(min(mean_est), max(mean_est))
  acov_lim <- c(min(acov_est), max(acov_est))
  acor_lim <- c(min(acor_est), max(acor_est))

  # figures by ggplot2
  mean_plot <- ggplot(data = data.frame(x = mean_lim), aes(x = x))
  mean_plot <- mean_plot + stat_function(fun = kdest, args = list(X = mean_est, h = mean_bw))
  mean_plot <- mean_plot + labs(x = "x", y = "")

  acov_plot <- ggplot(data = data.frame(x = acov_lim), aes(x = x))
  acov_plot <- acov_plot + stat_function(fun = kdest, args = list(X = acov_est, h = acov_bw))
  acov_plot <- acov_plot + labs(x = "x", y = "")

  acor_plot <- ggplot(data = data.frame(x = acor_lim), aes(x = x))
  acor_plot <- acor_plot + stat_function(fun = kdest, args = list(X = acor_est, h = acor_bw))
  acor_plot <- acor_plot + labs(x = "x", y = "")

  # functions
  mean_func <- function(x) {
    kdest(x = x, X = mean_est, h = mean_bw)
  }

  acov_func <- function(x) {
    kdest(x = x, X = acov_est, h = acov_bw)
  }

  acor_func <- function(x) {
    kdest(x = x, X = acor_est, h = acor_bw)
  }

  # results
  bandwidth <- c(mean_bw, acov_bw, acor_bw)
  names(bandwidth) <- c("mean", "autocovariance", "autocorrelation")
  quantity <- cbind(mean_est, acov_est, acor_est)
  colnames(quantity) <- c("mean", "autocovariance", "autocorrelation")
  result <- list(mean = mean_plot, acov = acov_plot, acor = acor_plot,
                 mean_func = mean_func, acov_func = acov_func, acor_func = acor_func,
                 bandwidth = bandwidth, quantity = quantity,
                 acov_order = acov_order, acor_order = acor_order, N = N, S = S)

  return(result)

}


#' computing kernel density estimate
#'
#' @param x point at which the density is estimated
#' @param X vector of cross-sectional data
#' @param h bandwidth
#'
kdest <- Vectorize(FUN = function(x, X, h) {

  # sample size
  N <- length(X)

  # kernel density estimate
  est <- sum(dnorm( (x - X) / h)) /  (N * h)

  return(est)

}, vectorize.args = "x")