#' The naive kernel density estimation
#'
#' The `nekd()` function enables to implement the naive kernel density
#' estimation without bias correction for the heterogeneous mean,
#' the autocovariance, and the autocorrelation.
#' The method is developed by Okui and Yanagi (2020).
#' For more details, see the package vignette with `vignette("panelhetero")`.
#'
#' @param data A matrix of panel data.
#' Each row corresponds to individual time series.
#' @param acov_order A non-negative integer of the order of autocovariance.
#' Default is 0.
#' @param acor_order A positive integer of the order of autocorrelation.
#' Default is 1.
#' @param mean_bw A scalar of bandwidth used for the estimation of
#' the denisty of mean.
#' Default is NULL, and the plug-in bandwidth is used.
#' @param acov_bw A scalar of bandwidth used for the estimation of
#' the denisty of autocovariance.
#' Default is NULL, and the plug-in bandwidth is used.
#' @param acor_bw A scalar of bandwidth used for the estimation of
#' the denisty of autocorrelation.
#' Default is NULL, and the plug-in bandwidth is used.
#'
#' @returns A list that contains the following elements:
#' \item{mean}{A plot of the corresponding density}
#' \item{acov}{A plot of the corresponding density}
#' \item{acor}{A plot of the corresponding density}
#' \item{mean_func}{A function that returns the corresponding density}
#' \item{acov_func}{A function that returns the corresponding density}
#' \item{acor_func}{A function that returns the corresponding density}
#' \item{bandwidth}{A Vector of the bandwidths}
#' \item{quantity}{A matrix of the estimated heterogeneous quantities}
#' \item{acov_order}{The order of autocovariance}
#' \item{acor_order}{The order of autocorrelation}
#' \item{N}{The number of cross-sectional units}
#' \item{S}{The length of time series}
#'
#' @examples
#' data <- panelhetero::simulation(N = 300, S = 50)
#' panelhetero::nekd(data = data)
#'
#' @references Okui, R. and Yanagi, T., 2020.
#' Kernel estimation for panel data with heterogeneous dynamics.
#' The Econometrics Journal, 23(1), pp.156-175.
#'
#' @export
#'
nekd <- function(data,
                 acov_order = 0,
                 acor_order = 1,
                 mean_bw = NULL,
                 acov_bw = NULL,
                 acor_bw = NULL) {

  # Error handling -------------------------------------------------------------

  error3(data = data,
         acov_order = acov_order,
         acor_order = acor_order,
         mean_bw = mean_bw,
         acov_bw = acov_bw,
         acor_bw = acor_bw)

  # variable definitions -------------------------------------------------------

  x <- NULL

  # Omit NA
  data <- stats::na.omit(data)

  # Sample size
  N <- nrow(data)
  S <- ncol(data)

  # Estimated means, autocovariances, autocorrelations
  mean_est <- rowMeans(data)
  acov_est <- apply(data, MARGIN = 1, acov, acov_order = acov_order)
  acor_est <- apply(data, MARGIN = 1, acor, acor_order = acor_order)

  # Plug-in bandwidth
  if (is.null(mean_bw)) {
    mean_bw <- KernSmooth::dpik(x = mean_est,
                                scalest = "minim",
                                kernel = "normal")
  }

  if (is.null(acov_bw)) {
    acov_bw <- KernSmooth::dpik(x = acov_est,
                                scalest = "minim",
                                kernel = "normal")
  }

  if (is.null(acor_bw)) {
    acor_bw <- KernSmooth::dpik(x = acor_est,
                                scalest = "minim",
                                kernel = "normal")
  }

  # Limits used for ggplot2
  mean_lim <- c(min(mean_est),
                max(mean_est))

  acov_lim <- c(min(acov_est),
                max(acov_est))

  acor_lim <- c(min(acor_est),
                max(acor_est))

  # Make figures using ggplot2
  mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_lim),
                               ggplot2::aes(x = x)) +
    ggplot2::stat_function(fun = kdest,
                           args = list(X = mean_est, h = mean_bw)) +
    ggplot2::labs(x = "x", y = "") +
    ggplot2::ggtitle("The heterogeneous mean") +
    ggplot2::theme_bw()

  acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_lim),
                               ggplot2::aes(x = x)) +
    ggplot2::stat_function(fun = kdest,
                           args = list(X = acov_est, h = acov_bw)) +
    ggplot2::labs(x = "x", y = "") +
    ggplot2::ggtitle("The heterogeneous autocovariance") +
    ggplot2::theme_bw()

  acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_lim),
                               ggplot2::aes(x = x)) +
    ggplot2::stat_function(fun = kdest,
                           args = list(X = acor_est, h = acor_bw)) +
    ggplot2::labs(x = "x", y = "") +
    ggplot2::ggtitle("The heterogeneous autocorrelation") +
    ggplot2::theme_bw()

  # Functions
  mean_func <- function(x) {
    kdest(x = x, X = mean_est, h = mean_bw)
  }

  acov_func <- function(x) {
    kdest(x = x, X = acov_est, h = acov_bw)
  }

  acor_func <- function(x) {
    kdest(x = x, X = acor_est, h = acor_bw)
  }

  # Results
  bandwidth <- c(mean_bw,
                 acov_bw,
                 acor_bw)

  quantity <- cbind(mean_est,
                    acov_est,
                    acor_est)

  names(bandwidth) <- colnames(quantity) <-
    c("mean", "autocovariance", "autocorrelation")

  return(list(mean = mean_plot,
              acov = acov_plot,
              acor = acor_plot,
              mean_func = mean_func,
              acov_func = acov_func,
              acor_func = acor_func,
              bandwidth = bandwidth,
              quantity = quantity,
              acov_order = acov_order,
              acor_order = acor_order,
              N = N,
              S = S)
  )

}

#' Compute kernel density estimates
#'
#' @param x An evaluation point
#' @param X A vector of cross-sectional data
#' @param h A scalar of bandwidth
#'
#' @returns A vector of kernel density estimates
#'
#' @noRd
#'
kdest <- Vectorize(FUN = function(x, X, h) {

  N <- length(X)
  est <- sum(stats::dnorm((x - X) / h)) / (N * h)
  return(est)

}, vectorize.args = "x")
