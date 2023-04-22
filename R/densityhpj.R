#' The HPJ bias-corrected kernel density estimation
#'
#' The `hpjkd()` function enables to implement the HPJ bias-corrected kernel
#' density estimation for the heterogeneous mean, the autocovariance,
#' and the autocorrelation.
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
#' panelhetero::hpjkd(data = data)
#'
#' @references Okui, R. and Yanagi, T., 2020.
#' Kernel estimation for panel data with heterogeneous dynamics.
#' The Econometrics Journal, 23(1), pp.156-175.
#'
#' @export
#'
hpjkd <- function(data,
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

  # Variable definitions -------------------------------------------------------

  # Initialization
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

  # HPJ bias-correction
  if (S %% 2 == 0) {

    # Half panel data for even T
    data1 <- data[, 1:(S / 2)]
    data2 <- data[, (S / 2 + 1):S]

    # Estimated quantities for half panel data
    mean_est1 <- rowMeans(data1)
    mean_est2 <- rowMeans(data2)

    acov_est1 <- apply(data1, MARGIN = 1, acov, acov_order = acov_order)
    acov_est2 <- apply(data2, MARGIN = 1, acov, acov_order = acov_order)

    acor_est1 <- apply(data1, MARGIN = 1, acor, acor_order = acor_order)
    acor_est2 <- apply(data2, MARGIN = 1, acor, acor_order = acor_order)

    # Make figures using ggplot2
    mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = hpjkdest1,
                             args = list(X = mean_est,
                                         X1 = mean_est1,
                                         X2 = mean_est2,
                                         h = mean_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous mean") +
      ggplot2::theme_bw()

    acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = hpjkdest1,
                             args = list(X = acov_est,
                                         X1 = acov_est1,
                                         X2 = acov_est2,
                                         h = acov_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocovariance") +
      ggplot2::theme_bw()

    acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = hpjkdest1,
                             args = list(X = acor_est,
                                         X1 = acor_est1,
                                         X2 = acor_est2,
                                         h = acor_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocorrelation") +
      ggplot2::theme_bw()

    # Functions
    mean_func <- function(x) {
      hpjkdest1(x = x,
                X = mean_est,
                X1 = mean_est1,
                X2 = mean_est2,
                h = mean_bw)
    }

    acov_func <- function(x) {
      hpjkdest1(x = x,
                X = acov_est,
                X1 = acov_est1,
                X2 = acov_est2,
                h = acov_bw)
    }

    acor_func <- function(x) {
      hpjkdest1(x = x,
                X = acor_est,
                X1 = acor_est1,
                X2 = acor_est2,
                h = acor_bw)
    }

  } else {

    # Half-panel data for odd T
    data1 <- data[, 1:floor(S / 2)]
    data2 <- data[, (floor(S / 2) + 1):S]
    data3 <- data[, 1:ceiling(S / 2)]
    data4 <- data[, (ceiling(S / 2) + 1):S]

    # Estimated quantities for half panel data
    mean_est1 <- rowMeans(data1)
    mean_est2 <- rowMeans(data2)
    mean_est3 <- rowMeans(data3)
    mean_est4 <- rowMeans(data4)

    acov_est1 <- apply(data1, MARGIN = 1, acov, acov_order = acov_order)
    acov_est2 <- apply(data2, MARGIN = 1, acov, acov_order = acov_order)
    acov_est3 <- apply(data3, MARGIN = 1, acov, acov_order = acov_order)
    acov_est4 <- apply(data4, MARGIN = 1, acov, acov_order = acov_order)

    acor_est1 <- apply(data1, MARGIN = 1, acor, acor_order = acor_order)
    acor_est2 <- apply(data2, MARGIN = 1, acor, acor_order = acor_order)
    acor_est3 <- apply(data3, MARGIN = 1, acor, acor_order = acor_order)
    acor_est4 <- apply(data4, MARGIN = 1, acor, acor_order = acor_order)

    # Make figures using ggplot2
    mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = hpjkdest2,
                             args = list(X = mean_est,
                                         X1 = mean_est1,
                                         X2 = mean_est2,
                                         X3 = mean_est3,
                                         X4 = mean_est4,
                                         h = mean_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous mean") +
      ggplot2::theme_bw()

    acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = hpjkdest2,
                             args = list(X = acov_est,
                                         X1 = acov_est1,
                                         X2 = acov_est2,
                                         X3 = acov_est3,
                                         X4 = acov_est4,
                                         h = acov_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocovariance") +
      ggplot2::theme_bw()

    acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = hpjkdest2,
                             args = list(X = acor_est,
                                         X1 = acor_est1,
                                         X2 = acor_est2,
                                         X3 = acor_est3,
                                         X4 = acor_est4,
                                         h = acor_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocorrelation") +
      ggplot2::theme_bw()

    # Functions
    mean_func <- function(x) {
      hpjkdest2(x = x,
                X = mean_est,
                X1 = mean_est1,
                X2 = mean_est2,
                X3 = mean_est3,
                X4 = mean_est4,
                h = mean_bw)
    }

    acov_func <- function(x) {
      hpjkdest2(x = x,
                X = acov_est,
                X1 = acov_est1,
                X2 = acov_est2,
                X3 = acov_est3,
                X4 = acov_est4,
                h = acov_bw)
    }

    acor_func <- function(x) {
      hpjkdest2(x = x,
                X = acor_est,
                X1 = acor_est1,
                X2 = acor_est2,
                X3 = acor_est3,
                X4 = acor_est4,
                h = acor_bw)
    }

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

#' Compute HPJ kernel density estimates for even T
#'
#' @param x An evaluation point
#' @param X A vector of cross-sectional data
#' @param X1 A vector of half-panel cross-sectional data 1
#' @param X2 A vector of half-panel cross-sectional data 2
#' @param h A scalar of bandwidth
#'
#' @returns A vector of kenrel density estimates
#'
#' @noRd
#'
hpjkdest1 <- Vectorize(FUN = function(x, X, X1, X2, h) {

  # Sample size
  N <- length(X)

  # Estimates
  est  <- sum(stats::dnorm((x - X) / h))  /  (N * h)
  est1 <- sum(stats::dnorm((x - X1) / h)) /  (N * h)
  est2 <- sum(stats::dnorm((x - X2) / h)) /  (N * h)

  # HPJ estimate
  hpjest <- 2 * est - (est1 + est2) / 2

  # Ensure non-negative estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)

  return(hpjest)

}, vectorize.args = "x")

#' Compute HPJ kernel density estimates for odd T
#'
#' @param x An evaluation point
#' @param X A vector of cross-sectional data
#' @param X1 A vector of half-panel cross-sectional data 1
#' @param X2 A vector of half-panel cross-sectional data 2
#' @param X3 A vector of half-panel cross-sectional data 3
#' @param X4 A vector of half-panel cross-sectional data 4
#' @param h A scalar of bandwidth
#'
#' @returns A vector of kernel density estimates
#'
#' @noRd
#'
hpjkdest2 <- Vectorize(FUN = function(x, X, X1, X2, X3, X4, h) {

  # Sample size
  N <- length(X)

  # Estimates
  est  <- sum(stats::dnorm((x - X)  / h)) /  (N * h)
  est1 <- sum(stats::dnorm((x - X1) / h)) /  (N * h)
  est2 <- sum(stats::dnorm((x - X2) / h)) /  (N * h)
  est3 <- sum(stats::dnorm((x - X3) / h)) /  (N * h)
  est4 <- sum(stats::dnorm((x - X4) / h)) /  (N * h)

  # HPJ estimate
  hpjest <- 2 * est - (est1 + est2 + est3 + est4) / 4

  # Ensure non-negative estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)

  return(hpjest)

}, vectorize.args = "x")
