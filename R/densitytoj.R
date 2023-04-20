#' The TOJ bias-corrected kernel density estimation
#'
#' The `tojkd()` function enables to implement the TOJ bias-corrected kernel
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
#' panelhetero::tojkd(data = data)
#'
#' @references Okui, R. and Yanagi, T., 2020.
#' Kernel estimation for panel data with heterogeneous dynamics.
#' The Econometrics Journal, 23(1), pp.156-175.
#'
#' @export
#'
tojkd <- function(data,
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

  acov_est <- apply(X = data,
                    MARGIN = 1,
                    FUN = acov,
                    acov_order = acov_order)

  acor_est <- apply(X = data,
                    MARGIN = 1,
                    FUN = acor,
                    acor_order = acor_order)

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

  # TOJ bias-correction
  if (S %% 6 == 0) {

    # Split  panel data for T equivalent to 0 modulo 6
    data21 <- data[, 1:(S / 2)]
    data22 <- data[, (S / 2 + 1):S]
    data31 <- data[, 1:(S / 3)]
    data32 <- data[, (S / 3 + 1):(2*S / 3)]
    data33 <- data[, (2 * S / 3 + 1):S]

    # Estimated quantities for split panel data
    mean_est21 <- rowMeans(data21)
    mean_est22 <- rowMeans(data22)
    mean_est31 <- rowMeans(data31)
    mean_est32 <- rowMeans(data32)
    mean_est33 <- rowMeans(data33)

    acov_est21 <- apply(data21, MARGIN = 1, acov, acov_order = acov_order)
    acov_est22 <- apply(data22, MARGIN = 1, acov, acov_order = acov_order)
    acov_est31 <- apply(data31, MARGIN = 1, acov, acov_order = acov_order)
    acov_est32 <- apply(data32, MARGIN = 1, acov, acov_order = acov_order)
    acov_est33 <- apply(data33, MARGIN = 1, acov, acov_order = acov_order)

    acor_est21 <- apply(data21, MARGIN = 1, acor, acor_order = acor_order)
    acor_est22 <- apply(data22, MARGIN = 1, acor, acor_order = acor_order)
    acor_est31 <- apply(data31, MARGIN = 1, acor, acor_order = acor_order)
    acor_est32 <- apply(data32, MARGIN = 1, acor, acor_order = acor_order)
    acor_est33 <- apply(data33, MARGIN = 1, acor, acor_order = acor_order)

    # Make figures using ggplot2
    mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest0,
                             args = list(X = mean_est,
                                         X21 = mean_est21,
                                         X22 = mean_est22,
                                         X31 = mean_est31,
                                         X32 = mean_est32,
                                         X33 = mean_est33,
                                         h = mean_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous mean") +
      ggplot2::theme_bw()

    acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest0,
                             args = list(X = acov_est,
                                         X21 = acov_est21,
                                         X22 = acov_est22,
                                         X31 = acov_est31,
                                         X32 = acov_est32,
                                         X33 = acov_est33,
                                         h = acov_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocovariance") +
      ggplot2::theme_bw()

    acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest0,
                             args = list(X = acor_est,
                                         X21 = acor_est21,
                                         X22 = acor_est22,
                                         X31 = acor_est31,
                                         X32 = acor_est32,
                                         X33 = acor_est33,
                                         h = acor_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocorrelation") +
      ggplot2::theme_bw()

    # Functions
    mean_func <- function(x) {
      tojkdest0(x   = x,
                X   = mean_est,
                X21 = mean_est21,
                X22 = mean_est22,
                X31 = mean_est31,
                X32 = mean_est32,
                X33 = mean_est33,
                h   = mean_bw)
    }

    acov_func <- function(x) {
      tojkdest0(x   = x,
                X   = acov_est,
                X21 = acov_est21,
                X22 = acov_est22,
                X31 = acov_est31,
                X32 = acov_est32,
                X33 = acov_est33,
                h   = acov_bw)
    }

    acor_func <- function(x) {
      tojkdest0(x   = x,
                X   = acor_est,
                X21 = acor_est21,
                X22 = acor_est22,
                X31 = acor_est31,
                X32 = acor_est32,
                X33 = acor_est33,
                h   = acor_bw)
    }

  } else if (S %% 6 == 1) {

    # Split panel data for T equivalent to 1 modulo 6
    data21 <- data[, 1:floor(S / 2)]
    data22 <- data[, (floor(S / 2) + 1):S]
    data23 <- data[, 1:ceiling(S / 2)]
    data24 <- data[, (ceiling(S / 2) + 1):S]
    data31 <- data[, 1:floor(S / 3)]
    data32 <- data[, (floor(S / 3) + 1):(2 * floor(S / 3))]
    data33 <- data[, (2 * floor(S / 3) + 1):S]
    data34 <- data[, 1:floor(S / 3)]
    data35 <- data[, (floor(S / 3) + 1):(2 * floor(S / 3) + 1)]
    data36 <- data[, (2 * floor(S / 3) + 2):S]
    data37 <- data[, 1:ceiling(S / 3)]
    data38 <- data[, (ceiling(S / 3) + 1):(2 * floor(S / 3) + 1)]
    data39 <- data[, (2 * floor(S / 3) + 2):S]

    # Estimated quantities for split panel data
    mean_est21 <- rowMeans(data21)
    mean_est22 <- rowMeans(data22)
    mean_est23 <- rowMeans(data23)
    mean_est24 <- rowMeans(data24)
    mean_est31 <- rowMeans(data31)
    mean_est32 <- rowMeans(data32)
    mean_est33 <- rowMeans(data33)
    mean_est34 <- rowMeans(data34)
    mean_est35 <- rowMeans(data35)
    mean_est36 <- rowMeans(data36)
    mean_est37 <- rowMeans(data37)
    mean_est38 <- rowMeans(data38)
    mean_est39 <- rowMeans(data39)

    acov_est21 <- apply(data21, MARGIN = 1, acov, acov_order = acov_order)
    acov_est22 <- apply(data22, MARGIN = 1, acov, acov_order = acov_order)
    acov_est23 <- apply(data23, MARGIN = 1, acov, acov_order = acov_order)
    acov_est24 <- apply(data24, MARGIN = 1, acov, acov_order = acov_order)
    acov_est31 <- apply(data31, MARGIN = 1, acov, acov_order = acov_order)
    acov_est32 <- apply(data32, MARGIN = 1, acov, acov_order = acov_order)
    acov_est33 <- apply(data33, MARGIN = 1, acov, acov_order = acov_order)
    acov_est34 <- apply(data34, MARGIN = 1, acov, acov_order = acov_order)
    acov_est35 <- apply(data35, MARGIN = 1, acov, acov_order = acov_order)
    acov_est36 <- apply(data36, MARGIN = 1, acov, acov_order = acov_order)
    acov_est37 <- apply(data37, MARGIN = 1, acov, acov_order = acov_order)
    acov_est38 <- apply(data38, MARGIN = 1, acov, acov_order = acov_order)
    acov_est39 <- apply(data39, MARGIN = 1, acov, acov_order = acov_order)

    acor_est21 <- apply(data21, MARGIN = 1, acor, acor_order = acor_order)
    acor_est22 <- apply(data22, MARGIN = 1, acor, acor_order = acor_order)
    acor_est23 <- apply(data23, MARGIN = 1, acor, acor_order = acor_order)
    acor_est24 <- apply(data24, MARGIN = 1, acor, acor_order = acor_order)
    acor_est31 <- apply(data31, MARGIN = 1, acor, acor_order = acor_order)
    acor_est32 <- apply(data32, MARGIN = 1, acor, acor_order = acor_order)
    acor_est33 <- apply(data33, MARGIN = 1, acor, acor_order = acor_order)
    acor_est34 <- apply(data34, MARGIN = 1, acor, acor_order = acor_order)
    acor_est35 <- apply(data35, MARGIN = 1, acor, acor_order = acor_order)
    acor_est36 <- apply(data36, MARGIN = 1, acor, acor_order = acor_order)
    acor_est37 <- apply(data37, MARGIN = 1, acor, acor_order = acor_order)
    acor_est38 <- apply(data38, MARGIN = 1, acor, acor_order = acor_order)
    acor_est39 <- apply(data39, MARGIN = 1, acor, acor_order = acor_order)

    # Make figures using ggplot2
    mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest1,
                             args = list(X   = mean_est,
                                         X21 = mean_est21,
                                         X22 = mean_est22,
                                         X23 = mean_est23,
                                         X24 = mean_est24,
                                         X31 = mean_est31,
                                         X32 = mean_est32,
                                         X33 = mean_est33,
                                         X34 = mean_est34,
                                         X35 = mean_est35,
                                         X36 = mean_est36,
                                         X37 = mean_est37,
                                         X38 = mean_est38,
                                         X39 = mean_est39,
                                         h   = mean_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous mean") +
      ggplot2::theme_bw()

    acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest1,
                             args = list(X   = acov_est,
                                         X21 = acov_est21,
                                         X22 = acov_est22,
                                         X23 = acov_est23,
                                         X24 = acov_est24,
                                         X31 = acov_est31,
                                         X32 = acov_est32,
                                         X33 = acov_est33,
                                         X34 = acov_est34,
                                         X35 = acov_est35,
                                         X36 = acov_est36,
                                         X37 = acov_est37,
                                         X38 = acov_est38,
                                         X39 = acov_est39,
                                         h   = acov_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocovariance") +
      ggplot2::theme_bw()

    acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest1,
                             args = list(X   = acor_est,
                                         X21 = acor_est21,
                                         X22 = acor_est22,
                                         X23 = acor_est23,
                                         X24 = acor_est24,
                                         X31 = acor_est31,
                                         X32 = acor_est32,
                                         X33 = acor_est33,
                                         X34 = acor_est34,
                                         X35 = acor_est35,
                                         X36 = acor_est36,
                                         X37 = acor_est37,
                                         X38 = acor_est38,
                                         X39 = acor_est39,
                                         h   = acor_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocorrelation") +
      ggplot2::theme_bw()

    # Functions
    mean_func <- function(x) {
      tojkdest1(x   = x,
                X   = mean_est,
                X21 = mean_est21,
                X22 = mean_est22,
                X23 = mean_est23,
                X24 = mean_est24,
                X31 = mean_est31,
                X32 = mean_est32,
                X33 = mean_est33,
                X34 = mean_est34,
                X35 = mean_est35,
                X36 = mean_est36,
                X37 = mean_est37,
                X38 = mean_est38,
                X39 = mean_est39,
                h = mean_bw)
    }

    acov_func <- function(x) {
      tojkdest1(x   = x,
                X   = acov_est,
                X21 = acov_est21,
                X22 = acov_est22,
                X23 = acov_est23,
                X24 = acov_est24,
                X31 = acov_est31,
                X32 = acov_est32,
                X33 = acov_est33,
                X34 = acov_est34,
                X35 = acov_est35,
                X36 = acov_est36,
                X37 = acov_est37,
                X38 = acov_est38,
                X39 = acov_est39,
                h   = acov_bw)
    }

    acor_func <- function(x) {
      tojkdest1(x   = x,
                X   = acor_est,
                X21 = acor_est21,
                X22 = acor_est22,
                X23 = acor_est23,
                X24 = acor_est24,
                X31 = acor_est31,
                X32 = acor_est32,
                X33 = acor_est33,
                X34 = acor_est34,
                X35 = acor_est35,
                X36 = acor_est36,
                X37 = acor_est37,
                X38 = acor_est38,
                X39 = acor_est39,
                h   = acor_bw)
    }

  } else if (S %% 6 == 2) {

    # Split  panel data for T equivalent to 2 modulo 6
    data21 <- data[, 1:(S / 2)]
    data22 <- data[, (S / 2 + 1):S]
    data31 <- data[, 1:floor(S / 3)]
    data32 <- data[, (floor(S / 3) + 1):(2 * floor(S / 3) + 1) ]
    data33 <- data[, (2 * ceiling(S / 3)):S]
    data34 <- data[, 1:ceiling(S / 3)]
    data35 <- data[, (ceiling(S / 3) + 1):(2 * floor(S / 3) + 1)]
    data36 <- data[, (2 * ceiling(S / 3)):S]
    data37 <- data[, 1:ceiling(S / 3)]
    data38 <- data[, (ceiling(S / 3) + 1):(2 * ceiling(S / 3))]
    data39 <- data[, (2 * ceiling(S / 3) + 1):S]

    # Estimated quantities for split panel data
    mean_est21 <- rowMeans(data21)
    mean_est22 <- rowMeans(data22)
    mean_est31 <- rowMeans(data31)
    mean_est32 <- rowMeans(data32)
    mean_est33 <- rowMeans(data33)
    mean_est34 <- rowMeans(data34)
    mean_est35 <- rowMeans(data35)
    mean_est36 <- rowMeans(data36)
    mean_est37 <- rowMeans(data37)
    mean_est38 <- rowMeans(data38)
    mean_est39 <- rowMeans(data39)

    acov_est21 <- apply(data21, MARGIN = 1, acov, acov_order = acov_order)
    acov_est22 <- apply(data22, MARGIN = 1, acov, acov_order = acov_order)
    acov_est31 <- apply(data31, MARGIN = 1, acov, acov_order = acov_order)
    acov_est32 <- apply(data32, MARGIN = 1, acov, acov_order = acov_order)
    acov_est33 <- apply(data33, MARGIN = 1, acov, acov_order = acov_order)
    acov_est34 <- apply(data34, MARGIN = 1, acov, acov_order = acov_order)
    acov_est35 <- apply(data35, MARGIN = 1, acov, acov_order = acov_order)
    acov_est36 <- apply(data36, MARGIN = 1, acov, acov_order = acov_order)
    acov_est37 <- apply(data37, MARGIN = 1, acov, acov_order = acov_order)
    acov_est38 <- apply(data38, MARGIN = 1, acov, acov_order = acov_order)
    acov_est39 <- apply(data39, MARGIN = 1, acov, acov_order = acov_order)

    acor_est21 <- apply(data21, MARGIN = 1, acor, acor_order = acor_order)
    acor_est22 <- apply(data22, MARGIN = 1, acor, acor_order = acor_order)
    acor_est31 <- apply(data31, MARGIN = 1, acor, acor_order = acor_order)
    acor_est32 <- apply(data32, MARGIN = 1, acor, acor_order = acor_order)
    acor_est33 <- apply(data33, MARGIN = 1, acor, acor_order = acor_order)
    acor_est34 <- apply(data34, MARGIN = 1, acor, acor_order = acor_order)
    acor_est35 <- apply(data35, MARGIN = 1, acor, acor_order = acor_order)
    acor_est36 <- apply(data36, MARGIN = 1, acor, acor_order = acor_order)
    acor_est37 <- apply(data37, MARGIN = 1, acor, acor_order = acor_order)
    acor_est38 <- apply(data38, MARGIN = 1, acor, acor_order = acor_order)
    acor_est39 <- apply(data39, MARGIN = 1, acor, acor_order = acor_order)

    # Make figures using ggplot2
    mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest2,
                             args = list(X   = mean_est,
                                         X21 = mean_est21,
                                         X22 = mean_est22,
                                         X31 = mean_est31,
                                         X32 = mean_est32,
                                         X33 = mean_est33,
                                         X34 = mean_est34,
                                         X35 = mean_est35,
                                         X36 = mean_est36,
                                         X37 = mean_est37,
                                         X38 = mean_est38,
                                         X39 = mean_est39,
                                         h   = mean_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous mean") +
      ggplot2::theme_bw()

    acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest2,
                             args = list(X   = acov_est,
                                         X21 = acov_est21,
                                         X22 = acov_est22,
                                         X31 = acov_est31,
                                         X32 = acov_est32,
                                         X33 = acov_est33,
                                         X34 = acov_est34,
                                         X35 = acov_est35,
                                         X36 = acov_est36,
                                         X37 = acov_est37,
                                         X38 = acov_est38,
                                         X39 = acov_est39,
                                         h   = acov_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocovariance") +
      ggplot2::theme_bw()

    acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest2,
                             args = list(X   = acor_est,
                                         X21 = acor_est21,
                                         X22 = acor_est22,
                                         X31 = acor_est31,
                                         X32 = acor_est32,
                                         X33 = acor_est33,
                                         X34 = acor_est34,
                                         X35 = acor_est35,
                                         X36 = acor_est36,
                                         X37 = acor_est37,
                                         X38 = acor_est38,
                                         X39 = acor_est39,
                                         h   = acor_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocorrelation") +
      ggplot2::theme_bw()

    # Functions
    mean_func <- function(x) {
      tojkdest2(x   = x,
                X   = mean_est,
                X21 = mean_est21,
                X22 = mean_est22,
                X31 = mean_est31,
                X32 = mean_est32,
                X33 = mean_est33,
                X34 = mean_est34,
                X35 = mean_est35,
                X36 = mean_est36,
                X37 = mean_est37,
                X38 = mean_est38,
                X39 = mean_est39,
                h   = mean_bw)
    }

    acov_func <- function(x) {
      tojkdest2(x   = x,
                X   = acov_est,
                X21 = acov_est21,
                X22 = acov_est22,
                X31 = acov_est31,
                X32 = acov_est32,
                X33 = acov_est33,
                X34 = acov_est34,
                X35 = acov_est35,
                X36 = acov_est36,
                X37 = acov_est37,
                X38 = acov_est38,
                X39 = acov_est39,
                h   = acov_bw)
    }

    acor_func <- function(x) {
      tojkdest2(x   = x,
                X   = acor_est,
                X21 = acor_est21,
                X22 = acor_est22,
                X31 = acor_est31,
                X32 = acor_est32,
                X33 = acor_est33,
                X34 = acor_est34,
                X35 = acor_est35,
                X36 = acor_est36,
                X37 = acor_est37,
                X38 = acor_est38,
                X39 = acor_est39,
                h = acor_bw)
    }

  } else if (S %% 6 == 3) {

    # Split  panel data for T equivalent to 3 modulo 6
    data21 <- data[, 1:floor(S / 2)]
    data22 <- data[, (floor(S / 2) + 1):S]
    data23 <- data[, 1:ceiling(S / 2)]
    data24 <- data[, (ceiling(S / 2) + 1):S]
    data31 <- data[, 1:(S / 3)]
    data32 <- data[, (S / 3 + 1):(2*S / 3)]
    data33 <- data[, (2 * S / 3 + 1):S]

    # Estimated quantities for split panel data
    mean_est21 <- rowMeans(data21)
    mean_est22 <- rowMeans(data22)
    mean_est23 <- rowMeans(data23)
    mean_est24 <- rowMeans(data24)
    mean_est31 <- rowMeans(data31)
    mean_est32 <- rowMeans(data32)
    mean_est33 <- rowMeans(data33)

    acov_est21 <- apply(data21, MARGIN = 1, acov, acov_order = acov_order)
    acov_est22 <- apply(data22, MARGIN = 1, acov, acov_order = acov_order)
    acov_est23 <- apply(data23, MARGIN = 1, acov, acov_order = acov_order)
    acov_est24 <- apply(data24, MARGIN = 1, acov, acov_order = acov_order)
    acov_est31 <- apply(data31, MARGIN = 1, acov, acov_order = acov_order)
    acov_est32 <- apply(data32, MARGIN = 1, acov, acov_order = acov_order)
    acov_est33 <- apply(data33, MARGIN = 1, acov, acov_order = acov_order)

    acor_est21 <- apply(data21, MARGIN = 1, acor, acor_order = acor_order)
    acor_est22 <- apply(data22, MARGIN = 1, acor, acor_order = acor_order)
    acor_est23 <- apply(data23, MARGIN = 1, acor, acor_order = acor_order)
    acor_est24 <- apply(data24, MARGIN = 1, acor, acor_order = acor_order)
    acor_est31 <- apply(data31, MARGIN = 1, acor, acor_order = acor_order)
    acor_est32 <- apply(data32, MARGIN = 1, acor, acor_order = acor_order)
    acor_est33 <- apply(data33, MARGIN = 1, acor, acor_order = acor_order)

    # Make figures using ggplot2
    mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest3,
                             args = list(X   = mean_est,
                                         X21 = mean_est21,
                                         X22 = mean_est22,
                                         X23 = mean_est23,
                                         X24 = mean_est24,
                                         X31 = mean_est31,
                                         X32 = mean_est32,
                                         X33 = mean_est33,
                                         h   = mean_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous mean") +
      ggplot2::theme_bw()

    acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest3,
                             args = list(X   = acov_est,
                                         X21 = acov_est21,
                                         X22 = acov_est22,
                                         X23 = acov_est23,
                                         X24 = acov_est24,
                                         X31 = acov_est31,
                                         X32 = acov_est32,
                                         X33 = acov_est33,
                                         h   = acov_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocovariance") +
      ggplot2::theme_bw()

    acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest3,
                             args = list(X   = acor_est,
                                         X21 = acor_est21,
                                         X22 = acor_est22,
                                         X23 = acor_est23,
                                         X24 = acor_est24,
                                         X31 = acor_est31,
                                         X32 = acor_est32,
                                         X33 = acor_est33,
                                         h   = acor_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocorrelation") +
      ggplot2::theme_bw()

    # Functions
    mean_func <- function(x) {
      tojkdest3(x   = x,
                X   = mean_est,
                X21 = mean_est21,
                X22 = mean_est22,
                X23 = mean_est23,
                X24 = mean_est24,
                X31 = mean_est31,
                X32 = mean_est32,
                X33 = mean_est33,
                h   = mean_bw)
    }

    acov_func <- function(x) {
      tojkdest3(x   = x,
                X   = acov_est,
                X21 = acov_est21,
                X22 = acov_est22,
                X23 = acov_est23,
                X24 = acov_est24,
                X31 = acov_est31,
                X32 = acov_est32,
                X33 = acov_est33,
                h   = acov_bw)
    }

    acor_func <- function(x) {
      tojkdest3(x   = x,
                X   = acor_est,
                X21 = acor_est21,
                X22 = acor_est22,
                X23 = acor_est23,
                X24 = acor_est24,
                X31 = acor_est31,
                X32 = acor_est32,
                X33 = acor_est33,
                h   = acor_bw)
    }

  } else if (S %% 6 == 4) {

    # Split  panel data for T equivalent to 4 modulo 6
    data21 <- data[, 1:(S / 2)]
    data22 <- data[, (S / 2 + 1):S]
    data31 <- data[, 1:floor(S / 3)]
    data32 <- data[, (floor(S / 3) + 1):(2 * floor(S / 3))]
    data33 <- data[, (2 * floor(S / 3) + 1):S]
    data34 <- data[, 1:floor(S / 3)]
    data35 <- data[, (floor(S / 3) + 1):(2 * floor(S / 3) + 1)]
    data36 <- data[, (2 * floor(S / 3) + 2):S]
    data37 <- data[, 1:ceiling(S / 3)]
    data38 <- data[, (ceiling(S / 3) + 1):(2 * floor(S / 3) + 1)]
    data39 <- data[, (2 * floor(S / 3) + 2):S]

    # Estimated quantities for split panel data
    mean_est21 <- rowMeans(data21)
    mean_est22 <- rowMeans(data22)
    mean_est31 <- rowMeans(data31)
    mean_est32 <- rowMeans(data32)
    mean_est33 <- rowMeans(data33)
    mean_est34 <- rowMeans(data34)
    mean_est35 <- rowMeans(data35)
    mean_est36 <- rowMeans(data36)
    mean_est37 <- rowMeans(data37)
    mean_est38 <- rowMeans(data38)
    mean_est39 <- rowMeans(data39)

    acov_est21 <- apply(data21, MARGIN = 1, acov, acov_order = acov_order)
    acov_est22 <- apply(data22, MARGIN = 1, acov, acov_order = acov_order)
    acov_est31 <- apply(data31, MARGIN = 1, acov, acov_order = acov_order)
    acov_est32 <- apply(data32, MARGIN = 1, acov, acov_order = acov_order)
    acov_est33 <- apply(data33, MARGIN = 1, acov, acov_order = acov_order)
    acov_est34 <- apply(data34, MARGIN = 1, acov, acov_order = acov_order)
    acov_est35 <- apply(data35, MARGIN = 1, acov, acov_order = acov_order)
    acov_est36 <- apply(data36, MARGIN = 1, acov, acov_order = acov_order)
    acov_est37 <- apply(data37, MARGIN = 1, acov, acov_order = acov_order)
    acov_est38 <- apply(data38, MARGIN = 1, acov, acov_order = acov_order)
    acov_est39 <- apply(data39, MARGIN = 1, acov, acov_order = acov_order)

    acor_est21 <- apply(data21, MARGIN = 1, acor, acor_order = acor_order)
    acor_est22 <- apply(data22, MARGIN = 1, acor, acor_order = acor_order)
    acor_est31 <- apply(data31, MARGIN = 1, acor, acor_order = acor_order)
    acor_est32 <- apply(data32, MARGIN = 1, acor, acor_order = acor_order)
    acor_est33 <- apply(data33, MARGIN = 1, acor, acor_order = acor_order)
    acor_est34 <- apply(data34, MARGIN = 1, acor, acor_order = acor_order)
    acor_est35 <- apply(data35, MARGIN = 1, acor, acor_order = acor_order)
    acor_est36 <- apply(data36, MARGIN = 1, acor, acor_order = acor_order)
    acor_est37 <- apply(data37, MARGIN = 1, acor, acor_order = acor_order)
    acor_est38 <- apply(data38, MARGIN = 1, acor, acor_order = acor_order)
    acor_est39 <- apply(data39, MARGIN = 1, acor, acor_order = acor_order)

    # Make figures using ggplot2
    mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest4,
                             args = list(X   = mean_est,
                                         X21 = mean_est21,
                                         X22 = mean_est22,
                                         X31 = mean_est31,
                                         X32 = mean_est32,
                                         X33 = mean_est33,
                                         X34 = mean_est34,
                                         X35 = mean_est35,
                                         X36 = mean_est36,
                                         X37 = mean_est37,
                                         X38 = mean_est38,
                                         X39 = mean_est39,
                                         h   = mean_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous mean") +
      ggplot2::theme_bw()

    acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest4,
                             args = list(X   = acov_est,
                                         X21 = acov_est21,
                                         X22 = acov_est22,
                                         X31 = acov_est31,
                                         X32 = acov_est32,
                                         X33 = acov_est33,
                                         X34 = acov_est34,
                                         X35 = acov_est35,
                                         X36 = acov_est36,
                                         X37 = acov_est37,
                                         X38 = acov_est38,
                                         X39 = acov_est39,
                                         h   = acov_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocovariance") +
      ggplot2::theme_bw()

    acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest4,
                             args = list(X   = acor_est,
                                         X21 = acor_est21,
                                         X22 = acor_est22,
                                         X31 = acor_est31,
                                         X32 = acor_est32,
                                         X33 = acor_est33,
                                         X34 = acor_est34,
                                         X35 = acor_est35,
                                         X36 = acor_est36,
                                         X37 = acor_est37,
                                         X38 = acor_est38,
                                         X39 = acor_est39,
                                         h   = acor_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocorrelation") +
      ggplot2::theme_bw()

    # Functions
    mean_func <- function(x) {
      tojkdest4(x   = x,
                X   = mean_est,
                X21 = mean_est21,
                X22 = mean_est22,
                X31 = mean_est31,
                X32 = mean_est32,
                X33 = mean_est33,
                X34 = mean_est34,
                X35 = mean_est35,
                X36 = mean_est36,
                X37 = mean_est37,
                X38 = mean_est38,
                X39 = mean_est39,
                h   = mean_bw)
    }

    acov_func <- function(x) {
      tojkdest4(x   = x,
                X   = acov_est,
                X21 = acov_est21,
                X22 = acov_est22,
                X31 = acov_est31,
                X32 = acov_est32,
                X33 = acov_est33,
                X34 = acov_est34,
                X35 = acov_est35,
                X36 = acov_est36,
                X37 = acov_est37,
                X38 = acov_est38,
                X39 = acov_est39,
                h   = acov_bw)
    }

    acor_func <- function(x) {
      tojkdest4(x   = x,
                X   = acor_est,
                X21 = acor_est21,
                X22 = acor_est22,
                X31 = acor_est31,
                X32 = acor_est32,
                X33 = acor_est33,
                X34 = acor_est34,
                X35 = acor_est35,
                X36 = acor_est36,
                X37 = acor_est37,
                X38 = acor_est38,
                X39 = acor_est39,
                h   = acor_bw)
    }

  } else {

    # Split  panel data for T equivalent to 5 modulo 6
    data21 <- data[, 1:floor(S / 2)]
    data22 <- data[, (floor(S / 2) + 1):S]
    data23 <- data[, 1:ceiling(S / 2)]
    data24 <- data[, (ceiling(S / 2) + 1):S]
    data31 <- data[, 1:floor(S / 3)]
    data32 <- data[, (floor(S / 3) + 1):(2 * floor(S / 3) + 1) ]
    data33 <- data[, (2 * ceiling(S / 3)):S]
    data34 <- data[, 1:ceiling(S / 3)]
    data35 <- data[, (ceiling(S / 3) + 1):(2 * floor(S / 3) + 1)]
    data36 <- data[, (2 * ceiling(S / 3)):S]
    data37 <- data[, 1:ceiling(S / 3)]
    data38 <- data[, (ceiling(S / 3) + 1):(2 * ceiling(S / 3))]
    data39 <- data[, (2 * ceiling(S / 3) + 1):S]

    # Estimated quantities for split panel data
    mean_est21 <- rowMeans(data21)
    mean_est22 <- rowMeans(data22)
    mean_est23 <- rowMeans(data23)
    mean_est24 <- rowMeans(data24)
    mean_est31 <- rowMeans(data31)
    mean_est32 <- rowMeans(data32)
    mean_est33 <- rowMeans(data33)
    mean_est34 <- rowMeans(data34)
    mean_est35 <- rowMeans(data35)
    mean_est36 <- rowMeans(data36)
    mean_est37 <- rowMeans(data37)
    mean_est38 <- rowMeans(data38)
    mean_est39 <- rowMeans(data39)

    acov_est21 <- apply(data21, MARGIN = 1, acov, acov_order = acov_order)
    acov_est22 <- apply(data22, MARGIN = 1, acov, acov_order = acov_order)
    acov_est23 <- apply(data23, MARGIN = 1, acov, acov_order = acov_order)
    acov_est24 <- apply(data24, MARGIN = 1, acov, acov_order = acov_order)
    acov_est31 <- apply(data31, MARGIN = 1, acov, acov_order = acov_order)
    acov_est32 <- apply(data32, MARGIN = 1, acov, acov_order = acov_order)
    acov_est33 <- apply(data33, MARGIN = 1, acov, acov_order = acov_order)
    acov_est34 <- apply(data34, MARGIN = 1, acov, acov_order = acov_order)
    acov_est35 <- apply(data35, MARGIN = 1, acov, acov_order = acov_order)
    acov_est36 <- apply(data36, MARGIN = 1, acov, acov_order = acov_order)
    acov_est37 <- apply(data37, MARGIN = 1, acov, acov_order = acov_order)
    acov_est38 <- apply(data38, MARGIN = 1, acov, acov_order = acov_order)
    acov_est39 <- apply(data39, MARGIN = 1, acov, acov_order = acov_order)

    acor_est21 <- apply(data21, MARGIN = 1, acor, acor_order = acor_order)
    acor_est22 <- apply(data22, MARGIN = 1, acor, acor_order = acor_order)
    acor_est23 <- apply(data23, MARGIN = 1, acor, acor_order = acor_order)
    acor_est24 <- apply(data24, MARGIN = 1, acor, acor_order = acor_order)
    acor_est31 <- apply(data31, MARGIN = 1, acor, acor_order = acor_order)
    acor_est32 <- apply(data32, MARGIN = 1, acor, acor_order = acor_order)
    acor_est33 <- apply(data33, MARGIN = 1, acor, acor_order = acor_order)
    acor_est34 <- apply(data34, MARGIN = 1, acor, acor_order = acor_order)
    acor_est35 <- apply(data35, MARGIN = 1, acor, acor_order = acor_order)
    acor_est36 <- apply(data36, MARGIN = 1, acor, acor_order = acor_order)
    acor_est37 <- apply(data37, MARGIN = 1, acor, acor_order = acor_order)
    acor_est38 <- apply(data38, MARGIN = 1, acor, acor_order = acor_order)
    acor_est39 <- apply(data39, MARGIN = 1, acor, acor_order = acor_order)

    # Make figures by ggplot2
    mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest5,
                             args = list(X   = mean_est,
                                         X21 = mean_est21,
                                         X22 = mean_est22,
                                         X23 = mean_est23,
                                         X24 = mean_est24,
                                         X31 = mean_est31,
                                         X32 = mean_est32,
                                         X33 = mean_est33,
                                         X34 = mean_est34,
                                         X35 = mean_est35,
                                         X36 = mean_est36,
                                         X37 = mean_est37,
                                         X38 = mean_est38,
                                         X39 = mean_est39,
                                         h   = mean_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous mean") +
      ggplot2::theme_bw()

    acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest5,
                             args = list(X   = acov_est,
                                         X21 = acov_est21,
                                         X22 = acov_est22,
                                         X23 = acov_est23,
                                         X24 = acov_est24,
                                         X31 = acov_est31,
                                         X32 = acov_est32,
                                         X33 = acov_est33,
                                         X34 = acov_est34,
                                         X35 = acov_est35,
                                         X36 = acov_est36,
                                         X37 = acov_est37,
                                         X38 = acov_est38,
                                         X39 = acov_est39,
                                         h   = acov_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocovariance") +
      ggplot2::theme_bw()

    acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = tojkdest5,
                             args = list(X   = acor_est,
                                         X21 = acor_est21,
                                         X22 = acor_est22,
                                         X23 = acor_est23,
                                         X24 = acor_est24,
                                         X31 = acor_est31,
                                         X32 = acor_est32,
                                         X33 = acor_est33,
                                         X34 = acor_est34,
                                         X35 = acor_est35,
                                         X36 = acor_est36,
                                         X37 = acor_est37,
                                         X38 = acor_est38,
                                         X39 = acor_est39,
                                         h   = acor_bw)) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocorrelation") +
      ggplot2::theme_bw()

    # Functions
    mean_func <- function(x) {
      tojkdest5(x   = x,
                X   = mean_est,
                X21 = mean_est21,
                X22 = mean_est22,
                X23 = mean_est23,
                X24 = mean_est24,
                X31 = mean_est31,
                X32 = mean_est32,
                X33 = mean_est33,
                X34 = mean_est34,
                X35 = mean_est35,
                X36 = mean_est36,
                X37 = mean_est37,
                X38 = mean_est38,
                X39 = mean_est39,
                h   = mean_bw)
    }

    acov_func <- function(x) {
      tojkdest5(x   = x,
                X   = acov_est,
                X21 = acov_est21,
                X22 = acov_est22,
                X23 = acov_est23,
                X24 = acov_est24,
                X31 = acov_est31,
                X32 = acov_est32,
                X33 = acov_est33,
                X34 = acov_est34,
                X35 = acov_est35,
                X36 = acov_est36,
                X37 = acov_est37,
                X38 = acov_est38,
                X39 = acov_est39,
                h   = acov_bw)
    }

    acor_func <- function(x) {
      tojkdest5(x   = x,
                X   = acor_est,
                X21 = acor_est21,
                X22 = acor_est22,
                X23 = acor_est23,
                X24 = acor_est24,
                X31 = acor_est31,
                X32 = acor_est32,
                X33 = acor_est33,
                X34 = acor_est34,
                X35 = acor_est35,
                X36 = acor_est36,
                X37 = acor_est37,
                X38 = acor_est38,
                X39 = acor_est39,
                h   = acor_bw)
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

#' Compute TOJ kernel density estimates for T equivalent to 0 modulo 6
#'
#' @param x An evaluation point
#' @param X A vector of cross-sectional data
#' @param X21 A vector of half-panel cross-sectional data 21
#' @param X22 A vector of half-panel cross-sectional data 22
#' @param X31 A vector of third-panel cross-sectional data 31
#' @param X32 A vector of third-panel cross-sectional data 32
#' @param X33 A vector of third-panel cross-sectional data 33
#' @param h A scalar of bandwidth
#'
#' @returns A vector of kernel density estimates
#'
tojkdest0 <- Vectorize(FUN = function(x, X, X21, X22, X31, X32, X33, h) {

  # Sample size
  N <- length(X)

  # Estimates
  est   <- sum(stats::dnorm((x - X)   / h)) / (N * h)
  est21 <- sum(stats::dnorm((x - X21) / h)) / (N * h)
  est22 <- sum(stats::dnorm((x - X22) / h)) / (N * h)
  est31 <- sum(stats::dnorm((x - X31) / h)) / (N * h)
  est32 <- sum(stats::dnorm((x - X32) / h)) / (N * h)
  est33 <- sum(stats::dnorm((x - X33) / h)) / (N * h)

  # TOJ estimates
  tojest <- 3.536 * est -
    4.072 * (est21 + est22) / 2 +
    1.536 * (est31 + est32 + est33) / 3

  # Ensure non-negative estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)

  return(tojest)

}, vectorize.args = "x")

#' Compute TOJ kernel density estimate for T equivalent to 1 modulo 6
#'
#' @param x An evaluation point
#' @param X A vector of cross-sectional data
#' @param X21 A vector of half-panel cross-sectional data 21
#' @param X22 A vector of half-panel cross-sectional data 22
#' @param X23 A vector of half-panel cross-sectional data 23
#' @param X24 A vector of half-panel cross-sectional data 24
#' @param X31 A vector of third-panel cross-sectional data 31
#' @param X32 A vector of third-panel cross-sectional data 32
#' @param X33 A vector of third-panel cross-sectional data 33
#' @param X34 A vector of third-panel cross-sectional data 34
#' @param X35 A vector of third-panel cross-sectional data 35
#' @param X36 A vector of third-panel cross-sectional data 36
#' @param X37 A vector of third-panel cross-sectional data 37
#' @param X38 A vector of third-panel cross-sectional data 38
#' @param X39 A vector of third-panel cross-sectional data 39
#' @param h A scalar of bandwidth
#'
#' @returns A vector of kernel density estimates
#'
tojkdest1 <- Vectorize(FUN = function(x,
                                      X,
                                      X21,
                                      X22,
                                      X23,
                                      X24,
                                      X31,
                                      X32,
                                      X33,
                                      X34,
                                      X35,
                                      X36,
                                      X37,
                                      X38,
                                      X39,
                                      h) {

  # Sample size
  N <- length(X)

  # Estimates
  est   <- sum(stats::dnorm((x - X)   / h)) / (N * h)
  est21 <- sum(stats::dnorm((x - X21) / h)) / (N * h)
  est22 <- sum(stats::dnorm((x - X22) / h)) / (N * h)
  est23 <- sum(stats::dnorm((x - X23) / h)) / (N * h)
  est24 <- sum(stats::dnorm((x - X24) / h)) / (N * h)
  est31 <- sum(stats::dnorm((x - X31) / h)) / (N * h)
  est32 <- sum(stats::dnorm((x - X32) / h)) / (N * h)
  est33 <- sum(stats::dnorm((x - X33) / h)) / (N * h)
  est34 <- sum(stats::dnorm((x - X34) / h)) / (N * h)
  est35 <- sum(stats::dnorm((x - X35) / h)) / (N * h)
  est36 <- sum(stats::dnorm((x - X36) / h)) / (N * h)
  est37 <- sum(stats::dnorm((x - X37) / h)) / (N * h)
  est38 <- sum(stats::dnorm((x - X38) / h)) / (N * h)
  est39 <- sum(stats::dnorm((x - X39) / h)) / (N * h)

  # TOJ estimates
  tojest <- 3.536 * est -
    4.072 * (est21 + est22 + est23 + est24) / 4 +
    1.536 * (est31 + est32 + est33 + est34 +
               est35 + est36 + est37 + est38 + est39) / 9

  # Ensure non-negative estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)

  return(tojest)

}, vectorize.args = "x")

#' Compute TOJ kernel density estimate for T equivalent to 2 modulo 6
#'
#' @param x An evaluation point
#' @param X A vector of cross-sectional data
#' @param X21 A vector of half-panel cross-sectional data 21
#' @param X22 A vector of half-panel cross-sectional data 22
#' @param X31 A vector of third-panel cross-sectional data 31
#' @param X32 A vector of third-panel cross-sectional data 32
#' @param X33 A vector of third-panel cross-sectional data 33
#' @param X34 A vector of third-panel cross-sectional data 34
#' @param X35 A vector of third-panel cross-sectional data 35
#' @param X36 A vector of third-panel cross-sectional data 36
#' @param X37 A vector of third-panel cross-sectional data 37
#' @param X38 A vector of third-panel cross-sectional data 38
#' @param X39 A vector of third-panel cross-sectional data 39
#' @param h A scalar of bandwidth
#'
#' @returns A vector of kernel density estimates
#'
tojkdest2 <- Vectorize(FUN = function(x,
                                      X,
                                      X21,
                                      X22,
                                      X31,
                                      X32,
                                      X33,
                                      X34,
                                      X35,
                                      X36,
                                      X37,
                                      X38,
                                      X39,
                                      h) {

  # Sample size
  N <- length(X)

  # Estimates
  est   <- sum(stats::dnorm((x - X)   / h)) / (N * h)
  est21 <- sum(stats::dnorm((x - X21) / h)) / (N * h)
  est22 <- sum(stats::dnorm((x - X22) / h)) / (N * h)
  est31 <- sum(stats::dnorm((x - X31) / h)) / (N * h)
  est32 <- sum(stats::dnorm((x - X32) / h)) / (N * h)
  est33 <- sum(stats::dnorm((x - X33) / h)) / (N * h)
  est34 <- sum(stats::dnorm((x - X34) / h)) / (N * h)
  est35 <- sum(stats::dnorm((x - X35) / h)) / (N * h)
  est36 <- sum(stats::dnorm((x - X36) / h)) / (N * h)
  est37 <- sum(stats::dnorm((x - X37) / h)) / (N * h)
  est38 <- sum(stats::dnorm((x - X38) / h)) / (N * h)
  est39 <- sum(stats::dnorm((x - X39) / h)) / (N * h)

  # TOJ estimates
  tojest <- 3.536 * est -
    4.072 * (est21 + est22) / 2 +
    1.536 * (est31 + est32 + est33 + est34 +
               est35 + est36 + est37 + est38 + est39) / 9

  # Ensure non-negative estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)

  return(tojest)

}, vectorize.args = "x")

#' Compute TOJ kernel density estimate for T equivalent to 3 modulo 6
#'
#' @param x An evaluation point
#' @param X A vector of cross-sectional data
#' @param X21 A vector of half-panel cross-sectional data 21
#' @param X22 A vector of half-panel cross-sectional data 22
#' @param X23 A vector of half-panel cross-sectional data 23
#' @param X24 A vector of half-panel cross-sectional data 24
#' @param X31 A vector of third-panel cross-sectional data 31
#' @param X32 A vector of third-panel cross-sectional data 32
#' @param X33 A vector of third-panel cross-sectional data 33
#' @param h A scalar of bandwidth
#'
#' @returns A vector of kernel density estimates
#'
tojkdest3 <- Vectorize(FUN = function(x,
                                      X,
                                      X21,
                                      X22,
                                      X23,
                                      X24,
                                      X31,
                                      X32,
                                      X33,
                                      h) {

  # Sample size
  N <- length(X)

  # Estimates
  est   <- sum(stats::dnorm((x - X)   / h)) / (N * h)
  est21 <- sum(stats::dnorm((x - X21) / h)) / (N * h)
  est22 <- sum(stats::dnorm((x - X22) / h)) / (N * h)
  est23 <- sum(stats::dnorm((x - X23) / h)) / (N * h)
  est24 <- sum(stats::dnorm((x - X24) / h)) / (N * h)
  est31 <- sum(stats::dnorm((x - X31) / h)) / (N * h)
  est32 <- sum(stats::dnorm((x - X32) / h)) / (N * h)
  est33 <- sum(stats::dnorm((x - X33) / h)) / (N * h)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22 + est23 + est24) / 4 +
    1.536 * (est31 + est32 + est33) / 3

  # Ensure non-negative estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)

  return(tojest)

}, vectorize.args = "x")

#' Compute TOJ kernel density estimate for T equivalent to 4 modulo 6
#'
#' @param x An evaluation point
#' @param X A vector of cross-sectional data
#' @param X21 A vector of half-panel cross-sectional data 21
#' @param X22 A vector of half-panel cross-sectional data 22
#' @param X31 A vector of third-panel cross-sectional data 31
#' @param X32 A vector of third-panel cross-sectional data 32
#' @param X33 A vector of third-panel cross-sectional data 33
#' @param X34 A vector of third-panel cross-sectional data 34
#' @param X35 A vector of third-panel cross-sectional data 35
#' @param X36 A vector of third-panel cross-sectional data 36
#' @param X37 A vector of third-panel cross-sectional data 37
#' @param X38 A vector of third-panel cross-sectional data 38
#' @param X39 A vector of third-panel cross-sectional data 39
#' @param h A scalar of bandwidth
#'
#' @returns A vector of kernel density estimates
#'
tojkdest4 <- Vectorize(FUN = function(x,
                                      X,
                                      X21,
                                      X22,
                                      X31,
                                      X32,
                                      X33,
                                      X34,
                                      X35,
                                      X36,
                                      X37,
                                      X38,
                                      X39,
                                      h) {

  # Sample size
  N <- length(X)

  # Estimates
  est   <- sum(stats::dnorm((x - X)   / h)) / (N * h)
  est21 <- sum(stats::dnorm((x - X21) / h)) / (N * h)
  est22 <- sum(stats::dnorm((x - X22) / h)) / (N * h)
  est31 <- sum(stats::dnorm((x - X31) / h)) / (N * h)
  est32 <- sum(stats::dnorm((x - X32) / h)) / (N * h)
  est33 <- sum(stats::dnorm((x - X33) / h)) / (N * h)
  est34 <- sum(stats::dnorm((x - X34) / h)) / (N * h)
  est35 <- sum(stats::dnorm((x - X35) / h)) / (N * h)
  est36 <- sum(stats::dnorm((x - X36) / h)) / (N * h)
  est37 <- sum(stats::dnorm((x - X37) / h)) / (N * h)
  est38 <- sum(stats::dnorm((x - X38) / h)) / (N * h)
  est39 <- sum(stats::dnorm((x - X39) / h)) / (N * h)

  # TOJ estimates
  tojest <- 3.536 * est -
    4.072 * (est21 + est22) / 2 +
    1.536 * (est31 + est32 + est33 + est34 +
               est35 + est36 + est37 + est38 + est39) / 9

  # Ensure non-negative estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)

  return(tojest)

}, vectorize.args = "x")

#' Compute TOJ kernel density estimate for T equivalent to 5 modulo 6
#'
#' @param x An evaluation point
#' @param X A vector of cross-sectional data
#' @param X21 A vector of half-panel cross-sectional data 21
#' @param X22 A vector of half-panel cross-sectional data 22
#' @param X23 A vector of half-panel cross-sectional data 23
#' @param X24 A vector of half-panel cross-sectional data 24
#' @param X31 A vector of third-panel cross-sectional data 31
#' @param X32 A vector of third-panel cross-sectional data 32
#' @param X33 A vector of third-panel cross-sectional data 33
#' @param X34 A vector of third-panel cross-sectional data 34
#' @param X35 A vector of third-panel cross-sectional data 35
#' @param X36 A vector of third-panel cross-sectional data 36
#' @param X37 A vector of third-panel cross-sectional data 37
#' @param X38 A vector of third-panel cross-sectional data 38
#' @param X39 A vector of third-panel cross-sectional data 39
#' @param h A scalar of bandwidth
#'
#' @returns A vector of kernel density estimates
#'
tojkdest5 <- Vectorize(FUN = function(x,
                                      X,
                                      X21,
                                      X22,
                                      X23,
                                      X24,
                                      X31,
                                      X32,
                                      X33,
                                      X34,
                                      X35,
                                      X36,
                                      X37,
                                      X38,
                                      X39,
                                      h) {

  # Sample size
  N <- length(X)

  # Estimates
  est   <- sum(stats::dnorm((x - X)   / h)) / (N * h)
  est21 <- sum(stats::dnorm((x - X21) / h)) / (N * h)
  est22 <- sum(stats::dnorm((x - X22) / h)) / (N * h)
  est23 <- sum(stats::dnorm((x - X23) / h)) / (N * h)
  est24 <- sum(stats::dnorm((x - X24) / h)) / (N * h)
  est31 <- sum(stats::dnorm((x - X31) / h)) / (N * h)
  est32 <- sum(stats::dnorm((x - X32) / h)) / (N * h)
  est33 <- sum(stats::dnorm((x - X33) / h)) / (N * h)
  est34 <- sum(stats::dnorm((x - X34) / h)) / (N * h)
  est35 <- sum(stats::dnorm((x - X35) / h)) / (N * h)
  est36 <- sum(stats::dnorm((x - X36) / h)) / (N * h)
  est37 <- sum(stats::dnorm((x - X37) / h)) / (N * h)
  est38 <- sum(stats::dnorm((x - X38) / h)) / (N * h)
  est39 <- sum(stats::dnorm((x - X39) / h)) / (N * h)

  # TOJ estimates
  tojest <- 3.536 * est -
    4.072 * (est21 + est22 + est23 + est24) / 4 +
    1.536 * (est31 + est32 + est33 + est34 +
               est35 + est36 + est37 + est38 + est39) / 9

  # Ensure non-negative estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)

  return(tojest)

}, vectorize.args = "x")
