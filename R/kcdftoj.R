#' TOJ bias-corrected kernel CDF estimation for heterogeneity in panel data
#'
#' \code{tojkcdf} implements the TOJ bias-corrected estimation of
#' the kernel density for the heterogeneous mean, autocovariance,
#' and autocorrelation.
#' The procedure is proposed in Okui and Yanagi (2019).
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
#'
#' @importFrom kerdiest PBbw
#' @importFrom Rearrangement rearrangement
#' @importFrom stats na.omit
#'
#' @return list that contains the following elements.
#' \item{mean}{graph of the TOJ bias-corrected kernel CDF estimation for the mean with rearrangement}
#' \item{acov}{graph of the TOJ bias-corrected kernel CDF estimation for the autocovariance with rearrangement}
#' \item{acor}{graph of the TOJ bias-corrected kernel CDF estimation for the autocorrelation with rearrangement}
#' \item{mean_func}{function that returns TOJ bias-corrected kernel CDF estimates for the mean without rearrangement}
#' \item{acov_func}{function that returns TOJ bias-corrected kernel CDF estimates for the autocovariance without rearrangement}
#' \item{acor_func}{function that returns TOJ bias-corrected kernel CDF estimates for the autocorrelation without rearrangement}
#' \item{bandwidth}{vector of the selected bandwidths}
#' \item{quantity}{matrix that contains the estimated quantities}
#' \item{acov_order}{the same as the argument}
#' \item{acor_order}{the same as the argument}
#' \item{N}{the number of cross-sectional units}
#' \item{S}{the length of time series}
#'
#' @export
#'
tojkcdf <- function(data, acov_order = 0, acor_order = 1, mean_bw = NULL, acov_bw = NULL, acor_bw = NULL) {

  # initialization
  x <- y <- NULL

  # handling errors
  error3(data = data, acov_order = acov_order, acor_order = acor_order, mean_bw = mean_bw, acov_bw = acov_bw, acor_bw = acor_bw)

  # omitting NA
  data <- na.omit(data)

  # sample sizes
  N <- nrow(data)
  S <- ncol(data)

  # estimated means, autocovariances, autocovariances
  mean_est <- rowMeans(data)
  acov_est <- apply(data, MARGIN = 1, acov, acov_order = acov_order)
  acor_est <- apply(data, MARGIN = 1, acor, acor_order = acor_order)

  # plug-in bandwidths
  if (is.null(mean_bw)) {
    mean_bw <- PBbw(vec_data = mean_est, type_kernel = "n", num_stage = 2)
  }
  if (is.null(acov_bw)) {
    acov_bw <- PBbw(vec_data = acov_est, type_kernel = "n", num_stage = 2)
  }
  if (is.null(acor_bw)) {
    acor_bw <- PBbw(vec_data = acor_est, type_kernel = "n", num_stage = 2)
  }

  # TOJ bias-correction
  if (S %% 6 == 0) {

    # split  panel data for T equivalent to 0 modulo 6
    data21 <- data[, 1:(S / 2)]
    data22 <- data[, (S / 2 + 1):S]
    data31 <- data[, 1:(S / 3)]
    data32 <- data[, (S / 3 + 1):(2*S / 3)]
    data33 <- data[, (2 * S / 3 + 1):S]

    # estimated quantities for split panel data
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

    # figures based on TOJ estimation with rearrangement
    mean_x <- seq(min(mean_est), max(mean_est), length = 1000)
    mean_y <- tojkcdfest0(x = mean_x, X = mean_est, X21 = mean_est21, X22 = mean_est22, X31 = mean_est31, X32 = mean_est32, X33 = mean_est33, h = mean_bw)
    mean_toj <- rearrangement(x = data.frame(x = mean_x), y = mean_y)
    mean_plot <- ggplot(data.frame(x = mean_x, y = mean_toj), aes(x = x, y = y))
    mean_plot <- mean_plot + geom_line() + xlim(min(mean_est), max(mean_est)) + ylim(0, 1)
    mean_plot <- mean_plot + labs(x = "x", y = "")

    acov_x <- seq(min(acov_est), max(acov_est), length = 1000)
    acov_y <- tojkcdfest0(x = acov_x, X = acov_est, X21 = acov_est21, X22 = acov_est22, X31 = acov_est31, X32 = acov_est32, X33 = acov_est33, h = acov_bw)
    acov_toj <- rearrangement(x = data.frame(x = acov_x), y = acov_y)
    acov_plot <- ggplot(data.frame(x = acov_x, y = acov_toj), aes(x = x, y = y))
    acov_plot <- acov_plot + geom_line() + xlim(min(acov_est), max(acov_est)) + ylim(0, 1)
    acov_plot <- acov_plot + labs(x = "x", y = "")

    acor_x <- seq(min(acor_est), max(acor_est), length = 1000)
    acor_y <- tojkcdfest0(x = acor_x, X = acor_est, X21 = acor_est21, X22 = acor_est22, X31 = acor_est31, X32 = acor_est32, X33 = acor_est33, h = acor_bw)
    acor_toj <- rearrangement(x = data.frame(x = acor_x), y = acor_y)
    acor_plot <- ggplot(data.frame(x = acor_x, y = acor_toj), aes(x = x, y = y))
    acor_plot <- acor_plot + geom_line() + xlim(min(acor_est), max(acor_est)) + ylim(0, 1)
    acor_plot <- acor_plot + labs(x = "x", y = "")

    # functions by TOJ estimation without rearrangement
    mean_func <- function(x) {
      tojkcdfest0(x = x, X = mean_est, X21 = mean_est21, X22 = mean_est22, X31 = mean_est31, X32 = mean_est32, X33 = mean_est33, h = mean_bw)
    }

    acov_func <- function(x) {
      tojkcdfest0(x = x, X = acov_est, X21 = acov_est21, X22 = acov_est22, X31 = acov_est31, X32 = acov_est32, X33 = acov_est33, h = acov_bw)
    }

    acor_func <- function(x) {
      tojkcdfest0(x = x, X = acor_est, X21 = acor_est21, X22 = acor_est22, X31 = acor_est31, X32 = acor_est32, X33 = acor_est33, h = acor_bw)
    }

  } else if (S %% 6 == 1) {

    # split  panel data for T equivalent to 1 modulo 6
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

    # estimated quantities for split panel data
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

    # figures based on TOJ estimation with rearrangement
    mean_x <- seq(min(mean_est), max(mean_est), length = 1000)
    mean_y <- tojkcdfest1(x = mean_x, X = mean_est, X21 = mean_est21, X22 = mean_est22, X23 = mean_est23, X24 = mean_est24, X31 = mean_est31, X32 = mean_est32, X33 = mean_est33, X34 = mean_est34, X35 = mean_est35, X36 = mean_est36, X37 = mean_est37, X38 = mean_est38, X39 = mean_est39, h = mean_bw)
    mean_toj <- rearrangement(x = data.frame(x = mean_x), y = mean_y)
    mean_plot <- ggplot(data.frame(x = mean_x, y = mean_toj), aes(x = x, y = y))
    mean_plot <- mean_plot + geom_line() + xlim(min(mean_est), max(mean_est)) + ylim(0, 1)
    mean_plot <- mean_plot + labs(x = "x", y = "")

    acov_x <- seq(min(acov_est), max(acov_est), length = 1000)
    acov_y <- tojkcdfest1(x = acov_x, X = acov_est, X21 = acov_est21, X22 = acov_est22, X23 = acov_est23, X24 = acov_est24, X31 = acov_est31, X32 = acov_est32, X33 = acov_est33, X34 = acov_est34, X35 = acov_est35, X36 = acov_est36, X37 = acov_est37, X38 = acov_est38, X39 = acov_est39, h = acov_bw)
    acov_toj <- rearrangement(x = data.frame(x = acov_x), y = acov_y)
    acov_plot <- ggplot(data.frame(x = acov_x, y = acov_toj), aes(x = x, y = y))
    acov_plot <- acov_plot + geom_line() + xlim(min(acov_est), max(acov_est)) + ylim(0, 1)
    acov_plot <- acov_plot + labs(x = "x", y = "")

    acor_x <- seq(min(acor_est), max(acor_est), length = 1000)
    acor_y <- tojkcdfest1(x = acor_x, X = acor_est, X21 = acor_est21, X22 = acor_est22, X23 = acor_est23, X24 = acor_est24, X31 = acor_est31, X32 = acor_est32, X33 = acor_est33, X34 = acor_est34, X35 = acor_est35, X36 = acor_est36, X37 = acor_est37, X38 = acor_est38, X39 = acor_est39, h = acor_bw)
    acor_toj <- rearrangement(x = data.frame(x = acor_x), y = acor_y)
    acor_plot <- ggplot(data.frame(x = acor_x, y = acor_toj), aes(x = x, y = y))
    acor_plot <- acor_plot + geom_line() + xlim(min(acor_est), max(acor_est)) + ylim(0, 1)
    acor_plot <- acor_plot + labs(x = "x", y = "")

    # functions by TOJ estimation without rearrangement
    mean_func <- function(x) {
      tojkcdfest1(x = x, X = mean_est, X21 = mean_est21, X22 = mean_est22, X23 = mean_est23, X24 = mean_est24, X31 = mean_est31, X32 = mean_est32, X33 = mean_est33, X34 = mean_est34, X35 = mean_est35, X36 = mean_est36, X37 = mean_est37, X38 = mean_est38, X39 = mean_est39, h = mean_bw)
    }

    acov_func <- function(x) {
      tojkcdfest1(x = x, X = acov_est, X21 = acov_est21, X22 = acov_est22, X23 = acov_est23, X24 = acov_est24, X31 = acov_est31, X32 = acov_est32, X33 = acov_est33, X34 = acov_est34, X35 = acov_est35, X36 = acov_est36, X37 = acov_est37, X38 = acov_est38, X39 = acov_est39, h = acov_bw)
    }

    acor_func <- function(x) {
      tojkcdfest1(x = x, X = acor_est, X21 = acor_est21, X22 = acor_est22, X23 = acor_est23, X24 = acor_est24, X31 = acor_est31, X32 = acor_est32, X33 = acor_est33, X34 = acor_est34, X35 = acor_est35, X36 = acor_est36, X37 = acor_est37, X38 = acor_est38, X39 = acor_est39, h = acor_bw)
    }


  } else if (S %% 6 == 2) {

    # split  panel data for T equivalent to 2 modulo 6
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

    # estimated quantities for split panel data
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

    # figures based on TOJ estimation with rearrangement
    mean_x <- seq(min(mean_est), max(mean_est), length = 1000)
    mean_y <- tojkcdfest2(x = mean_x, X = mean_est, X21 = mean_est21, X22 = mean_est22, X31 = mean_est31, X32 = mean_est32, X33 = mean_est33, X34 = mean_est34, X35 = mean_est35, X36 = mean_est36, X37 = mean_est37, X38 = mean_est38, X39 = mean_est39, h = mean_bw)
    mean_toj <- rearrangement(x = data.frame(x = mean_x), y = mean_y)
    mean_plot <- ggplot(data.frame(x = mean_x, y = mean_toj), aes(x = x, y = y))
    mean_plot <- mean_plot + geom_line() + xlim(min(mean_est), max(mean_est)) + ylim(0, 1)
    mean_plot <- mean_plot + labs(x = "x", y = "")

    acov_x <- seq(min(acov_est), max(acov_est), length = 1000)
    acov_y <- tojkcdfest2(x = acov_x, X = acov_est, X21 = acov_est21, X22 = acov_est22, X31 = acov_est31, X32 = acov_est32, X33 = acov_est33, X34 = acov_est34, X35 = acov_est35, X36 = acov_est36, X37 = acov_est37, X38 = acov_est38, X39 = acov_est39, h = acov_bw)
    acov_toj <- rearrangement(x = data.frame(x = acov_x), y = acov_y)
    acov_plot <- ggplot(data.frame(x = acov_x, y = acov_toj), aes(x = x, y = y))
    acov_plot <- acov_plot + geom_line() + xlim(min(acov_est), max(acov_est)) + ylim(0, 1)
    acov_plot <- acov_plot + labs(x = "x", y = "")

    acor_x <- seq(min(acor_est), max(acor_est), length = 1000)
    acor_y <- tojkcdfest2(x = acor_x, X = acor_est, X21 = acor_est21, X22 = acor_est22, X31 = acor_est31, X32 = acor_est32, X33 = acor_est33, X34 = acor_est34, X35 = acor_est35, X36 = acor_est36, X37 = acor_est37, X38 = acor_est38, X39 = acor_est39, h = acor_bw)
    acor_toj <- rearrangement(x = data.frame(x = acor_x), y = acor_y)
    acor_plot <- ggplot(data.frame(x = acor_x, y = acor_toj), aes(x = x, y = y))
    acor_plot <- acor_plot + geom_line() + xlim(min(acor_est), max(acor_est)) + ylim(0, 1)
    acor_plot <- acor_plot + labs(x = "x", y = "")

    # functions by TOJ estimation without rearrangement
    mean_func <- function(x) {
      tojkcdfest2(x = x, X = mean_est, X21 = mean_est21, X22 = mean_est22, X31 = mean_est31, X32 = mean_est32, X33 = mean_est33, X34 = mean_est34, X35 = mean_est35, X36 = mean_est36, X37 = mean_est37, X38 = mean_est38, X39 = mean_est39, h = mean_bw)
    }

    acov_func <- function(x) {
      tojkcdfest2(x = x, X = acov_est, X21 = acov_est21, X22 = acov_est22, X31 = acov_est31, X32 = acov_est32, X33 = acov_est33, X34 = acov_est34, X35 = acov_est35, X36 = acov_est36, X37 = acov_est37, X38 = acov_est38, X39 = acov_est39, h = acov_bw)
    }

    acor_func <- function(x) {
      tojkcdfest2(x = x, X = acor_est, X21 = acor_est21, X22 = acor_est22, X31 = acor_est31, X32 = acor_est32, X33 = acor_est33, X34 = acor_est34, X35 = acor_est35, X36 = acor_est36, X37 = acor_est37, X38 = acor_est38, X39 = acor_est39, h = acor_bw)
    }


  } else if (S %% 6 == 3) {

    # split  panel data for T equivalent to 3 modulo 6
    data21 <- data[, 1:floor(S / 2)]
    data22 <- data[, (floor(S / 2) + 1):S]
    data23 <- data[, 1:ceiling(S / 2)]
    data24 <- data[, (ceiling(S / 2) + 1):S]
    data31 <- data[, 1:(S / 3)]
    data32 <- data[, (S / 3 + 1):(2*S / 3)]
    data33 <- data[, (2 * S / 3 + 1):S]

    # estimated quantities for split panel data
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

    # figures based on TOJ estimation with rearrangement
    mean_x <- seq(min(mean_est), max(mean_est), length = 1000)
    mean_y <- tojkcdfest3(x = mean_x, X = mean_est, X21 = mean_est21, X22 = mean_est22, X23 = mean_est23, X24 = mean_est24, X31 = mean_est31, X32 = mean_est32, X33 = mean_est33, h = mean_bw)
    mean_toj <- rearrangement(x = data.frame(x = mean_x), y = mean_y)
    mean_plot <- ggplot(data.frame(x = mean_x, y = mean_toj), aes(x = x, y = y))
    mean_plot <- mean_plot + geom_line() + xlim(min(mean_est), max(mean_est)) + ylim(0, 1)
    mean_plot <- mean_plot + labs(x = "x", y = "")

    acov_x <- seq(min(acov_est), max(acov_est), length = 1000)
    acov_y <- tojkcdfest3(x = acov_x, X = acov_est, X21 = acov_est21, X22 = acov_est22, X23 = acov_est23, X24 = acov_est24, X31 = acov_est31, X32 = acov_est32, X33 = acov_est33, h = acov_bw)
    acov_toj <- rearrangement(x = data.frame(x = acov_x), y = acov_y)
    acov_plot <- ggplot(data.frame(x = acov_x, y = acov_toj), aes(x = x, y = y))
    acov_plot <- acov_plot + geom_line() + xlim(min(acov_est), max(acov_est)) + ylim(0, 1)
    acov_plot <- acov_plot + labs(x = "x", y = "")

    acor_x <- seq(min(acor_est), max(acor_est), length = 1000)
    acor_y <- tojkcdfest3(x = acor_x, X = acor_est, X21 = acor_est21, X22 = acor_est22, X23 = acor_est23, X24 = acor_est24, X31 = acor_est31, X32 = acor_est32, X33 = acor_est33, h = acor_bw)
    acor_toj <- rearrangement(x = data.frame(x = acor_x), y = acor_y)
    acor_plot <- ggplot(data.frame(x = acor_x, y = acor_toj), aes(x = x, y = y))
    acor_plot <- acor_plot + geom_line() + xlim(min(acor_est), max(acor_est)) + ylim(0, 1)
    acor_plot <- acor_plot + labs(x = "x", y = "")

    # functions by TOJ estimation without rearrangement
    mean_func <- function(x) {
      tojkcdfest3(x = x, X = mean_est, X21 = mean_est21, X22 = mean_est22, X23 = mean_est23, X24 = mean_est24, X31 = mean_est31, X32 = mean_est32, X33 = mean_est33, h = mean_bw)
    }

    acov_func <- function(x) {
      tojkcdfest3(x = x, X = acov_est, X21 = acov_est21, X22 = acov_est22, X23 = acov_est23, X24 = acov_est24, X31 = acov_est31, X32 = acov_est32, X33 = acov_est33, h = acov_bw)
    }

    acor_func <- function(x) {
      tojkcdfest3(x = x, X = acor_est, X21 = acor_est21, X22 = acor_est22, X23 = acor_est23, X24 = acor_est24, X31 = acor_est31, X32 = acor_est32, X33 = acor_est33, h = acor_bw)
    }


  } else if (S %% 6 == 4) {

    # split  panel data for T equivalent to 4 modulo 6
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

    # estimated quantities for split panel data
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

    # figures based on TOJ estimation with rearrangement
    mean_x <- seq(min(mean_est), max(mean_est), length = 1000)
    mean_y <- tojkcdfest4(x = mean_x, X = mean_est, X21 = mean_est21, X22 = mean_est22, X31 = mean_est31, X32 = mean_est32, X33 = mean_est33, X34 = mean_est34, X35 = mean_est35, X36 = mean_est36, X37 = mean_est37, X38 = mean_est38, X39 = mean_est39, h = mean_bw)
    mean_toj <- rearrangement(x = data.frame(x = mean_x), y = mean_y)
    mean_plot <- ggplot(data.frame(x = mean_x, y = mean_toj), aes(x = x, y = y))
    mean_plot <- mean_plot + geom_line() + xlim(min(mean_est), max(mean_est)) + ylim(0, 1)
    mean_plot <- mean_plot + labs(x = "x", y = "")

    acov_x <- seq(min(acov_est), max(acov_est), length = 1000)
    acov_y <- tojkcdfest4(x = acov_x, X = acov_est, X21 = acov_est21, X22 = acov_est22, X31 = acov_est31, X32 = acov_est32, X33 = acov_est33, X34 = acov_est34, X35 = acov_est35, X36 = acov_est36, X37 = acov_est37, X38 = acov_est38, X39 = acov_est39, h = acov_bw)
    acov_toj <- rearrangement(x = data.frame(x = acov_x), y = acov_y)
    acov_plot <- ggplot(data.frame(x = acov_x, y = acov_toj), aes(x = x, y = y))
    acov_plot <- acov_plot + geom_line() + xlim(min(acov_est), max(acov_est)) + ylim(0, 1)
    acov_plot <- acov_plot + labs(x = "x", y = "")

    acor_x <- seq(min(acor_est), max(acor_est), length = 1000)
    acor_y <- tojkcdfest4(x = acor_x, X = acor_est, X21 = acor_est21, X22 = acor_est22, X31 = acor_est31, X32 = acor_est32, X33 = acor_est33, X34 = acor_est34, X35 = acor_est35, X36 = acor_est36, X37 = acor_est37, X38 = acor_est38, X39 = acor_est39, h = acor_bw)
    acor_toj <- rearrangement(x = data.frame(x = acor_x), y = acor_y)
    acor_plot <- ggplot(data.frame(x = acor_x, y = acor_toj), aes(x = x, y = y))
    acor_plot <- acor_plot + geom_line() + xlim(min(acor_est), max(acor_est)) + ylim(0, 1)
    acor_plot <- acor_plot + labs(x = "x", y = "")

    # functions by TOJ estimation without rearrangement
    mean_func <- function(x) {
      tojkcdfest4(x = x, X = mean_est, X21 = mean_est21, X22 = mean_est22, X31 = mean_est31, X32 = mean_est32, X33 = mean_est33, X34 = mean_est34, X35 = mean_est35, X36 = mean_est36, X37 = mean_est37, X38 = mean_est38, X39 = mean_est39, h = mean_bw)
    }

    acov_func <- function(x) {
      tojkcdfest4(x = x, X = acov_est, X21 = acov_est21, X22 = acov_est22, X31 = acov_est31, X32 = acov_est32, X33 = acov_est33, X34 = acov_est34, X35 = acov_est35, X36 = acov_est36, X37 = acov_est37, X38 = acov_est38, X39 = acov_est39, h = acov_bw)
    }

    acor_func <- function(x) {
      tojkcdfest4(x = x, X = acor_est, X21 = acor_est21, X22 = acor_est22, X31 = acor_est31, X32 = acor_est32, X33 = acor_est33, X34 = acor_est34, X35 = acor_est35, X36 = acor_est36, X37 = acor_est37, X38 = acor_est38, X39 = acor_est39, h = acor_bw)
    }


  } else {

    # split  panel data for T equivalent to 5 modulo 6
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

    # estimated quantities for split panel data
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

    # figures based on TOJ estimation with rearrangement
    mean_x <- seq(min(mean_est), max(mean_est), length = 1000)
    mean_y <- tojkcdfest5(x = mean_x, X = mean_est, X21 = mean_est21, X22 = mean_est22, X23 = mean_est23, X24 = mean_est24, X31 = mean_est31, X32 = mean_est32, X33 = mean_est33, X34 = mean_est34, X35 = mean_est35, X36 = mean_est36, X37 = mean_est37, X38 = mean_est38, X39 = mean_est39, h = mean_bw)
    mean_toj <- rearrangement(x = data.frame(x = mean_x), y = mean_y)
    mean_plot <- ggplot(data.frame(x = mean_x, y = mean_toj), aes(x = x, y = y))
    mean_plot <- mean_plot + geom_line() + xlim(min(mean_est), max(mean_est)) + ylim(0, 1)
    mean_plot <- mean_plot + labs(x = "x", y = "")

    acov_x <- seq(min(acov_est), max(acov_est), length = 1000)
    acov_y <- tojkcdfest5(x = acov_x, X = acov_est, X21 = acov_est21, X22 = acov_est22, X23 = acov_est23, X24 = acov_est24, X31 = acov_est31, X32 = acov_est32, X33 = acov_est33, X34 = acov_est34, X35 = acov_est35, X36 = acov_est36, X37 = acov_est37, X38 = acov_est38, X39 = acov_est39, h = acov_bw)
    acov_toj <- rearrangement(x = data.frame(x = acov_x), y = acov_y)
    acov_plot <- ggplot(data.frame(x = acov_x, y = acov_toj), aes(x = x, y = y))
    acov_plot <- acov_plot + geom_line() + xlim(min(acov_est), max(acov_est)) + ylim(0, 1)
    acov_plot <- acov_plot + labs(x = "x", y = "")

    acor_x <- seq(min(acor_est), max(acor_est), length = 1000)
    acor_y <- tojkcdfest5(x = acor_x, X = acor_est, X21 = acor_est21, X22 = acor_est22, X23 = acor_est23, X24 = acor_est24, X31 = acor_est31, X32 = acor_est32, X33 = acor_est33, X34 = acor_est34, X35 = acor_est35, X36 = acor_est36, X37 = acor_est37, X38 = acor_est38, X39 = acor_est39, h = acor_bw)
    acor_toj <- rearrangement(x = data.frame(x = acor_x), y = acor_y)
    acor_plot <- ggplot(data.frame(x = acor_x, y = acor_toj), aes(x = x, y = y))
    acor_plot <- acor_plot + geom_line() + xlim(min(acor_est), max(acor_est)) + ylim(0, 1)
    acor_plot <- acor_plot + labs(x = "x", y = "")

    # functions by TOJ estimation without rearrangement
    mean_func <- function(x) {
      tojkcdfest5(x = x, X = mean_est, X21 = mean_est21, X22 = mean_est22, X23 = mean_est23, X24 = mean_est24, X31 = mean_est31, X32 = mean_est32, X33 = mean_est33, X34 = mean_est34, X35 = mean_est35, X36 = mean_est36, X37 = mean_est37, X38 = mean_est38, X39 = mean_est39, h = mean_bw)
    }

    acov_func <- function(x) {
      tojkcdfest5(x = x, X = acov_est, X21 = acov_est21, X22 = acov_est22, X23 = acov_est23, X24 = acov_est24, X31 = acov_est31, X32 = acov_est32, X33 = acov_est33, X34 = acov_est34, X35 = acov_est35, X36 = acov_est36, X37 = acov_est37, X38 = acov_est38, X39 = acov_est39, h = acov_bw)
    }

    acor_func <- function(x) {
      tojkcdfest5(x = x, X = acor_est, X21 = acor_est21, X22 = acor_est22, X23 = acor_est23, X24 = acor_est24, X31 = acor_est31, X32 = acor_est32, X33 = acor_est33, X34 = acor_est34, X35 = acor_est35, X36 = acor_est36, X37 = acor_est37, X38 = acor_est38, X39 = acor_est39, h = acor_bw)
    }


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


#' computing TOJ kernel CDF estimate for T equivalent to 0 modulo 6
#'
#' @param x point at which the CDF is estimated
#' @param X vector of original cross-sectional data
#' @param X21 vector of half-panel cross-sectional data based on time series 1 ~ T/2
#' @param X22 vector of half-panel cross-sectional data based on time series (T/2 + 1) ~ T
#' @param X31 vector of one-third-panel cross-sectional data based on time series 1 ~ T/3
#' @param X32 vector of one-third-panel cross-sectional data based on time series (T/3 + 1) ~ 2 * T/3
#' @param X33 vector of one-third-panel cross-sectional data based on time series 2 * T/3 + 1 ~ T
#' @param h bandwidth
#'
tojkcdfest0 <- Vectorize(FUN = function(x, X, X21, X22, X31, X32, X33, h) {

  # estimates
  est <- mean(pnorm( (x - X) / h))
  est21 <- mean(pnorm( (x - X21) / h))
  est22 <- mean(pnorm( (x - X22) / h))
  est31 <- mean(pnorm( (x - X31) / h))
  est32 <- mean(pnorm( (x - X32) / h))
  est33 <- mean(pnorm( (x - X33) / h))

  # TOJ estimate
  tojest <- 3.536 * est - 4.072 * (est21 + est22) / 2 + 1.536 * (est31 + est32 + est33) / 3

  # correction to ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}, vectorize.args = "x")


#' computing TOJ kernel CDF estimate for T equivalent to 1 modulo 6
#'
#' @param x point at which the CDF is estimated
#' @param X vector of original cross-sectional data
#' @param X21 vector of half-panel cross-sectional data based on time series 1 ~ floor(T/2)
#' @param X22 vector of half-panel cross-sectional data based on time series (floor(T/2) + 1) ~ T
#' @param X23 vector of half-panel cross-sectional data based on time series 1 ~ ceiling(T/2)
#' @param X24 vector of half-panel cross-sectional data based on time series (ceiling(T/2) + 1) ~ T
#' @param X31 vector of one-third-panel cross-sectional data based on time series 1 ~ floor(T/3)
#' @param X32 vector of one-third-panel cross-sectional data based on time series (floor(T/3) + 1) ~ (2 * floor(T/3))
#' @param X33 vector of one-third-panel cross-sectional data based on time series (2 * floor(T/3) + 1) ~ T
#' @param X34 vector of one-third-panel cross-sectional data based on time series 1 ~ floor(T/3)
#' @param X35 vector of one-third-panel cross-sectional data based on time series (floor(T/3) + 1) ~ (2 * floor(T/3) + 1)
#' @param X36 vector of one-third-panel cross-sectional data based on time series (2 * floor(T/3) + 2) ~ T
#' @param X37 vector of one-third-panel cross-sectional data based on time series 1 ~ ceiling(T/3)
#' @param X38 vector of one-third-panel cross-sectional data based on time series (ceiling(T/3) + 1) ~ (2 * floor(T/3) + 1)
#' @param X39 vector of one-third-panel cross-sectional data based on time series (2 * floor(T/3) + 2) ~ T
#' @param h bandwidth
#'
tojkcdfest1 <- Vectorize(FUN = function(x, X, X21, X22, X23, X24, X31, X32, X33, X34, X35, X36, X37, X38, X39, h) {

  # estimates
  est <- mean(pnorm( (x - X) / h))
  est21 <- mean(pnorm( (x - X21) / h))
  est22 <- mean(pnorm( (x - X22) / h))
  est23 <- mean(pnorm( (x - X23) / h))
  est24 <- mean(pnorm( (x - X24) / h))
  est31 <- mean(pnorm( (x - X31) / h))
  est32 <- mean(pnorm( (x - X32) / h))
  est33 <- mean(pnorm( (x - X33) / h))
  est34 <- mean(pnorm( (x - X34) / h))
  est35 <- mean(pnorm( (x - X35) / h))
  est36 <- mean(pnorm( (x - X36) / h))
  est37 <- mean(pnorm( (x - X37) / h))
  est38 <- mean(pnorm( (x - X38) / h))
  est39 <- mean(pnorm( (x - X39) / h))

  # TOJ estimate
  tojest <- 3.536 * est - 4.072 * (est21 + est22 + est23 + est24) / 4 + 1.536 * (est31 + est32 + est33 + est34 + est35 + est36 + est37 + est38 + est39) / 9

  # correction to ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}, vectorize.args = "x")

#' computing TOJ kernel CDF estimate for T equivalent to 2 modulo 6
#'
#' @param x point at which the CDF is estimated
#' @param X vector of original cross-sectional data
#' @param X21 vector of half-panel cross-sectional data based on time series 1 ~ T/2
#' @param X22 vector of half-panel cross-sectional data based on time series (T/2 + 1) ~ T
#' @param X31 vector of one-third-panel cross-sectional data based on time series 1 ~ floor(T/3)
#' @param X32 vector of one-third-panel cross-sectional data based on time series (floor(T/3) + 1) ~ (2 * floor(T/3) + 1)
#' @param X33 vector of one-third-panel cross-sectional data based on time series (2 * ceiling(T/3)) ~ T
#' @param X34 vector of one-third-panel cross-sectional data based on time series 1 ~ ceiling(T/3)
#' @param X35 vector of one-third-panel cross-sectional data based on time series (ceiling(T/3) + 1) ~ (2 * floor(T/3) + 1)
#' @param X36 vector of one-third-panel cross-sectional data based on time series (2 * ceiling(T/3)) ~ T
#' @param X37 vector of one-third-panel cross-sectional data based on time series 1 ~ ceiling(T/3)
#' @param X38 vector of one-third-panel cross-sectional data based on time series (ceiling(T/3) + 1) ~ (2 * ceiling(T/3))
#' @param X39 vector of one-third-panel cross-sectional data based on time series (2 * ceiling(T/3) + 1) ~ T
#' @param h bandwidth
#'
tojkcdfest2 <- Vectorize(FUN = function(x, X, X21, X22, X31, X32, X33, X34, X35, X36, X37, X38, X39, h) {

  # estimates
  est <- mean(pnorm( (x - X) / h))
  est21 <- mean(pnorm( (x - X21) / h))
  est22 <- mean(pnorm( (x - X22) / h))
  est31 <- mean(pnorm( (x - X31) / h))
  est32 <- mean(pnorm( (x - X32) / h))
  est33 <- mean(pnorm( (x - X33) / h))
  est34 <- mean(pnorm( (x - X34) / h))
  est35 <- mean(pnorm( (x - X35) / h))
  est36 <- mean(pnorm( (x - X36) / h))
  est37 <- mean(pnorm( (x - X37) / h))
  est38 <- mean(pnorm( (x - X38) / h))
  est39 <- mean(pnorm( (x - X39) / h))

  # TOJ estimate
  tojest <- 3.536 * est - 4.072 * (est21 + est22) / 2 + 1.536 * (est31 + est32 + est33 + est34 + est35 + est36 + est37 + est38 + est39) / 9

  # correction to ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}, vectorize.args = "x")

#' computing TOJ kernel CDF estimate for T equivalent to 3 modulo 6
#'
#' @param x point at which the CDF is estimated
#' @param X vector of original cross-sectional data
#' @param X21 vector of half-panel cross-sectional data based on time series 1 ~ floor(T/2)
#' @param X22 vector of half-panel cross-sectional data based on time series (floor(T/2) + 1) ~ T
#' @param X23 vector of half-panel cross-sectional data based on time series 1 ~ ceiling(T/2)
#' @param X24 vector of half-panel cross-sectional data based on time series (ceiling(T/2) + 1) ~ T
#' @param X31 vector of one-third-panel cross-sectional data based on time series 1 ~ T/3
#' @param X32 vector of one-third-panel cross-sectional data based on time series (T/3 + 1) ~ 2 * T/3
#' @param X33 vector of one-third-panel cross-sectional data based on time series 2 * T/3 + 1 ~ T
#' @param h bandwidth
#'
tojkcdfest3 <- Vectorize(FUN = function(x, X, X21, X22, X23, X24, X31, X32, X33, h) {

  # estimates
  est <- mean(pnorm( (x - X) / h))
  est21 <- mean(pnorm( (x - X21) / h))
  est22 <- mean(pnorm( (x - X22) / h))
  est23 <- mean(pnorm( (x - X23) / h))
  est24 <- mean(pnorm( (x - X24) / h))
  est31 <- mean(pnorm( (x - X31) / h))
  est32 <- mean(pnorm( (x - X32) / h))
  est33 <- mean(pnorm( (x - X33) / h))

  # TOJ estimate
  tojest <- 3.536 * est - 4.072 * (est21 + est22 + est23 + est24) / 4 + 1.536 * (est31 + est32 + est33) / 3

  # correction to ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}, vectorize.args = "x")

#' computing TOJ kernel CDF estimate for T equivalent to 4 modulo 6
#'
#' @param x point at which the CDF is estimated
#' @param X vector of original cross-sectional data
#' @param X21 vector of half-panel cross-sectional data based on time series 1 ~ T/2
#' @param X22 vector of half-panel cross-sectional data based on time series (T/2 + 1) ~ T
#' @param X31 vector of one-third-panel cross-sectional data based on time series 1 ~ floor(T/3)
#' @param X32 vector of one-third-panel cross-sectional data based on time series (floor(T/3) + 1) ~ (2 * floor(T/3))
#' @param X33 vector of one-third-panel cross-sectional data based on time series (2 * floor(T/3) + 1) ~ T
#' @param X34 vector of one-third-panel cross-sectional data based on time series 1 ~ floor(T/3)
#' @param X35 vector of one-third-panel cross-sectional data based on time series (floor(T/3) + 1) ~ (2 * floor(T/3) + 1)
#' @param X36 vector of one-third-panel cross-sectional data based on time series (2 * floor(T/3) + 2) ~ T
#' @param X37 vector of one-third-panel cross-sectional data based on time series 1 ~ ceiling(T/3)
#' @param X38 vector of one-third-panel cross-sectional data based on time series (ceiling(T/3) + 1) ~ (2 * floor(T/3) + 1)
#' @param X39 vector of one-third-panel cross-sectional data based on time series (2 * floor(T/3) + 2) ~ T
#' @param h bandwidth
#'
tojkcdfest4 <- Vectorize(FUN = function(x, X, X21, X22, X31, X32, X33, X34, X35, X36, X37, X38, X39, h) {

  # estimates
  est <- mean(pnorm( (x - X) / h))
  est21 <- mean(pnorm( (x - X21) / h))
  est22 <- mean(pnorm( (x - X22) / h))
  est31 <- mean(pnorm( (x - X31) / h))
  est32 <- mean(pnorm( (x - X32) / h))
  est33 <- mean(pnorm( (x - X33) / h))
  est34 <- mean(pnorm( (x - X34) / h))
  est35 <- mean(pnorm( (x - X35) / h))
  est36 <- mean(pnorm( (x - X36) / h))
  est37 <- mean(pnorm( (x - X37) / h))
  est38 <- mean(pnorm( (x - X38) / h))
  est39 <- mean(pnorm( (x - X39) / h))

  # TOJ estimate
  tojest <- 3.536 * est - 4.072 * (est21 + est22) / 2 + 1.536 * (est31 + est32 + est33 + est34 + est35 + est36 + est37 + est38 + est39) / 9

  # correction to ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}, vectorize.args = "x")

#' computing TOJ kernel CDF estimate for T equivalent to 5 modulo 6
#'
#' @param x point at which the CDF is estimated
#' @param X vector of original cross-sectional data
#' @param X21 vector of half-panel cross-sectional data based on time series 1 ~ floor(T/2)
#' @param X22 vector of half-panel cross-sectional data based on time series (floor(T/2) + 1) ~ T
#' @param X23 vector of half-panel cross-sectional data based on time series 1 ~ ceiling(T/2)
#' @param X24 vector of half-panel cross-sectional data based on time series (ceiling(T/2) + 1) ~ T
#' @param X31 vector of one-third-panel cross-sectional data based on time series 1 ~ floor(T/3)
#' @param X32 vector of one-third-panel cross-sectional data based on time series (floor(T/3) + 1) ~ (2 * floor(T/3) + 1)
#' @param X33 vector of one-third-panel cross-sectional data based on time series (2 * ceiling(T/3)) ~ T
#' @param X34 vector of one-third-panel cross-sectional data based on time series 1 ~ ceiling(T/3)
#' @param X35 vector of one-third-panel cross-sectional data based on time series (ceiling(T/3) + 1) ~ (2 * floor(T/3) + 1)
#' @param X36 vector of one-third-panel cross-sectional data based on time series (2 * ceiling(T/3)) ~ T
#' @param X37 vector of one-third-panel cross-sectional data based on time series 1 ~ ceiling(T/3)
#' @param X38 vector of one-third-panel cross-sectional data based on time series (ceiling(T/3) + 1) ~ (2 * ceiling(T/3))
#' @param X39 vector of one-third-panel cross-sectional data based on time series (2 * ceiling(T/3) + 1) ~ T
#' @param h bandwidth
#'
tojkcdfest5 <- Vectorize(FUN = function(x, X, X21, X22, X23, X24, X31, X32, X33, X34, X35, X36, X37, X38, X39, h) {

  # estimates
  est <- mean(pnorm( (x - X) / h))
  est21 <- mean(pnorm( (x - X21) / h))
  est22 <- mean(pnorm( (x - X22) / h))
  est23 <- mean(pnorm( (x - X23) / h))
  est24 <- mean(pnorm( (x - X24) / h))
  est31 <- mean(pnorm( (x - X31) / h))
  est32 <- mean(pnorm( (x - X32) / h))
  est33 <- mean(pnorm( (x - X33) / h))
  est34 <- mean(pnorm( (x - X34) / h))
  est35 <- mean(pnorm( (x - X35) / h))
  est36 <- mean(pnorm( (x - X36) / h))
  est37 <- mean(pnorm( (x - X37) / h))
  est38 <- mean(pnorm( (x - X38) / h))
  est39 <- mean(pnorm( (x - X39) / h))

  # TOJ estimate
  tojest <- 3.536 * est - 4.072 * (est21 + est22 + est23 + est24) / 4 + 1.536 * (est31 + est32 + est33 + est34 + est35 + est36 + est37 + est38 + est39) / 9

  # correction to ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}, vectorize.args = "x")