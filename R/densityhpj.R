#' HPJ bias-corrected kernel density estimation for heterogeneity in panel data
#'
#' \code{hpjkd} implements the HPJ bias-corrected estimation of
#' the kernel density for the heterogeneous mean, autocovariance,
#' and autocorrelation.
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
#' \item{mean}{graph of the HPJ bias-corrected kernel density estimation for the mean}
#' \item{acov}{graph of the HPJ bias-corrected kernel density estimation for the autocovariance}
#' \item{acor}{graph of the HPJ bias-corrected kernel density estimation for the autocorrelation}
#' \item{mean_func}{function that returns HPJ bias-corrected kernel density estimates for the mean}
#' \item{acov_func}{function that returns HPJ bias-corrected kernel density estimates for the autocovariance}
#' \item{acor_func}{function that returns HPJ bias-corrected kernel density estimates for the autocorrelation}
#' \item{bandwidth}{vector of the selected bandwidths}
#' \item{quantity}{matrix that contains the estimated quantities}
#' \item{acov_order}{the same as the argument}
#' \item{acor_order}{the same as the argument}
#' \item{N}{the number of cross-sectional units}
#' \item{S}{the length of time series}
#'
#' @export
#'
hpjkd <- function(data, acov_order = 0, acor_order = 1, mean_bw = NULL, acov_bw = NULL, acor_bw = NULL) {

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

  # limits for figures by ggplot2
  mean_lim <- c(min(mean_est), max(mean_est))
  acov_lim <- c(min(acov_est), max(acov_est))
  acor_lim <- c(min(acor_est), max(acor_est))

  # HPJ bias-correction
  if (S %% 2 == 0) {

    # half panel data for even T
    data1 <- data[, 1:(S / 2)]
    data2 <- data[, (S / 2 + 1):S]

    # estimated quantities for half panel data
    mean_est1 <- rowMeans(data1)
    mean_est2 <- rowMeans(data2)

    acov_est1 <- apply(data1, MARGIN = 1, acov, acov_order = acov_order)
    acov_est2 <- apply(data2, MARGIN = 1, acov, acov_order = acov_order)

    acor_est1 <- apply(data1, MARGIN = 1, acor, acor_order = acor_order)
    acor_est2 <- apply(data2, MARGIN = 1, acor, acor_order = acor_order)

    # figures by ggplot2
    mean_plot <- ggplot(data = data.frame(x = mean_lim), aes(x = x))
    mean_plot <- mean_plot + stat_function(fun = hpjkdest1, args = list(X = mean_est, X1 = mean_est1, X2 = mean_est2, h = mean_bw))
    mean_plot <- mean_plot + labs(x = "x", y = "")

    acov_plot <- ggplot(data = data.frame(x = acov_lim), aes(x = x))
    acov_plot <- acov_plot + stat_function(fun = hpjkdest1, args = list(X = acov_est, X1 = acov_est1, X2 = acov_est2, h = acov_bw))
    acov_plot <- acov_plot + labs(x = "x", y = "")

    acor_plot <- ggplot(data = data.frame(x = acor_lim), aes(x = x))
    acor_plot <- acor_plot + stat_function(fun = hpjkdest1, args = list(X = acor_est, X1 = acor_est1, X2 = acor_est2, h = acor_bw))
    acor_plot <- acor_plot + labs(x = "x", y = "")

    # functions
    mean_func <- function(x) {
      hpjkdest1(x = x, X = mean_est, X1 = mean_est1, X2 = mean_est2, h = mean_bw)
    }

    acov_func <- function(x) {
      hpjkdest1(x = x, X = acov_est, X1 = acov_est1, X2 = acov_est2, h = acov_bw)
    }

    acor_func <- function(x) {
      hpjkdest1(x = x, X = acor_est, X1 = acor_est1, X2 = acor_est2, h = acor_bw)
    }

  } else {

    # half-panel data for odd T
    data1 <- data[, 1:floor(S / 2)]
    data2 <- data[, (floor(S / 2) + 1):S]
    data3 <- data[, 1:ceiling(S / 2)]
    data4 <- data[, (ceiling(S / 2) + 1):S]

    # estimated quantities for half panel data
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

    # figures by ggplot2
    mean_plot <- ggplot(data = data.frame(x = mean_lim), aes(x = x))
    mean_plot <- mean_plot + stat_function(fun = hpjkdest2, args = list(X = mean_est, X1 = mean_est1, X2 = mean_est2, X3 = mean_est3, X4 = mean_est4, h = mean_bw))
    mean_plot <- mean_plot + labs(x = "x", y = "")

    acov_plot <- ggplot(data = data.frame(x = acov_lim), aes(x = x))
    acov_plot <- acov_plot + stat_function(fun = hpjkdest2, args = list(X = acov_est, X1 = acov_est1, X2 = acov_est2, X3 = acov_est3, X4 = acov_est4, h = acov_bw))
    acov_plot <- acov_plot + labs(x = "x", y = "")

    acor_plot <- ggplot(data = data.frame(x = acor_lim), aes(x = x))
    acor_plot <- acor_plot + stat_function(fun = hpjkdest2, args = list(X = acor_est, X1 = acor_est1, X2 = acor_est2, X3 = acor_est3, X4 = acor_est4, h = acor_bw))
    acor_plot <- acor_plot + labs(x = "x", y = "")

    # functions
    mean_func <- function(x) {
      hpjkdest2(x = x, X = mean_est, X1 = mean_est1, X2 = mean_est2, X3 = mean_est3, X4 = mean_est4, h = mean_bw)
    }

    acov_func <- function(x) {
      hpjkdest2(x = x, X = acov_est, X1 = acov_est1, X2 = acov_est2, X3 = acov_est3, X4 = acov_est4, h = acov_bw)
    }

    acor_func <- function(x) {
      hpjkdest2(x = x, X = acor_est, X1 = acor_est1, X2 = acor_est2, X3 = acor_est3, X4 = acor_est4, h = acor_bw)
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


#' computing HPJ kernel density estimate for even T
#'
#' @param x point at which the density is estimated
#' @param X vector of original cross-sectional data
#' @param X1 vector of half-panel cross-sectional data based on time series 1 ~ T/2
#' @param X2 vector of half-panel cross-sectional data based on time series (T/2 + 1) ~ T
#' @param h bandwidth
#'
hpjkdest1 <- Vectorize(FUN = function(x, X, X1, X2, h) {

  # sample size
  N <- length(X)

  # estimates
  est <- sum(dnorm( (x - X) / h)) /  (N * h)
  est1 <- sum(dnorm( (x - X1) / h)) /  (N * h)
  est2 <- sum(dnorm( (x - X2) / h)) /  (N * h)

  # HPJ estimate
  hpjest <- 2 * est - (est1 + est2) / 2

  # correction to ensure non-negative estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)

  return(hpjest)

}, vectorize.args = "x")


#' computing HPJ kernel density estimate for odd T
#'
#' @param x point at which the density is estimated
#' @param X vector of original cross-sectional data
#' @param X1 vector of half-panel cross-sectional data based on time series 1 ~ floor(T/2)
#' @param X2 vector of half-panel cross-sectional data based on time series (floor(T/2) + 1) ~ T
#' @param X3 vector of half-panel cross-sectional data based on time series 1 ~ ceiling(T/2)
#' @param X4 vector of half-panel cross-sectional data based on time series (ceiling(T/2) + 1) ~ T
#' @param h bandwidth
#'
hpjkdest2 <- Vectorize(FUN = function(x, X, X1, X2, X3, X4, h) {

  # sample size
  N <- length(X)

  # estimates
  est <- sum(dnorm( (x - X) / h)) /  (N * h)
  est1 <- sum(dnorm( (x - X1) / h)) /  (N * h)
  est2 <- sum(dnorm( (x - X2) / h)) /  (N * h)
  est3 <- sum(dnorm( (x - X3) / h)) /  (N * h)
  est4 <- sum(dnorm( (x - X4) / h)) /  (N * h)

  # HPJ estimate
  hpjest <- 2 * est - (est1 + est2 + est3 + est4) / 4

  # correction to ensure non-negative estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)

  return(hpjest)

}, vectorize.args = "x")
