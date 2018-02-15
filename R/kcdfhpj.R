#' HPJ bias-corrected kernel CDF estimation for heterogeneity in panel data
#'
#' \code{hpjkcdf} implements the HPJ bias-corrected estimation of
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
#'
#' @importFrom kerdiest PBbw
#' @importFrom Rearrangement rearrangement
#' @importFrom stats na.omit
#'
#' @return list that contains the following elements.
#' \item{mean}{graph of the HPJ bias-corrected kernel CDF estimation for the mean with rearrangement}
#' \item{acov}{graph of the HPJ bias-corrected kernel CDF estimation for the autocovariance with rearrangement}
#' \item{acor}{graph of the HPJ bias-corrected kernel CDF estimation for the autocorrelation with rearrangement}
#' \item{mean_func}{function that returns HPJ bias-corrected kernel CDF estimates for the mean without rearrangement}
#' \item{acov_func}{function that returns HPJ bias-corrected kernel CDF estimates for the autocovariance without rearrangement}
#' \item{acor_func}{function that returns HPJ bias-corrected kernel CDF estimates for the autocorrelation without rearrangement}
#' \item{bandwidth}{vector of the selected bandwidths}
#' \item{quantity}{matrix that contains the estimated quantities}
#' \item{acov_order}{the same as the argument}
#' \item{acor_order}{the same as the argument}
#' \item{N}{the number of cross-sectional units}
#' \item{S}{the length of time series}
#'
#' @export
#'
hpjkcdf <- function(data, acov_order = 0, acor_order = 1, mean_bw = NULL, acov_bw = NULL, acor_bw = NULL) {

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

  # HPJ bias-correction
  if (S %% 2 == 0) {

    # half-panel data for even T
    data1 <- data[, 1:(S / 2)]
    data2 <- data[, (S / 2 + 1):S]

    # estimated quantities for half-panel data
    mean_est1 <- rowMeans(data1)
    mean_est2 <- rowMeans(data2)

    acov_est1 <- apply(data1, MARGIN = 1, acov, acov_order = acov_order)
    acov_est2 <- apply(data2, MARGIN = 1, acov, acov_order = acov_order)

    acor_est1 <- apply(data1, MARGIN = 1, acor, acor_order = acor_order)
    acor_est2 <- apply(data2, MARGIN = 1, acor, acor_order = acor_order)

    # figures based on HPJ estimation with rearrangement
    mean_x <- seq(min(mean_est), max(mean_est), length = 1000)
    mean_y <- hpjkcdfest1(x = mean_x, X = mean_est, X1 = mean_est1, X2 = mean_est2, h = mean_bw)
    mean_hpj <- rearrangement(x = data.frame(x = mean_x), y = mean_y)
    mean_plot <- ggplot(data.frame(x = mean_x, y = mean_hpj), aes(x = x, y = y))
    mean_plot <- mean_plot + geom_line() + xlim(min(mean_est), max(mean_est)) + ylim(0, 1)
    mean_plot <- mean_plot + labs(x = "x", y = "")

    acov_x <- seq(min(acov_est), max(acov_est), length = 1000)
    acov_y <- hpjkcdfest1(x = acov_x, X = acov_est, X1 = acov_est1, X2 = acov_est2, h = acov_bw)
    acov_hpj <- rearrangement(x = data.frame(x = acov_x), y = acov_y)
    acov_plot <- ggplot(data.frame(x = acov_x, y = acov_hpj), aes(x = x, y = y))
    acov_plot <- acov_plot + geom_line() + xlim(min(acov_est), max(acov_est)) + ylim(0, 1)
    acov_plot <- acov_plot + labs(x = "x", y = "")

    acor_x <- seq(min(acor_est), max(acor_est), length = 1000)
    acor_y <- hpjkcdfest1(x = acor_x, X = acor_est, X1 = acor_est1, X2 = acor_est2, h = acor_bw)
    acor_hpj <- rearrangement(x = data.frame(x = acor_x), y = acor_y)
    acor_plot <- ggplot(data.frame(x = acor_x, y = acor_hpj), aes(x = x, y = y))
    acor_plot <- acor_plot + geom_line() + xlim(min(acor_est), max(acor_est)) + ylim(0, 1)
    acor_plot <- acor_plot + labs(x = "x", y = "")

    # functions by HPJ estimation without rearrangement
    mean_func <- function(x) {
      hpjkcdfest1(x = x, X = mean_est, X1 = mean_est1, X2 = mean_est2, h = mean_bw)
    }

    acov_func <- function(x) {
      hpjkcdfest1(x = x, X = acov_est, X1 = acov_est1, X2 = acov_est2, h = acov_bw)
    }

    acor_func <- function(x) {
      hpjkcdfest1(x = x, X = acor_est, X1 = acor_est1, X2 = acor_est2, h = acor_bw)
    }

  } else {

    # half-panel data for odd T
    data1 <- data[, 1:floor(S / 2)]
    data2 <- data[, (floor(S / 2) + 1):S]
    data3 <- data[, 1:ceiling(S / 2)]
    data4 <- data[, (ceiling(S / 2) + 1):S]

    # estimated quantities for half-panel data
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

    # figures based on HPJ estimation with rearrangement
    mean_x <- seq(min(mean_est), max(mean_est), length = 1000)
    mean_y <- hpjkcdfest2(x = mean_x, X = mean_est, X1 = mean_est1, X2 = mean_est2, X3 = mean_est3, X4 = mean_est4, h = mean_bw)
    mean_hpj <- rearrangement(x = data.frame(x = mean_x), y = mean_y)
    mean_plot <- ggplot(data.frame(x = mean_x, y = mean_hpj), aes(x = x, y = y))
    mean_plot <- mean_plot + geom_line() + xlim(min(mean_est), max(mean_est)) + ylim(0, 1)
    mean_plot <- mean_plot + labs(x = "x", y = "")

    acov_x <- seq(min(acov_est), max(acov_est), length = 1000)
    acov_y <- hpjkcdfest2(x = acov_x, X = acov_est, X1 = acov_est1, X2 = acov_est2, X3 = acov_est3, X4 = acov_est4, h = acov_bw)
    acov_hpj <- rearrangement(x = data.frame(x = acov_x), y = acov_y)
    acov_plot <- ggplot(data.frame(x = acov_x, y = acov_hpj), aes(x = x, y = y))
    acov_plot <- acov_plot + geom_line() + xlim(min(acov_est), max(acov_est)) + ylim(0, 1)
    acov_plot <- acov_plot + labs(x = "x", y = "")

    acor_x <- seq(min(acor_est), max(acor_est), length = 1000)
    acor_y <- hpjkcdfest2(x = acor_x, X = acor_est, X1 = acor_est1, X2 = acor_est2, X3 = acor_est3, X4 = acor_est4, h = acor_bw)
    acor_hpj <- rearrangement(x = data.frame(x = acor_x), y = acor_y)
    acor_plot <- ggplot(data.frame(x = acor_x, y = acor_hpj), aes(x = x, y = y))
    acor_plot <- acor_plot + geom_line() + xlim(min(acor_est), max(acor_est)) + ylim(0, 1)
    acor_plot <- acor_plot + labs(x = "x", y = "")

    # functions by HPJ estimation without rearrangement
    mean_func <- function(x) {
      hpjkcdfest2(x = x, X = mean_est, X1 = mean_est1, X2 = mean_est2, X3 = mean_est3, X4 = mean_est4, h = mean_bw)
    }

    acov_func <- function(x) {
      hpjkcdfest2(x = x, X = acov_est, X1 = acov_est1, X2 = acov_est2, X3 = acov_est3, X4 = acov_est4, h = acov_bw)
    }

    acor_func <- function(x) {
      hpjkcdfest2(x = x, X = acor_est, X1 = acor_est1, X2 = acor_est2, X3 = acor_est3, X4 = acor_est4, h = acor_bw)
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


#' computing HPJ kernel CDF estimate for even T
#'
#' @param x point at which the CDF is estimated
#' @param X vector of original cross-sectional data
#' @param X1 vector of half-panel cross-sectional data based on time series 1 ~ T/2
#' @param X2 vector of half-panel cross-sectional data based on time series (T/2 + 1) ~ T
#' @param h bandwidth
#'
hpjkcdfest1 <- Vectorize(FUN = function(x, X, X1, X2, h) {

  # estimates
  est <- mean(pnorm( (x - X) / h))
  est1 <- mean(pnorm( (x - X1) / h))
  est2 <- mean(pnorm( (x - X2) / h))

  # HPJ estimate
  hpjest <- 2 * est - (est1 + est2) / 2

  # correction to ensure valid estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)
  hpjest <- ifelse(hpjest <= 1, hpjest, 1)

  return(hpjest)

}, vectorize.args = "x")


#' computing HPJ kernel CDF estimate for odd T
#'
#' @param x point at which the CDF is estimated
#' @param X vector of original cross-sectional data
#' @param X1 vector of half-panel cross-sectional data based on time series 1 ~ floor(T/2)
#' @param X2 vector of half-panel cross-sectional data based on time series (floor(T/2) + 1) ~ T
#' @param X3 vector of half-panel cross-sectional data based on time series 1 ~ ceiling(T/2)
#' @param X4 vector of half-panel cross-sectional data based on time series (ceiling(T/2) + 1) ~ T
#' @param h bandwidth
#'
hpjkcdfest2 <- Vectorize(FUN = function(x, X, X1, X2, X3, X4, h) {

  # estimates
  est <- mean(pnorm( (x - X) / h))
  est1 <- mean(pnorm( (x - X1) / h))
  est2 <- mean(pnorm( (x - X2) / h))
  est3 <- mean(pnorm( (x - X3) / h))
  est4 <- mean(pnorm( (x - X4) / h))

  # HPJ estimate
  hpjest <- 2 * est - (est1 + est2 + est3 + est4) / 4

  # correction to ensure valid estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)
  hpjest <- ifelse(hpjest <= 1, hpjest, 1)

  return(hpjest)

}, vectorize.args = "x")
