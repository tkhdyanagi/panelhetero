#' HPJ bias-corrected empirical CDF estimation for heterogeneity in panel data
#'
#' \code{hpjecdf} implements the HPJ bias-corrected estimation of
#' the empirical CDF for the heterogeneous mean, autocovariance,
#' and autocorrelation.
#' The procedure is proposed in Okui and Yanagi (2019).
#' See the package vignette via \code{vignette("panelhetero")} for details.
#'
#' @param data matrix of panel data in which each row is individual time series
#' @param acov_order non-negative integer for the order of the autocovariance
#' @param acor_order positive integer for the order of the autocorrelation
#' @param ci bool for the confidence interval
#' @param R positive integer for the number of bootstrap replications
#'
#' @import ggplot2
#' @importFrom Rearrangement rearrangement
#' @importFrom stats na.omit
#'
#' @return list that contains the following elements.
#' \item{mean}{graph of the HPJ bias-corrected empirical CDF estimation for the mean with rearrangement}
#' \item{acov}{graph of the HPJ bias-corrected empirical CDF estimation for the autocovariance with rearrangement}
#' \item{acor}{graph of the HPJ bias-corrected empirical CDF estimation for the autocorrelation with rearrangement}
#' \item{mean_func}{function that returns HPJ bias-corrected empirical CDF estimates for the mean without rearrangement}
#' \item{acov_func}{function that returns HPJ bias-corrected empirical CDF estimates for the autocovariance without rearrangement}
#' \item{acor_func}{function that returns HPJ bias-corrected empirical CDF estimates for the autocorrelation without rearrangement}
#' \item{mean_ci_func}{function that returns 95 percent bootstrap confidence interval for HPJ empirical CDF estimates for the mean}
#' \item{acov_ci_func}{function that returns 95 percent bootstrap confidence interval for HPJ empirical CDF estimates for the autocovariance}
#' \item{acor_ci_func}{function that returns 95 percent bootstrap confidence interval for HPJ empirical CDF estimates for the autocorrelation}
#' \item{quantity}{matrix that contains the estimated quantities}
#' \item{acov_order}{the same as the argument}
#' \item{acor_order}{the same as the argument}
#' \item{N}{the number of cross-sectional units}
#' \item{S}{the length of time series}
#' \item{R}{the number of bootstrap replications}
#'
#' @export
#'
hpjecdf <- function(data, acov_order = 0, acor_order = 1, R = 1000, ci = TRUE) {

  # initialization
  x <- y <- NULL

  # handling errors
  error1(data = data, acov_order = acov_order, acor_order = acor_order, R = R)

  # omitting NA
  data <- na.omit(data)

  # sample sizes
  N <- nrow(data)
  S <- ncol(data)

  # estimated means, autocovariances, autocorrelations
  mean_est <- rowMeans(data)
  acov_est <- apply(data, MARGIN = 1, acov, acov_order = acov_order)
  acor_est <- apply(data, MARGIN = 1, acor, acor_order = acor_order)

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

    # function for bootstrap confidence interval
    mean_ci_func <- Vectorize(FUN = function(x){

    # hpj estimation with bootstrap
    bootstrap <- boot(data = cbind(mean_est, mean_est1, mean_est2) - x, statistic = hpj1_boot, R = R)
    estimate <- bootstrap$t0
    temp <- bootstrap$t- bootstrap$t0

    # confidence interval
    quantile <- apply(temp, MARGIN = 2, quantile, probs = c(0.025, 0.975))
    l_ci <- max(0, estimate + quantile[1])
    u_ci <- min(1, estimate + quantile[2])
    ci <- cbind(l_ci, u_ci)
    colnames(ci) <- c("95% CI lower", "95% CI upper")
    ci

    }, vectorize.args = "x")

    acov_ci_func <- Vectorize(FUN = function(x){

    # hpj estimation with bootstrap
    bootstrap <- boot(data = cbind(acov_est, acov_est1, acov_est2) - x, statistic = hpj1_boot, R = R)
    estimate <- bootstrap$t0
    temp <- bootstrap$t- bootstrap$t0

    # confidence interval
    quantile <- apply(temp, MARGIN = 2, quantile, probs = c(0.025, 0.975))
    l_ci <- max(0, estimate + quantile[1])
    u_ci <- min(1, estimate + quantile[2])

    ci <- cbind(l_ci, u_ci)
    colnames(ci) <- c("95% CI lower", "95% CI upper")
    ci

    }, vectorize.args = "x")

    acor_ci_func <- Vectorize(FUN = function(x){

    # hpj estimation with bootstrap
    bootstrap <- boot(data = cbind(acor_est, acor_est1, acor_est2) - x, statistic = hpj1_boot, R = R)
    estimate <- bootstrap$t0
    temp <- bootstrap$t- bootstrap$t0

    # confidence interval
    quantile <- apply(temp, MARGIN = 2, quantile, probs = c(0.025, 0.975))
    l_ci <- max(0, estimate + quantile[1])
    u_ci <- min(1, estimate + quantile[2])

    ci <- cbind(l_ci, u_ci)
    colnames(ci) <- c("95% CI lower", "95% CI upper")
    ci

    }, vectorize.args = "x")

    # limits for figures by ggplot2
    mean_lim <- c(min(mean_est), max(mean_est))
    acov_lim <- c(min(acov_est), max(acov_est))
    acor_lim <- c(min(acor_est), max(acor_est))

    # computation for confidence intervals
    if (ci) {
    mean_grid <- seq(mean_lim[1], mean_lim[2], length.out = 101)
    acov_grid <- seq(acov_lim[1], acov_lim[2], length.out = 101)
    acor_grid <- seq(acor_lim[1], acor_lim[2], length.out = 101)

    mean_ci <- mean_ci_func(mean_grid)
    acov_ci <- acov_ci_func(acov_grid)
    acor_ci <- acor_ci_func(acor_grid)
    }

    # figures based on HPJ estimation with rearrangement
    if (!ci){
      mean_x <- seq(min(mean_est), max(mean_est), length = 101)
      mean_y <- hpjecdfest1(x = mean_x, X = mean_est, X1 = mean_est1, X2 = mean_est2)
      mean_hpj <- rearrangement(x = data.frame(x = mean_x), y = mean_y)
      mean_plot <- ggplot(data.frame(x = mean_x, y = mean_hpj), aes(x = x, y = y))
      mean_plot <- mean_plot + geom_line() + xlim(min(mean_est), max(mean_est)) + ylim(0, 1)
      mean_plot <- mean_plot + labs(x = "x", y = "")
    }
    if (ci){
      mean_plot <- ggplot(data = data.frame(x = mean_grid), aes(x = mean_grid))
      mean_hpj <- cbind(hpjecdfest1(x = mean_grid, X = mean_est, X1 = mean_est1, X2 = mean_est2), t(mean_ci))
      mean_plot <- mean_plot + geom_line(aes(x = mean_grid, y = rearrangement(x = data.frame(x = mean_grid), y = mean_hpj[, 1])))
      mean_plot <- mean_plot + geom_ribbon(aes(x = mean_grid, ymin = rearrangement(x = data.frame(x = mean_grid), y = mean_hpj[, 2]), ymax = rearrangement(x = data.frame(x = mean_grid), y = mean_hpj[, 3])), alpha = 0.1)
      mean_plot <- mean_plot + labs(x = "x", y = "") + ylim(0, 1)
    }

    if (!ci){
      acov_x <- seq(min(acov_est), max(acov_est), length = 101)
      acov_y <- hpjecdfest1(x = acov_x, X = acov_est, X1 = acov_est1, X2 = acov_est2)
      acov_hpj <- rearrangement(x = data.frame(x = acov_x), y = acov_y)
      acov_plot <- ggplot(data.frame(x = acov_x, y = acov_hpj), aes(x = x, y = y))
      acov_plot <- acov_plot + geom_line() + xlim(min(acov_est), max(acov_est)) + ylim(0, 1)
      acov_plot <- acov_plot + labs(x = "x", y = "")
    }
    if (ci){
      acov_plot <- ggplot(data = data.frame(x = acov_grid), aes(x = acov_grid))
      acov_hpj <- cbind(hpjecdfest1(x = acov_grid, X = acov_est, X1 = acov_est1, X2 = acov_est2), t(acov_ci))
      acov_plot <- acov_plot + geom_line(aes(x = acov_grid, y = rearrangement(x = data.frame(x = acov_grid), y = acov_hpj[, 1])))
      acov_plot <- acov_plot + geom_ribbon(aes(x = acov_grid, ymin = rearrangement(x = data.frame(x = acov_grid), y = acov_hpj[, 2]), ymax = rearrangement(x = data.frame(x = acov_grid), y = acov_hpj[, 3])), alpha = 0.1)
      acov_plot <- acov_plot + labs(x = "x", y = "") + ylim(0, 1)
    }

    if (!ci){
      acor_x <- seq(min(acor_est), max(acor_est), length = 101)
      acor_y <- hpjecdfest1(x = acor_x, X = acor_est, X1 = acor_est1, X2 = acor_est2)
      acor_hpj <- rearrangement(x = data.frame(x = acor_x), y = acor_y)
      acor_plot <- ggplot(data.frame(x = acor_x, y = acor_hpj), aes(x = x, y = y))
      acor_plot <- acor_plot + geom_line() + xlim(min(acor_est), max(acor_est)) + ylim(0, 1)
      acor_plot <- acor_plot + labs(x = "x", y = "")
    }
    if (ci){
      acor_plot <- ggplot(data = data.frame(x = acor_grid), aes(x = acor_grid))
      acor_hpj <- cbind(hpjecdfest1(x = acor_grid, X = acor_est, X1 = acor_est1, X2 = acor_est2), t(acor_ci))
      acor_plot <- acor_plot + geom_line(aes(x = acor_grid, y = rearrangement(x = data.frame(x = acor_grid), y = acor_hpj[, 1])))
      acor_plot <- acor_plot + geom_ribbon(aes(x = acor_grid, ymin = rearrangement(x = data.frame(x = acor_grid), y = acor_hpj[, 2]), ymax = rearrangement(x = data.frame(x = acor_grid), y = acor_hpj[, 3])), alpha = 0.1)
      acor_plot <- acor_plot + labs(x = "x", y = "") + ylim(0, 1)
    }


    # functions without rearrangement
    mean_func <- function(x) {
      hpjecdfest1(x = x, X = mean_est, X1 = mean_est1, X2 = mean_est2)
    }

    acov_func <- function(x) {
      hpjecdfest1(x = x, X = acov_est, X1 = acov_est1, X2 = acov_est2)
    }

    acor_func <- function(x) {
      hpjecdfest1(x = x, X = acor_est, X1 = acor_est1, X2 = acor_est2)
    }

  } else {

    # half panel data for odd T
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

    # function for bootstrap confidence interval
    mean_ci_func <- Vectorize(FUN = function(x){

    # hpj estimation with bootstrap
    bootstrap <- boot(data = cbind(mean_est, mean_est1, mean_est2, mean_est3, mean_est4) - x, statistic = hpj2_boot, R = R)
    estimate <- bootstrap$t0
    temp <- bootstrap$t- bootstrap$t0

    # confidence interval
    quantile <- apply(temp, MARGIN = 2, quantile, probs = c(0.025, 0.975))
    l_ci <- max(0, estimate + quantile[1])
    u_ci <- min(1, estimate + quantile[2])
    ci <- cbind(l_ci, u_ci)
    colnames(ci) <- c("95% CI lower", "95% CI upper")
    ci

    }, vectorize.args = "x")

    acov_ci_func <- Vectorize(FUN = function(x){

    # hpj estimation with bootstrap
    bootstrap <- boot(data = cbind(acov_est, acov_est1, acov_est2, acov_est3, acov_est4) - x, statistic = hpj2_boot, R = R)
    estimate <- bootstrap$t0
    temp <- bootstrap$t- bootstrap$t0

    # confidence interval
    quantile <- apply(temp, MARGIN = 2, quantile, probs = c(0.025, 0.975))
    l_ci <- max(0, estimate + quantile[1])
    u_ci <- min(1, estimate + quantile[2])

    ci <- cbind(l_ci, u_ci)
    colnames(ci) <- c("95% CI lower", "95% CI upper")
    ci

    }, vectorize.args = "x")

    acor_ci_func <- Vectorize(FUN = function(x){

    # hpj estimation with bootstrap
    bootstrap <- boot(data = cbind(acor_est, acor_est1, acor_est2, acor_est3, acor_est4) - x, statistic = hpj2_boot, R = R)
    estimate <- bootstrap$t0
    temp <- bootstrap$t- bootstrap$t0

    # confidence interval
    quantile <- apply(temp, MARGIN = 2, quantile, probs = c(0.025, 0.975))
    l_ci <- max(0, estimate + quantile[1])
    u_ci <- min(1, estimate + quantile[2])

    ci <- cbind(l_ci, u_ci)
    colnames(ci) <- c("95% CI lower", "95% CI upper")
    ci

    }, vectorize.args = "x")

    # limits for figures by ggplot2
    mean_lim <- c(min(mean_est), max(mean_est))
    acov_lim <- c(min(acov_est), max(acov_est))
    acor_lim <- c(min(acor_est), max(acor_est))

    # computation for confidence intervals
    if (ci) {
    mean_grid <- seq(mean_lim[1], mean_lim[2], length.out = 101)
    acov_grid <- seq(acov_lim[1], acov_lim[2], length.out = 101)
    acor_grid <- seq(acor_lim[1], acor_lim[2], length.out = 101)

    mean_ci <- mean_ci_func(mean_grid)
    acov_ci <- acov_ci_func(acov_grid)
    acor_ci <- acor_ci_func(acor_grid)
    }

    # figures based on HPJ estimation with rearrangement
    if (!ci){
      mean_x <- seq(min(mean_est), max(mean_est), length = 101)
      mean_y <- hpjecdfest2(x = mean_x, X = mean_est, X1 = mean_est1, X2 = mean_est2, X3 = mean_est3, X4 = mean_est4)
      mean_hpj <- rearrangement(x = data.frame(x = mean_x), y = mean_y)
      mean_plot <- ggplot(data.frame(x = mean_x, y = mean_hpj), aes(x = x, y = y))
      mean_plot <- mean_plot + geom_line() + xlim(min(mean_est), max(mean_est)) + ylim(0, 1)
      mean_plot <- mean_plot + labs(x = "x", y = "")
    }
    if (ci){
      mean_plot <- ggplot(data = data.frame(x = mean_grid), aes(x = mean_grid))
      mean_hpj <- cbind(hpjecdfest2(x = mean_grid, X = mean_est, X1 = mean_est1, X2 = mean_est2, X3 = mean_est3, X4 = mean_est4), t(mean_ci))
      mean_plot <- mean_plot + geom_line(aes(x = mean_grid, y = rearrangement(x = data.frame(x = mean_grid), y = mean_hpj[, 1])))
      mean_plot <- mean_plot + geom_ribbon(aes(x = mean_grid, ymin = rearrangement(x = data.frame(x = mean_grid), y = mean_hpj[, 2]), ymax = rearrangement(x = data.frame(x = mean_grid), y = mean_hpj[, 3])), alpha = 0.1)
      mean_plot <- mean_plot + labs(x = "x", y = "") + ylim(0, 1)
    }

    if (!ci){
      acov_x <- seq(min(acov_est), max(acov_est), length = 101)
      acov_y <- hpjecdfest2(x = acov_x, X = acov_est, X1 = acov_est1, X2 = acov_est2, X3 = acov_est3, X4 = acov_est4)
      acov_hpj <- rearrangement(x = data.frame(x = acov_x), y = acov_y)
      acov_plot <- ggplot(data.frame(x = acov_x, y = acov_hpj), aes(x = x, y = y))
      acov_plot <- acov_plot + geom_line() + xlim(min(acov_est), max(acov_est)) + ylim(0, 1)
      acov_plot <- acov_plot + labs(x = "x", y = "")
    }
    if (ci){
      acov_plot <- ggplot(data = data.frame(x = acov_grid), aes(x = acov_grid))
      acov_hpj <- cbind(hpjecdfest2(x = acov_grid, X = acov_est, X1 = acov_est1, X2 = acov_est2, X3 = acov_est3, X4 = acov_est4), t(acov_ci))
      acov_plot <- acov_plot + geom_line(aes(x = acov_grid, y = rearrangement(x = data.frame(x = acov_grid), y = acov_hpj[, 1])))
      acov_plot <- acov_plot + geom_ribbon(aes(x = acov_grid, ymin = rearrangement(x = data.frame(x = acov_grid), y = acov_hpj[, 2]), ymax = rearrangement(x = data.frame(x = acov_grid), y = acov_hpj[, 3])), alpha = 0.1)
      acov_plot <- acov_plot + labs(x = "x", y = "") + ylim(0, 1)
    }

    if (!ci){
      acor_x <- seq(min(acor_est), max(acor_est), length = 101)
      acor_y <- hpjecdfest2(x = acor_x, X = acor_est, X1 = acor_est1, X2 = acor_est2, X3 = acor_est3, X4 = acor_est4)
      acor_hpj <- rearrangement(x = data.frame(x = acor_x), y = acor_y)
      acor_plot <- ggplot(data.frame(x = acor_x, y = acor_hpj), aes(x = x, y = y))
      acor_plot <- acor_plot + geom_line() + xlim(min(acor_est), max(acor_est)) + ylim(0, 1)
      acor_plot <- acor_plot + labs(x = "x", y = "")
    }
    if (ci){
      acor_plot <- ggplot(data = data.frame(x = acor_grid), aes(x = acor_grid))
      acor_hpj <- cbind(hpjecdfest2(x = acor_grid, X = acor_est, X1 = acor_est1, X2 = acor_est2, X3 = acor_est3, X4 = acor_est4), t(acor_ci))
      acor_plot <- acor_plot + geom_line(aes(x = acor_grid, y = rearrangement(x = data.frame(x = acor_grid), y = acor_hpj[, 1])))
      acor_plot <- acor_plot + geom_ribbon(aes(x = acor_grid, ymin = rearrangement(x = data.frame(x = acor_grid), y = acor_hpj[, 2]), ymax = rearrangement(x = data.frame(x = acor_grid), y = acor_hpj[, 3])), alpha = 0.1)
      acor_plot <- acor_plot + labs(x = "x", y = "") + ylim(0, 1)
    }


    # functions without rearrangement
    mean_func <- function(x) {
      hpjecdfest2(x = x, X = mean_est, X1 = mean_est1, X2 = mean_est2, X3 = mean_est3, X4 = mean_est4)
    }

    acov_func <- function(x) {
      hpjecdfest2(x = x, X = acov_est, X1 = acov_est1, X2 = acov_est2, X3 = acov_est3, X4 = acov_est4)
    }

    acor_func <- function(x) {
      hpjecdfest2(x = x, X = acor_est, X1 = acor_est1, X2 = acor_est2, X3 = acor_est3, X4 = acor_est4)
    }

  }

  # results
  quantity <- cbind(mean_est, acov_est, acor_est)
  colnames(quantity) <- c("mean", "autocovariance", "autocorrelation")
  result <- list(mean = mean_plot, acov = acov_plot, acor = acor_plot,
                 mean_func = mean_func, acov_func = acov_func, acor_func = acor_func,
                 mean_ci_func = mean_ci_func, acov_ci_func = acov_ci_func, acor_ci_func = acor_ci_func,
                 quantity = quantity, acov_order = acov_order,
                 acor_order = acor_order, N = N, S = S, R = R)

  return(result)

}


#' computing HPJ empirical CDF estimate for even T
#'
#' @param x point at which the CDF is estimated
#' @param X vector of original cross-sectional data
#' @param X1 vector of half-panel cross-sectional data based on time series 1 ~ T/2
#' @param X2 vector of half-panel cross-sectional data based on time series (T/2 + 1) ~ T
#'
hpjecdfest1 <- Vectorize(FUN = function(x, X, X1, X2) {

  # estimates
  est <- mean(ifelse(X <= x, 1, 0))
  est1 <- mean(ifelse(X1 <= x, 1, 0))
  est2 <- mean(ifelse(X2 <= x, 1, 0))

  # HPJ estimate
  hpjest <- 2 * est - (est1 + est2) / 2

  # correction to ensure valid estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)
  hpjest <- ifelse(hpjest <= 1, hpjest, 1)

  return(hpjest)

}, vectorize.args = "x")


#' computing HPJ empirical CDF estimate for odd T
#'
#' @param x point at which the CDF is estimated
#' @param X vector of original cross-sectional data
#' @param X1 vector of half-panel cross-sectional data based on time series 1 ~ floor(T/2)
#' @param X2 vector of half-panel cross-sectional data based on time series (floor(T/2) + 1) ~ T
#' @param X3 vector of half-panel cross-sectional data based on time series 1 ~ ceiling(T/2)
#' @param X4 vector of half-panel cross-sectional data based on time series (ceiling(T/2) + 1) ~ T
#'
hpjecdfest2 <- Vectorize(FUN = function(x, X, X1, X2, X3, X4) {

  # estimates
  est <- mean(ifelse(X <= x, 1, 0))
  est1 <- mean(ifelse(X1 <= x, 1, 0))
  est2 <- mean(ifelse(X2 <= x, 1, 0))
  est3 <- mean(ifelse(X3 <= x, 1, 0))
  est4 <- mean(ifelse(X4 <= x, 1, 0))

  # HPJ estimate
  hpjest <- 2 * est - (est1 + est2 + est3 + est4) / 4

  # correction to ensure valid estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)
  hpjest <- ifelse(hpjest <= 1, hpjest, 1)

  return(hpjest)

}, vectorize.args = "x")

#' computing bootstrap HPJ empirical CDF estimate for odd T
#'
#' @param quantity N * 3 matrix of estimates
#' @param indices indices for bootstrap replications
#'
hpj1_boot <- function(quantity, indices){

  # estimates
  est1 <- mean(ifelse(quantity[indices, 1] <= 0, 1, 0))
  est2 <- mean(ifelse(quantity[indices, 2] <= 0, 1, 0))
  est3 <- mean(ifelse(quantity[indices, 3] <= 0, 1, 0))

  # HPJ estimate
  hpjest <- 2 * est1 - (est2 + est3) / 2

  # correction to ensure valid estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)
  hpjest <- ifelse(hpjest <= 1, hpjest, 1)

  return(hpjest)

}

#' computing bootstrap HPJ empirical CDF estimate for even T
#'
#' @param quantity N * 5 matrix of estimates
#' @param indices indices for bootstrap replications
#'
hpj2_boot <- function(quantity, indices){

  est1 <- mean(ifelse(quantity[indices, 1] <= 0, 1, 0))
  est2 <- mean(ifelse(quantity[indices, 2] <= 0, 1, 0))
  est3 <- mean(ifelse(quantity[indices, 3] <= 0, 1, 0))
  est4 <- mean(ifelse(quantity[indices, 4] <= 0, 1, 0))
  est5 <- mean(ifelse(quantity[indices, 5] <= 0, 1, 0))

  # HPJ estimate
  hpjest <- 2 * est1 - (est2 + est3 + est4 + est5) / 4

  # correction to ensure valid estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)
  hpjest <- ifelse(hpjest <= 1, hpjest, 1)

  return(hpjest)

}