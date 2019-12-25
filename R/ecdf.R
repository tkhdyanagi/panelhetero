#' Empirical CDF estimation for heterogeneity in panel data without bias-correction
#'
#' \code{neecdf} implements the naive estimation of the empirical CDF
#' for the heterogeneous mean, autocovariance, and autocorrelation.
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
#' @importFrom stats na.omit
#'
#' @return list that contains the following elements.
#' \item{mean}{graph of the naive empirical CDF estimation for the mean}
#' \item{acov}{graph of the naive empirical CDF estimation for the autocovariance}
#' \item{acor}{graph of the naive empirical CDF estimation for the autocorrelation}
#' \item{mean_func}{function that returns naive empirical CDF estimates for the mean}
#' \item{acov_func}{function that returns naive empirical CDF estimates for the autocovariance}
#' \item{acor_func}{function that returns naive empirical CDF estimates for the autocorrelation}
#' \item{mean_ci_func}{function that returns 95 percent bootstrap confidence interval for naive empirical CDF estimates for the mean}
#' \item{acov_ci_func}{function that returns 95 percent bootstrap confidence interval for naive empirical CDF estimates for the autocovariance}
#' \item{acor_ci_func}{function that returns 95 percent bootstrap confidence interval for naive empirical CDF estimates for the autocorrelation}
#' \item{quantity}{matrix that contains the estimated quantities}
#' \item{acov_order}{the same as the argument}
#' \item{acor_order}{the same as the argument}
#' \item{N}{the number of cross-sectional units}
#' \item{S}{the length of time series}
#' \item{R}{the number of bootstrap replications}
#'
#' @export
#'
neecdf <- function(data, acov_order = 0, acor_order = 1, R = 1000, ci = TRUE) {

  # initialization
  x <- NULL

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

  # function for bootstrap confidence interval
  mean_ci_func <- Vectorize(FUN = function(x){

    # naive estimation with bootstrap
    bootstrap <- boot(data = mean_est - x, statistic = ne_boot, R = R)
    estimate <- bootstrap$t0
    temp <- bootstrap$t - bootstrap$t0

    # confidence interval
    quantile <- apply(temp, MARGIN = 2, quantile, probs = c(0.025, 0.975))
    l_ci <- max(0, estimate + quantile[1])
    u_ci <- min(1, estimate + quantile[2])
    ci <- cbind(l_ci, u_ci)
    colnames(ci) <- c("95% CI lower", "95% CI upper")
    ci

  }, vectorize.args = "x")

  acov_ci_func <- Vectorize(FUN = function(x){

    # naive estimation with bootstrap
    bootstrap <- boot(data = acov_est - x, statistic = ne_boot, R = R)
    estimate <- bootstrap$t0
    temp <- bootstrap$t - bootstrap$t0

    # confidence interval
    quantile <- apply(temp, MARGIN = 2, quantile, probs = c(0.025, 0.975))
    l_ci <- max(0, estimate + quantile[1])
    u_ci <- min(1, estimate + quantile[2])

    ci <- cbind(l_ci, u_ci)
    colnames(ci) <- c("95% CI lower", "95% CI upper")
    ci

  }, vectorize.args = "x")

  acor_ci_func <- Vectorize(FUN = function(x){

    # naive estimation with bootstrap
    bootstrap <- boot(data = acor_est - x, statistic = ne_boot, R = R)
    estimate <- bootstrap$t0
    temp <- bootstrap$t - bootstrap$t0

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

  # figures by ggplot2
  if (!ci){
    mean_plot <- ggplot(data = data.frame(x = mean_lim), aes(x = x))
    mean_plot <- mean_plot + stat_function(fun = ecdfest, args = list(X = mean_est), geom = "step")
  }

  if (ci){
    mean_plot <- ggplot(data = data.frame(x = mean_grid), aes(x = mean_grid))
    mean_plot <- mean_plot + geom_step(aes(x = mean_grid, y = ecdfest(mean_grid, mean_est)))
    mean_plot <- mean_plot + geom_ribbon(aes(x = mean_grid, ymin = rearrangement(x = data.frame(x = mean_grid), y = mean_ci[1, ]), ymax = rearrangement(x = data.frame(x = mean_grid), y = mean_ci[2, ])), alpha = 0.1)
  }
  mean_plot <- mean_plot + labs(x = "x", y = "") + ylim(0, 1)

  if (!ci){
    acov_plot <- ggplot(data = data.frame(x = acov_lim), aes(x = x))
    acov_plot <- acov_plot + stat_function(fun = ecdfest, args = list(X = acov_est), geom = "step")
  }

  if (ci){
    acov_plot <- ggplot(data = data.frame(x = acov_grid), aes(x = acov_grid))
    acov_plot <- acov_plot + geom_step(aes(x = acov_grid, y = ecdfest(acov_grid, acov_est)))
    acov_plot <- acov_plot + geom_ribbon(aes(x = acov_grid, ymin = rearrangement(x = data.frame(x = acov_grid), y = acov_ci[1, ]), ymax = rearrangement(x = data.frame(x = acov_grid), y = acov_ci[2, ])), alpha = 0.1)
  }
  acov_plot <- acov_plot + labs(x = "x", y = "") + ylim(0, 1)

  if (!ci){
    acor_plot <- ggplot(data = data.frame(x = acor_lim), aes(x = x))
    acor_plot <- acor_plot + stat_function(fun = ecdfest, args = list(X = acor_est), geom = "step")
  }
  if (ci){
    acor_plot <- ggplot(data = data.frame(x = acor_grid), aes(x = acor_grid))
    acor_plot <- acor_plot + geom_step(aes(x = acor_grid, y = ecdfest(acor_grid, acor_est)))
    acor_plot <- acor_plot + geom_ribbon(aes(x = acor_grid, ymin = rearrangement(x = data.frame(x = acor_grid), y = acor_ci[1, ]), ymax = rearrangement(x = data.frame(x = acor_grid), y = acor_ci[2, ])), alpha = 0.1)
  }
  acor_plot <- acor_plot + labs(x = "x", y = "") + ylim(0, 1)

  # functions
  mean_func <- function(x) {
    ecdfest(x = x, X = mean_est)
  }

  acov_func <- function(x) {
    ecdfest(x = x, X = acov_est)
  }

  acor_func <- function(x) {
    ecdfest(x = x, X = acor_est)
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


#' computing empirical CDF estimate
#'
#' @param x point at which the CDF is estimated
#' @param X vector of cross-sectional data
#'
ecdfest <- Vectorize(FUN = function(x, X) {

  est <- mean(ifelse(X <= x, 1, 0))

  return(est)

}, vectorize.args = "x")

#' computing bootstrap empirical CDF estimate
#'
#' @param quantity N * 1 matrix of estimate
#' @param indices indices for bootstrap replications
#'
ne_boot <- function(quantity, indices){

  est <- mean(ifelse(quantity[indices] <= 0, 1, 0))

  return(est)

}