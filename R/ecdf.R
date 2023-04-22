#' The naive empirical CDF estimation without bias correction
#'
#' The `neecdf()` function enables to implement the naive estimation of
#' the cumulative distribution function (CDF) of the heterogeneous mean,
#' the heterogeneous autocovariance, and the heterogeneous autocorrelation.
#' The method is developed by Okui and Yanagi (2019).
#' For more details, see the package vignette with `vignette("panelhetero")`.
#'
#' @param data A matrix of panel data.
#' Each row corresponds to individual time series.
#' @param acov_order A non-negative integer of the order of autocovariance.
#' Default is 0.
#' @param acor_order A positive integer of the order of autocorrelation.
#' Default is 1.
#' @param ci A logical whether to estimate the confidence interval.
#' Default is TRUE.
#' @param R A positive integer of the number of bootstrap repetitions.
#' Default is 1000.
#'
#' @returns A list that contains the following elements.
#' \item{mean}{A plot of the corresponding CDF}
#' \item{acov}{A plot of the corresponding CDF}
#' \item{acor}{A plot of the corresponding CDF}
#' \item{mean_func}{A function that returns the corresponding CDF}
#' \item{acov_func}{A function that returns the corresponding CDF}
#' \item{acor_func}{A function that returns the corresponding CDF}
#' \item{mean_ci_func}{A function that returns the 95 percent confidence
#' interval for the corresponding CDF}
#' \item{acov_ci_func}{A function that returns the 95 percent confidence
#' interval for the corresponding CDF}
#' \item{acor_ci_func}{A function that returns the 95 percent confidence
#' interval for the corresponding CDF}
#' \item{quantity}{A matrix of the estimated heterogeneous quantities}
#' \item{acov_order}{The order of autocovariance}
#' \item{acor_order}{The order of autocorrelation}
#' \item{N}{The number of cross-sectional units}
#' \item{S}{The length of time series}
#' \item{R}{The number of bootstrap repetitions}
#'
#' @examples
#' data <- panelhetero::simulation(N = 300, S = 50)
#' panelhetero::neecdf(data = data, R = 50)
#'
#' @references Okui, R. and Yanagi, T., 2019.
#' Panel data analysis with heterogeneous dynamics.
#' Journal of Econometrics, 212(2), pp.451-475.
#'
#' @export
#'
neecdf <- function(data,
                   acov_order = 0,
                   acor_order = 1,
                   R = 1000,
                   ci = TRUE) {

  # Error handling -------------------------------------------------------------

  error1(data = data,
         acov_order = acov_order,
         acor_order = acor_order,
         R = R)

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

  # Function for the confidence interval
  mean_ci_func <- Vectorize(FUN = function(x) {

    # naive estimation with bootstrap
    bootstrap <- boot::boot(data = mean_est - x, statistic = ne_boot, R = R)
    estimate <- bootstrap$t0
    temp <- bootstrap$t - bootstrap$t0

    # confidence interval
    quantile <- apply(temp,
                      MARGIN = 2,
                      stats::quantile,
                      probs = c(0.025, 0.975))

    l_ci <- max(0, estimate + quantile[1])
    u_ci <- min(1, estimate + quantile[2])

    ci <- cbind(l_ci, u_ci)
    colnames(ci) <- c("95% CI lower", "95% CI upper")

    return(ci)

  }, vectorize.args = "x")

  acov_ci_func <- Vectorize(FUN = function(x) {

    # naive estimation with bootstrap
    bootstrap <- boot::boot(data = acov_est - x, statistic = ne_boot, R = R)
    estimate <- bootstrap$t0
    temp <- bootstrap$t - bootstrap$t0

    # confidence interval
    quantile <- apply(temp,
                      MARGIN = 2,
                      stats::quantile,
                      probs = c(0.025, 0.975))

    l_ci <- max(0, estimate + quantile[1])
    u_ci <- min(1, estimate + quantile[2])

    ci <- cbind(l_ci, u_ci)
    colnames(ci) <- c("95% CI lower", "95% CI upper")

    return(ci)

  }, vectorize.args = "x")

  acor_ci_func <- Vectorize(FUN = function(x) {

    # naive estimation with bootstrap
    bootstrap <- boot::boot(data = acor_est - x, statistic = ne_boot, R = R)
    estimate <- bootstrap$t0
    temp <- bootstrap$t - bootstrap$t0

    # confidence interval
    quantile <- apply(temp,
                      MARGIN = 2,
                      stats::quantile,
                      probs = c(0.025, 0.975))

    l_ci <- max(0, estimate + quantile[1])
    u_ci <- min(1, estimate + quantile[2])

    ci <- cbind(l_ci, u_ci)
    colnames(ci) <- c("95% CI lower", "95% CI upper")

    return(ci)

  }, vectorize.args = "x")

  # Limits used for ggplot2
  mean_lim <- c(min(mean_est),
                max(mean_est))

  acov_lim <- c(min(acov_est),
                max(acov_est))

  acor_lim <- c(min(acor_est),
                max(acor_est))

  # Compute confidence intervals
  if (ci) {

    mean_grid <- seq(mean_lim[1], mean_lim[2], length.out = 101)
    acov_grid <- seq(acov_lim[1], acov_lim[2], length.out = 101)
    acor_grid <- seq(acor_lim[1], acor_lim[2], length.out = 101)

    mean_ci <- mean_ci_func(mean_grid)
    acov_ci <- acov_ci_func(acov_grid)
    acor_ci <- acor_ci_func(acor_grid)
  }

  # Make figures using ggplot2
  if (!ci) {

    # Mean
    mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = ecdfest,
                             args = list(X = mean_est),
                             geom = "step") +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous mean") +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_bw()

    # Autocovariance
    acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = ecdfest,
                             args = list(X = acov_est),
                             geom = "step") +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocovariance") +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_bw()

    # Autocorrelation
    acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_lim),
                                 ggplot2::aes(x = x)) +
      ggplot2::stat_function(fun = ecdfest,
                             args = list(X = acor_est),
                             geom = "step") +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocorrelation") +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_bw()

  }

  if (ci) {

    # Mean
    mean_plot <- ggplot2::ggplot(
      data = data.frame(x = mean_grid),
      ggplot2::aes(x = mean_grid)) +
      ggplot2::geom_step(
        ggplot2::aes(x = mean_grid,
                     y = ecdfest(mean_grid, mean_est))) +
      ggplot2::geom_ribbon(
        ggplot2::aes(x = mean_grid,
                     ymin = Rearrangement::rearrangement(
                       x = data.frame(x = mean_grid),
                       y = mean_ci[1, ]
                     ),
                     ymax = Rearrangement::rearrangement(
                       x = data.frame(x = mean_grid),
                       y = mean_ci[2, ]
                     )
        ),
        alpha = 0.1) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous mean") +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_bw()

    # Autocovariance
    acov_plot <- ggplot2::ggplot(
      data = data.frame(x = acov_grid),
      ggplot2::aes(x = acov_grid)) +
      ggplot2::geom_step(
        ggplot2::aes(x = acov_grid,
                     y = ecdfest(acov_grid, acov_est))) +
      ggplot2::geom_ribbon(
        ggplot2::aes(x = acov_grid,
                     ymin = Rearrangement::rearrangement(
                       x = data.frame(x = acov_grid),
                       y = acov_ci[1, ]
                     ),
                     ymax = Rearrangement::rearrangement(
                       x = data.frame(x = acov_grid),
                       y = acov_ci[2, ]
                     )
        ),
        alpha = 0.1) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocovariance") +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_bw()

    # Autocorrelation
    acor_plot <- ggplot2::ggplot(
      data = data.frame(x = acor_grid),
      ggplot2::aes(x = acor_grid)) +
      ggplot2::geom_step(
        ggplot2::aes(x = acor_grid,
                     y = ecdfest(acor_grid, acor_est))) +
      ggplot2::geom_ribbon(
        ggplot2::aes(x = acor_grid,
                     ymin = Rearrangement::rearrangement(
                       x = data.frame(x = acor_grid),
                       y = acor_ci[1, ]
                     ),
                     ymax = Rearrangement::rearrangement(
                       x = data.frame(x = acor_grid),
                       y = acor_ci[2, ]
                     )
        ),
        alpha = 0.1) +
      ggplot2::labs(x = "x", y = "") +
      ggplot2::ggtitle("The heterogeneous autocorrelation") +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_bw()

  }

  # Functions
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
  quantity <- cbind(mean_est,
                    acov_est,
                    acor_est)

  colnames(quantity) <- c("mean",
                          "autocovariance",
                          "autocorrelation")

  return(list(mean = mean_plot,
              acov = acov_plot,
              acor = acor_plot,
              mean_func = mean_func,
              acov_func = acov_func,
              acor_func = acor_func,
              mean_ci_func = mean_ci_func,
              acov_ci_func = acov_ci_func,
              acor_ci_func = acor_ci_func,
              quantity = quantity,
              acov_order = acov_order,
              acor_order = acor_order,
              N = N,
              S = S,
              R = R)
  )

}

#' Compute empirical CDF estimate
#'
#' @param x An evaluation point
#' @param X A vector of cross-sectional data
#'
#' @noRd
#'
ecdfest <- Vectorize(FUN = function(x, X) {

  return(mean(X <= x))

}, vectorize.args = "x")

#' Compute bootstrap empirical CDF estimate
#'
#' @param quantity An N * 1 matrix of estimates
#' @param indices A vector of indices for bootstrap repetitions
#'
#' @noRd
#'
ne_boot <- function(quantity, indices) {

  return(mean(quantity[indices] <= 0))

}