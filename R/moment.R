#' The naive estimation of the moments
#'
#' The `nemoment()` function enables to implement the naive estimation of
#' the moments of the heterogeneous mean, the heterogeneous autocovariance,
#' and the heterogeneous autocorrelation.
#' The method is developed by Okui and Yanagi (2019).
#' For more details, see the package vignette with `vignette("panelhetero")`.
#'
#' @param data A matrix of panel data.
#' Each row corresponds to individual time series.
#' @param acov_order A non-negative integer of the order of autocovariance.
#' Default is 0.
#' @param acor_order A positive integer of the order of autocorrelation
#' Default is 1.
#' @param R A positive integer of the number of bootstrap repetitions.
#' Default is 1000.
#'
#' @returns A list that contains the following elements.
#' \item{estimate}{A vector of the parameter estimates}
#' \item{se}{A vector of the standard errors}
#' \item{ci}{A matrix of the 95 percent confidence intervals}
#' \item{quantity}{A matrix of the estimated heterogeneous quantities}
#' \item{acov_order}{The order of autocovariance}
#' \item{acor_order}{The order of autocovariance}
#' \item{N}{The number of cross-sectional units}
#' \item{S}{The length of time series}
#' \item{R}{The number of bootstrap repetitions}
#'
#' @examples
#' data <- panelhetero::simulation(N = 300, S = 50)
#' panelhetero::nemoment(data = data)
#'
#' @references Okui, R. and Yanagi, T., 2019.
#' Panel data analysis with heterogeneous dynamics.
#' Journal of Econometrics, 212(2), pp.451-475.
#'
#' @export
#'
nemoment <- function(data,
                     acov_order = 0,
                     acor_order = 1,
                     R = 1000) {

  # Error handling

  error1(data = data, acov_order = acov_order, acor_order = acor_order, R = R)

  # Omit NA
  data <- stats::na.omit(data)

  # Sample size
  N <- nrow(data)
  S <- ncol(data)

  # Estimated means, autocovariances, autocovariances
  mean_est <- rowMeans(data)
  acov_est <- apply(data, MARGIN = 1, acov, acov_order = acov_order)
  acor_est <- apply(data, MARGIN = 1, acor, acor_order = acor_order)
  equantity <- cbind(mean_est, acov_est, acor_est)

  # Naive estimation with bootstrap
  bootstrap <- boot::boot(data = equantity, statistic = momentest, R = R)

  # Estimates
  estimate <- bootstrap$t0

  # Standard errors
  se <- apply(bootstrap$t, MARGIN = 2, stats::sd)

  # Confidence intervals
  temp <- t(t(bootstrap$t) - bootstrap$t0)
  quantiles <- apply(temp,
                     MARGIN = 2,
                     stats::quantile,
                     probs = c(0.025, 0.975))
  ci <- cbind(estimate + quantiles[1, ], estimate + quantiles[2, ])

  # Names
  names(estimate) <- names(se) <- rownames(ci) <-
    c("E(mean)", "E(acov)", "E(acor)",
      "var(mean)", "var(acov)", "var(acor)",
      "cor(mean, acov)", "cor(mean, acor)", "cor(acov, acor)")

  colnames(ci) <- c("95% CI lower", "95% CI upper")

  # result
  quantity <- cbind(mean_est,
                    acov_est,
                    acor_est)

  colnames(quantity) <- c("mean",
                          "autocovariance",
                          "autocorrelation")

  return(list(estimate = estimate,
              se = se,
              ci = ci,
              quantity = quantity,
              acov_order = acov_order,
              acor_order = acor_order,
              N = N,
              S = S,
              R = R)
  )

}

#' Compute the estimates of the moments
#'
#' @param quantity An N * 3 matrix of the estimated quantities
#' @param indices A vector of indices for bootstrap repetitions
#'
momentest <- function(quantity, indices) {

  # Estimated quantities
  mean_est <- quantity[indices, 1]
  acov_est <- quantity[indices, 2]
  acor_est <- quantity[indices, 3]

  # Means
  mean_mean <- mean(mean_est)
  acov_mean <- mean(acov_est)
  acor_mean <- mean(acor_est)

  # Variances
  mean_var <- stats::var(mean_est)
  acov_var <- stats::var(acov_est)
  acor_var <- stats::var(acor_est)

  # Correlations
  mean_acov_cor <- stats::cor(mean_est, acov_est)
  mean_acor_cor <- stats::cor(mean_est, acor_est)
  acov_acor_cor <- stats::cor(acov_est, acor_est)

  # Return
  return(c(mean_mean,
           acov_mean,
           acor_mean,
           mean_var,
           acov_var,
           acor_var,
           mean_acov_cor,
           mean_acor_cor,
           acov_acor_cor)
  )

}