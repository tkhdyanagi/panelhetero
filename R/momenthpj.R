#' The HPJ bias-corrected estimation of the moments
#'
#' The `hpjmoment()` function enables to implement the HPJ bias-corrected
#' estimation of the moments of the heterogeneous mean,
#' the heterogeneous autocovariance, and the heterogeneous autocorrelation.
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
#' panelhetero::hpjmoment(data = data)
#'
#' @references Okui, R. and Yanagi, T., 2019.
#' Panel data analysis with heterogeneous dynamics.
#' Journal of Econometrics, 212(2), pp.451-475.
#'
#' @export
#'
hpjmoment <- function(data,
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

  # HPJ bias-corrected estimation
  if (S %% 2 == 0) {

    # Half-panel data for even T
    data1 <- data[, 1:(S / 2)]
    data2 <- data[, (S / 2 + 1):S]

    # Estimated quangtities
    mean_est1 <- rowMeans(data1)
    mean_est2 <- rowMeans(data2)

    acov_est1 <- apply(data1, MARGIN = 1, acov, acov_order = acov_order)
    acov_est2 <- apply(data2, MARGIN = 1, acov, acov_order = acov_order)

    acor_est1 <- apply(data1, MARGIN = 1, acor, acor_order = acor_order)
    acor_est2 <- apply(data2, MARGIN = 1, acor, acor_order = acor_order)

    equantity <- cbind(mean_est,
                       mean_est1,
                       mean_est2,
                       acov_est,
                       acov_est1,
                       acov_est2,
                       acor_est,
                       acor_est1,
                       acor_est2)

    # Estimation with bootstrap
    bootstrap <- boot::boot(data = equantity,
                            statistic = hpjmomentest1,
                            R = R)

  } else {

    # Half-panel data for odd T
    data1 <- data[, 1:floor(S / 2)]
    data2 <- data[, (floor(S / 2) + 1):S]
    data3 <- data[, 1:ceiling(S / 2)]
    data4 <- data[, (ceiling(S / 2) + 1):S]

    # Estimated quantities
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

    equantity <- cbind(mean_est,
                       mean_est1,
                       mean_est2,
                       mean_est3,
                       mean_est4,
                       acov_est,
                       acov_est1,
                       acov_est2,
                       acov_est3,
                       acov_est4,
                       acor_est,
                       acor_est1,
                       acor_est2,
                       acor_est3,
                       acor_est4)

    # Estimation with bootstrap
    bootstrap <- boot::boot(data = equantity,
                            statistic = hpjmomentest2,
                            R = R)

  }

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
  names(estimate) <- names(se) <- rownames(ci) <- c("E(mean)",
                                                    "E(acov)",
                                                    "E(acor)",
                                                    "var(mean)",
                                                    "var(acov)",
                                                    "var(acor)",
                                                    "cor(mean, acov)",
                                                    "cor(mean, acor)",
                                                    "cor(acov, acor)")

  colnames(ci) <- c("95% CI lower", "95% CI upper")

  # Results
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

#' Compute HPJ bias-corrected estimates of the moments for even T
#'
#' @param quantity An N * 3 matrix of the estimated quantities
#' @param indices A vector of indices for bootstrap repetitions
#'
#' @noRd
#'
hpjmomentest1 <- function(quantity, indices) {

  # Estimated quantities
  mean_est  <- quantity[indices, 1]
  mean_est1 <- quantity[indices, 2]
  mean_est2 <- quantity[indices, 3]

  acov_est  <- quantity[indices, 4]
  acov_est1 <- quantity[indices, 5]
  acov_est2 <- quantity[indices, 6]

  acor_est  <- quantity[indices, 7]
  acor_est1 <- quantity[indices, 8]
  acor_est2 <- quantity[indices, 9]

  # Means
  mean_mean  <- mean(mean_est)
  mean_mean1 <- mean(mean_est1)
  mean_mean2 <- mean(mean_est2)

  acov_mean  <- mean(acov_est)
  acov_mean1 <- mean(acov_est1)
  acov_mean2 <- mean(acov_est2)

  acor_mean  <- mean(acor_est)
  acor_mean1 <- mean(acor_est1)
  acor_mean2 <- mean(acor_est2)

  # Variances
  mean_var  <- stats::var(mean_est)
  mean_var1 <- stats::var(mean_est1)
  mean_var2 <- stats::var(mean_est2)

  acov_var  <- stats::var(acov_est)
  acov_var1 <- stats::var(acov_est1)
  acov_var2 <- stats::var(acov_est2)

  acor_var  <- stats::var(acor_est)
  acor_var1 <- stats::var(acor_est1)
  acor_var2 <- stats::var(acor_est2)

  # Correlations
  mean_acov_cor  <- stats::cor(mean_est, acov_est)
  mean_acov_cor1 <- stats::cor(mean_est1, acov_est1)
  mean_acov_cor2 <- stats::cor(mean_est2, acov_est2)

  mean_acor_cor  <- stats::cor(mean_est, acor_est)
  mean_acor_cor1 <- stats::cor(mean_est1, acor_est1)
  mean_acor_cor2 <- stats::cor(mean_est2, acor_est2)

  acov_acor_cor  <- stats::cor(acov_est, acor_est)
  acov_acor_cor1 <- stats::cor(acov_est1, acor_est1)
  acov_acor_cor2 <- stats::cor(acov_est2, acor_est2)

  # Estimates
  estimate  <- c(mean_mean,
                 acov_mean,
                 acor_mean,
                 mean_var,
                 acov_var,
                 acor_var,
                 mean_acov_cor,
                 mean_acor_cor,
                 acov_acor_cor)

  estimate1 <- c(mean_mean1,
                 acov_mean1,
                 acor_mean1,
                 mean_var1,
                 acov_var1,
                 acor_var1,
                 mean_acov_cor1,
                 mean_acor_cor1,
                 acov_acor_cor1)

  estimate2 <- c(mean_mean2,
                 acov_mean2,
                 acor_mean2,
                 mean_var2,
                 acov_var2,
                 acor_var2,
                 mean_acov_cor2,
                 mean_acor_cor2,
                 acov_acor_cor2)

  # HPJ bias-corrected estimate
  hpjestimate <- 2 * estimate - (estimate1 + estimate2) / 2

  return(hpjestimate)

}

#' Compute HPJ bias-corrected estimates of the moments for odd T
#'
#' @param quantity An N * 3 matrix of the estimated quantities
#' @param indices A vector of indices for bootstrap repetitions
#'
#' @noRd
#'
hpjmomentest2 <- function(quantity, indices) {

  # Estimated quantities
  mean_est  <- quantity[indices, 1]
  mean_est1 <- quantity[indices, 2]
  mean_est2 <- quantity[indices, 3]
  mean_est3 <- quantity[indices, 4]
  mean_est4 <- quantity[indices, 5]

  acov_est  <- quantity[indices, 6]
  acov_est1 <- quantity[indices, 7]
  acov_est2 <- quantity[indices, 8]
  acov_est3 <- quantity[indices, 9]
  acov_est4 <- quantity[indices, 10]

  acor_est  <- quantity[indices, 11]
  acor_est1 <- quantity[indices, 12]
  acor_est2 <- quantity[indices, 13]
  acor_est3 <- quantity[indices, 14]
  acor_est4 <- quantity[indices, 15]

  # Means
  mean_mean  <- mean(mean_est)
  mean_mean1 <- mean(mean_est1)
  mean_mean2 <- mean(mean_est2)
  mean_mean3 <- mean(mean_est3)
  mean_mean4 <- mean(mean_est4)

  acov_mean  <- mean(acov_est)
  acov_mean1 <- mean(acov_est1)
  acov_mean2 <- mean(acov_est2)
  acov_mean3 <- mean(acov_est3)
  acov_mean4 <- mean(acov_est4)

  acor_mean  <- mean(acor_est)
  acor_mean1 <- mean(acor_est1)
  acor_mean2 <- mean(acor_est2)
  acor_mean3 <- mean(acor_est3)
  acor_mean4 <- mean(acor_est4)

  # Variances
  mean_var  <- stats::var(mean_est)
  mean_var1 <- stats::var(mean_est1)
  mean_var2 <- stats::var(mean_est2)
  mean_var3 <- stats::var(mean_est3)
  mean_var4 <- stats::var(mean_est4)

  acov_var  <- stats::var(acov_est)
  acov_var1 <- stats::var(acov_est1)
  acov_var2 <- stats::var(acov_est2)
  acov_var3 <- stats::var(acov_est3)
  acov_var4 <- stats::var(acov_est4)

  acor_var  <- stats::var(acor_est)
  acor_var1 <- stats::var(acor_est1)
  acor_var2 <- stats::var(acor_est2)
  acor_var3 <- stats::var(acor_est3)
  acor_var4 <- stats::var(acor_est4)

  # Correlations
  mean_acov_cor  <- stats::cor(mean_est,  acov_est)
  mean_acov_cor1 <- stats::cor(mean_est1, acov_est1)
  mean_acov_cor2 <- stats::cor(mean_est2, acov_est2)
  mean_acov_cor3 <- stats::cor(mean_est3, acov_est3)
  mean_acov_cor4 <- stats::cor(mean_est4, acov_est4)

  mean_acor_cor  <- stats::cor(mean_est,  acor_est)
  mean_acor_cor1 <- stats::cor(mean_est1, acor_est1)
  mean_acor_cor2 <- stats::cor(mean_est2, acor_est2)
  mean_acor_cor3 <- stats::cor(mean_est3, acor_est3)
  mean_acor_cor4 <- stats::cor(mean_est4, acor_est4)

  acov_acor_cor  <- stats::cor(acov_est,  acor_est)
  acov_acor_cor1 <- stats::cor(acov_est1, acor_est1)
  acov_acor_cor2 <- stats::cor(acov_est2, acor_est2)
  acov_acor_cor3 <- stats::cor(acov_est3, acor_est3)
  acov_acor_cor4 <- stats::cor(acov_est4, acor_est4)

  # Estimates
  estimate  <- c(mean_mean,
                 acov_mean,
                 acor_mean,
                 mean_var,
                 acov_var,
                 acor_var,
                 mean_acov_cor,
                 mean_acor_cor,
                 acov_acor_cor)

  estimate1 <- c(mean_mean1,
                 acov_mean1,
                 acor_mean1,
                 mean_var1,
                 acov_var1,
                 acor_var1,
                 mean_acov_cor1,
                 mean_acor_cor1,
                 acov_acor_cor1)

  estimate2 <- c(mean_mean2,
                 acov_mean2,
                 acor_mean2,
                 mean_var2,
                 acov_var2,
                 acor_var2,
                 mean_acov_cor2,
                 mean_acor_cor2,
                 acov_acor_cor2)

  estimate3 <- c(mean_mean3,
                 acov_mean3,
                 acor_mean3,
                 mean_var3,
                 acov_var3,
                 acor_var3,
                 mean_acov_cor3,
                 mean_acor_cor3,
                 acov_acor_cor3)

  estimate4 <- c(mean_mean4,
                 acov_mean4,
                 acor_mean4,
                 mean_var4,
                 acov_var4,
                 acor_var4,
                 mean_acov_cor4,
                 mean_acor_cor4,
                 acov_acor_cor4)

  # HPJ bias-corrected estimate
  hpjestimate <- 2 * estimate - (estimate1 + estimate2 + estimate3 + estimate4) / 4

  return(hpjestimate)

}