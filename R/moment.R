#' Estimation of moments for heterogeneity in panel data without bias-correction
#'
#' \code{nemoment} implements the naive estimation of the moments
#' for the heterogeneous mean, autocovariance, and autocorrelation.
#' The procedure is proposed in Okui and Yanagi (2017).
#' See the package vignette via \code{vignette("panelhetero")} for details.
#'
#' @param data matrix of panel data in which each row is individual time series
#' @param acov_order non-negative integer for the order of the autocovariance
#' @param acor_order positive integer for the order of the autocorrelation
#' @param R positive integer for the number of bootstrap replications
#'
#' @importFrom boot boot
#' @importFrom stats na.omit quantile
#'
#' @return list that contains the following elements.
#' \item{estimate}{vector of estimates of the parameters}
#' \item{se}{vector of standard errors of the estimators}
#' \item{ci}{matrix of 95 percent confidence intervals of the parameters}
#' \item{quantity}{matrix that contains the estimated quantities}
#' \item{acov_order}{the same as the argument}
#' \item{acor_order}{the same as the argument}
#' \item{N}{the number of cross-sectional units}
#' \item{S}{the length of time series}
#' \item{R}{the number of bootstrap replications}
#'
#' @export
#'
nemoment <- function(data, acov_order = 0, acor_order = 1, R = 1000) {

  # handling errors
  error1(data = data, acov_order = acov_order, acor_order = acor_order, R = R)

  # omitting NA
  data <- na.omit(data)

  # sample sizes
  N <- nrow(data)
  S <- ncol(data)

  # estimated means, autocovariances, autocovariances
  mean_est <- rowMeans(data)
  acov_est <- apply(data, MARGIN = 1, acov, acov_order = acov_order)
  acor_est <- apply(data, MARGIN = 1, acor, acor_order = acor_order)
  equantity <- cbind(mean_est, acov_est, acor_est)

  # naive estimation with bootstrap
  bootstrap <- boot(data = equantity, statistic = momentest, R = R)

  # estimate
  estimate <- bootstrap$t0
  names(estimate) <- c("E(mean)", "E(acov)", "E(acor)",
                       "var(mean)", "var(acov)", "var(acor)",
                       "cor(mean, acov)", "cor(mean, acor)", "cor(acov, acor)")

  # standard errors
  temp <- t(t(bootstrap$t) - bootstrap$t0)
  se <- sqrt(colMeans(temp * temp))
  names(se) <- c("E(mean)", "E(acov)", "E(acor)",
                 "var(mean)", "var(acov)", "var(acor)",
                 "cor(mean, acov)", "cor(mean, acor)", "cor(acov, acor)")

  # confidence intervals
  quantiles <- apply(temp, MARGIN = 2, quantile, probs = c(0.025, 0.975))
  ci <- cbind(estimate + quantiles[1, ], estimate + quantiles[2, ])
  rownames(ci) <- c("E(mean)", "E(acov)", "E(acor)",
                    "var(mean)", "var(acov)", "var(acor)",
                    "cor(mean, acov)", "cor(mean, acor)", "cor(acov, acor)")
  colnames(ci) <- c("95% CI lower", "95% CI upper")

  # result
  quantity <- cbind(mean_est, acov_est, acor_est)
  colnames(quantity) <- c("mean", "autocovariance", "autocorrelation")

  result <- list(estimate = estimate, se = se, ci = ci,
                 quantity = quantity, acov_order = acov_order,
                 acor_order = acor_order, N = N, S = S, R = R)

  return(result)

}


#' computing estimates for the moments
#'
#' @param quantity N * 3 matrix of the means, autocovariances, autocorrelations
#' @param indices indices for bootstrap replications
#' @importFrom stats var cor
#'
momentest <- function(quantity, indices) {

  # estimated quantities
  mean_est <- quantity[indices, 1]
  acov_est <- quantity[indices, 2]
  acor_est <- quantity[indices, 3]

  # means
  mean_mean <- mean(mean_est)
  acov_mean <- mean(acov_est)
  acor_mean <- mean(acor_est)

  # variances
  mean_var <- var(mean_est)
  acov_var <- var(acov_est)
  acor_var <- var(acor_est)

  # correlations
  mean_acov_cor <- cor(mean_est, acov_est)
  mean_acor_cor <- cor(mean_est, acor_est)
  acov_acor_cor <- cor(acov_est, acor_est)

  # estimates
  estimate <- c(mean_mean, acov_mean, acor_mean, mean_var, acov_var, acor_var, mean_acov_cor, mean_acor_cor, acov_acor_cor)

  return(estimate)

}
