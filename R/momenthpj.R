#' HPJ bias-corrected estimation of moments for heterogeneity in panel data
#'
#' \code{hpjmoment} implements the HPJ bias-corrected estimation of the moments
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
hpjmoment <- function(data, acov_order = 0, acor_order = 1, R = 1000) {

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

  # HPJ bias-correction
  if (S %% 2 == 0) {

    # half-panel data for even T
    data1 <- data[, 1:(S / 2)]
    data2 <- data[, (S / 2 + 1):S]

    # estimated quangtities
    mean_est1 <- rowMeans(data1)
    mean_est2 <- rowMeans(data2)

    acov_est1 <- apply(data1, MARGIN = 1, acov, acov_order = acov_order)
    acov_est2 <- apply(data2, MARGIN = 1, acov, acov_order = acov_order)

    acor_est1 <- apply(data1, MARGIN = 1, acor, acor_order = acor_order)
    acor_est2 <- apply(data2, MARGIN = 1, acor, acor_order = acor_order)

    equantity <- cbind(mean_est, mean_est1, mean_est2,
                       acov_est, acov_est1, acov_est2,
                       acor_est, acor_est1, acor_est2)

    # estimation with bootstrap
    bootstrap <- boot(data = equantity, statistic = hpjmomentest1, R = R)

  } else {

    # half-panel data for odd T
    data1 <- data[, 1:floor(S / 2)]
    data2 <- data[, (floor(S / 2) + 1):S]
    data3 <- data[, 1:ceiling(S / 2)]
    data4 <- data[, (ceiling(S / 2) + 1):S]

    # estimated quantities
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

    equantity <- cbind(mean_est, mean_est1, mean_est2, mean_est3, mean_est4,
                       acov_est, acov_est1, acov_est2, acov_est3, acov_est4,
                       acor_est, acor_est1, acor_est2, acor_est3, acor_est4)

    # estimation with bootstrap
    bootstrap <- boot(data = equantity, statistic = hpjmomentest2, R = R)

  }

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

  # results
  quantity <- cbind(mean_est, acov_est, acor_est)
  colnames(quantity) <- c("mean", "autocovariance", "autocorrelation")

  result <- list(estimate = estimate, se = se, ci = ci, quantity = quantity,
                 acov_order = acov_order, acor_order = acor_order,
                 N = N, S = S, R = R)

  return(result)

}


#' computing HPJ bias-corrected estimates of the moments for even T
#'
#' @param quantity N * 3 matrix of the means, autocovariances, autocorrelations
#' @param indices indices for bootstrap replications
#'
hpjmomentest1 <- function(quantity, indices) {

  # estimated quantities
  mean_est <- quantity[indices, 1]
  mean_est1 <- quantity[indices, 2]
  mean_est2 <- quantity[indices, 3]

  acov_est <- quantity[indices, 4]
  acov_est1 <- quantity[indices, 5]
  acov_est2 <- quantity[indices, 6]

  acor_est <- quantity[indices, 7]
  acor_est1 <- quantity[indices, 8]
  acor_est2 <- quantity[indices, 9]

  # the means
  mean_mean <- mean(mean_est)
  mean_mean1 <- mean(mean_est1)
  mean_mean2 <- mean(mean_est2)

  acov_mean <- mean(acov_est)
  acov_mean1 <- mean(acov_est1)
  acov_mean2 <- mean(acov_est2)

  acor_mean <- mean(acor_est)
  acor_mean1 <- mean(acor_est1)
  acor_mean2 <- mean(acor_est2)

  # the variances
  mean_var <- var(mean_est)
  mean_var1 <- var(mean_est1)
  mean_var2 <- var(mean_est2)

  acov_var <- var(acov_est)
  acov_var1 <- var(acov_est1)
  acov_var2 <- var(acov_est2)

  acor_var <- var(acor_est)
  acor_var1 <- var(acor_est1)
  acor_var2 <- var(acor_est2)

  # the correlations
  mean_acov_cor <- cor(mean_est, acov_est)
  mean_acov_cor1 <- cor(mean_est1, acov_est1)
  mean_acov_cor2 <- cor(mean_est2, acov_est2)

  mean_acor_cor <- cor(mean_est, acor_est)
  mean_acor_cor1 <- cor(mean_est1, acor_est1)
  mean_acor_cor2 <- cor(mean_est2, acor_est2)

  acov_acor_cor <- cor(acov_est, acor_est)
  acov_acor_cor1 <- cor(acov_est1, acor_est1)
  acov_acor_cor2 <- cor(acov_est2, acor_est2)

  # estimates
  estimate <- c(mean_mean, acov_mean, acor_mean, mean_var, acov_var, acor_var, mean_acov_cor, mean_acor_cor, acov_acor_cor)
  estimate1 <- c(mean_mean1, acov_mean1, acor_mean1, mean_var1, acov_var1, acor_var1, mean_acov_cor1, mean_acor_cor1, acov_acor_cor1)
  estimate2 <- c(mean_mean2, acov_mean2, acor_mean2, mean_var2, acov_var2, acor_var2, mean_acov_cor2, mean_acor_cor2, acov_acor_cor2)

  # HPJ estimate
  hpjestimate <- 2 * estimate - (estimate1 + estimate2) / 2

  return(hpjestimate)

}


#' computing HPJ bias-corrected estimates of the moments for odd T
#'
#' @param quantity N * 3 matrix of the estimated means, autocovariances, autocorrelations
#' @param indices indices for bootstrap replications
#'
#' @importFrom stats var cor
#'
hpjmomentest2 <- function(quantity, indices) {

  # estimated quantities
  mean_est <- quantity[indices, 1]
  mean_est1 <- quantity[indices, 2]
  mean_est2 <- quantity[indices, 3]
  mean_est3 <- quantity[indices, 4]
  mean_est4 <- quantity[indices, 5]

  acov_est <- quantity[indices, 6]
  acov_est1 <- quantity[indices, 7]
  acov_est2 <- quantity[indices, 8]
  acov_est3 <- quantity[indices, 9]
  acov_est4 <- quantity[indices, 10]

  acor_est <- quantity[indices, 11]
  acor_est1 <- quantity[indices, 12]
  acor_est2 <- quantity[indices, 13]
  acor_est3 <- quantity[indices, 14]
  acor_est4 <- quantity[indices, 15]

  # the means
  mean_mean <- mean(mean_est)
  mean_mean1 <- mean(mean_est1)
  mean_mean2 <- mean(mean_est2)
  mean_mean3 <- mean(mean_est3)
  mean_mean4 <- mean(mean_est4)

  acov_mean <- mean(acov_est)
  acov_mean1 <- mean(acov_est1)
  acov_mean2 <- mean(acov_est2)
  acov_mean3 <- mean(acov_est3)
  acov_mean4 <- mean(acov_est4)

  acor_mean <- mean(acor_est)
  acor_mean1 <- mean(acor_est1)
  acor_mean2 <- mean(acor_est2)
  acor_mean3 <- mean(acor_est3)
  acor_mean4 <- mean(acor_est4)

  # the variances
  mean_var <- var(mean_est)
  mean_var1 <- var(mean_est1)
  mean_var2 <- var(mean_est2)
  mean_var3 <- var(mean_est3)
  mean_var4 <- var(mean_est4)

  acov_var <- var(acov_est)
  acov_var1 <- var(acov_est1)
  acov_var2 <- var(acov_est2)
  acov_var3 <- var(acov_est3)
  acov_var4 <- var(acov_est4)

  acor_var <- var(acor_est)
  acor_var1 <- var(acor_est1)
  acor_var2 <- var(acor_est2)
  acor_var3 <- var(acor_est3)
  acor_var4 <- var(acor_est4)

  # the correlations
  mean_acov_cor <- cor(mean_est, acov_est)
  mean_acov_cor1 <- cor(mean_est1, acov_est1)
  mean_acov_cor2 <- cor(mean_est2, acov_est2)
  mean_acov_cor3 <- cor(mean_est3, acov_est3)
  mean_acov_cor4 <- cor(mean_est4, acov_est4)

  mean_acor_cor <- cor(mean_est, acor_est)
  mean_acor_cor1 <- cor(mean_est1, acor_est1)
  mean_acor_cor2 <- cor(mean_est2, acor_est2)
  mean_acor_cor3 <- cor(mean_est3, acor_est3)
  mean_acor_cor4 <- cor(mean_est4, acor_est4)

  acov_acor_cor <- cor(acov_est, acor_est)
  acov_acor_cor1 <- cor(acov_est1, acor_est1)
  acov_acor_cor2 <- cor(acov_est2, acor_est2)
  acov_acor_cor3 <- cor(acov_est3, acor_est3)
  acov_acor_cor4 <- cor(acov_est4, acor_est4)

  # estimates
  estimate <- c(mean_mean, acov_mean, acor_mean, mean_var, acov_var, acor_var, mean_acov_cor, mean_acor_cor, acov_acor_cor)
  estimate1 <- c(mean_mean1, acov_mean1, acor_mean1, mean_var1, acov_var1, acor_var1, mean_acov_cor1, mean_acor_cor1, acov_acor_cor1)
  estimate2 <- c(mean_mean2, acov_mean2, acor_mean2, mean_var2, acov_var2, acor_var2, mean_acov_cor2, mean_acor_cor2, acov_acor_cor2)
  estimate3 <- c(mean_mean3, acov_mean3, acor_mean3, mean_var3, acov_var3, acor_var3, mean_acov_cor3, mean_acor_cor3, acov_acor_cor3)
  estimate4 <- c(mean_mean4, acov_mean4, acor_mean4, mean_var4, acov_var4, acor_var4, mean_acov_cor4, mean_acor_cor4, acov_acor_cor4)

  # HPJ estimate
  hpjestimate <- 2 * estimate - (estimate1 + estimate2 + estimate3 + estimate4) / 4

  return(hpjestimate)

}
