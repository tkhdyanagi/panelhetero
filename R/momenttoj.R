#' TOJ bias-corrected estimation of moments for heterogeneity in panel data
#'
#' \code{tojmoment} implements the TOJ bias-corrected estimation of the moments
#' for the heterogeneous mean, autocovariance, and autocorrelation.
#' The procedure is proposed in Okui and Yanagi (2019).
#' See the package vignette via \code{vignette("panelhetero")} for details.
#'
#' @param data matrix of panel data in which each row is individual time series
#' @param acov_order non-negative integer for the order of the autocovariance
#' @param acor_order positive integer for the order of the autocorrelation
#' @param R positive integer for the number of bootstrap replications
#'
#' @importFrom boot boot
#' @importFrom stats na.omit quantile sd
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
tojmoment <- function(data, acov_order = 0, acor_order = 1, R = 1000) {

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

  # TOJ bias-correction
  if (S %% 6 == 0) {

    # split panel data for T equivalent to 0 modulo 6
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

    equantity <- cbind(mean_est, mean_est21, mean_est22, mean_est31, mean_est32, mean_est33,
                       acov_est, acov_est21, acov_est22, acov_est31, acov_est32, acov_est33,
                       acor_est, acor_est21, acor_est22, acor_est31, acor_est32, acor_est33)

    # estimation with bootstrap
    bootstrap <- boot(data = equantity, statistic = tojmomentest0, R = R)

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

    equantity <- cbind(mean_est, mean_est21, mean_est22, mean_est23, mean_est24, mean_est31, mean_est32, mean_est33, mean_est34, mean_est35, mean_est36, mean_est37, mean_est38, mean_est39,
                       acov_est, acov_est21, acov_est22, acov_est23, acov_est24, acov_est31, acov_est32, acov_est33, acov_est34, acov_est35, acov_est36, acov_est37, acov_est38, acov_est39,
                       acor_est, acor_est21, acor_est22, acor_est23, acor_est24, acor_est31, acor_est32, acor_est33, acor_est34, acor_est35, acor_est36, acor_est37, acor_est38, acor_est39)

    # estimation with bootstrap
    bootstrap <- boot(data = equantity, statistic = tojmomentest1, R = R)


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

    equantity <- cbind(mean_est, mean_est21, mean_est22, mean_est31, mean_est32, mean_est33, mean_est34, mean_est35, mean_est36, mean_est37, mean_est38, mean_est39,
                       acov_est, acov_est21, acov_est22, acov_est31, acov_est32, acov_est33, acov_est34, acov_est35, acov_est36, acov_est37, acov_est38, acov_est39,
                       acor_est, acor_est21, acor_est22, acor_est31, acor_est32, acor_est33, acor_est34, acor_est35, acor_est36, acor_est37, acor_est38, acor_est39)

    # estimation with bootstrap
    bootstrap <- boot(data = equantity, statistic = tojmomentest2, R = R)


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

    equantity <- cbind(mean_est, mean_est21, mean_est22, mean_est23, mean_est24, mean_est31, mean_est32, mean_est33,
                       acov_est, acov_est21, acov_est22, acov_est23, acov_est24, acov_est31, acov_est32, acov_est33,
                       acor_est, acor_est21, acor_est22, acor_est23, acor_est24, acor_est31, acor_est32, acor_est33)

    # estimation with bootstrap
    bootstrap <- boot(data = equantity, statistic = tojmomentest3, R = R)


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

    equantity <- cbind(mean_est, mean_est21, mean_est22, mean_est31, mean_est32, mean_est33, mean_est34, mean_est35, mean_est36, mean_est37, mean_est38, mean_est39,
                       acov_est, acov_est21, acov_est22, acov_est31, acov_est32, acov_est33, acov_est34, acov_est35, acov_est36, acov_est37, acov_est38, acov_est39,
                       acor_est, acor_est21, acor_est22, acor_est31, acor_est32, acor_est33, acor_est34, acor_est35, acor_est36, acor_est37, acor_est38, acor_est39)

    # estimation with bootstrap
    bootstrap <- boot(data = equantity, statistic = tojmomentest4, R = R)


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

    equantity <- cbind(mean_est, mean_est21, mean_est22, mean_est23, mean_est24, mean_est31, mean_est32, mean_est33, mean_est34, mean_est35, mean_est36, mean_est37, mean_est38, mean_est39,
                       acov_est, acov_est21, acov_est22, acov_est23, acov_est24, acov_est31, acov_est32, acov_est33, acov_est34, acov_est35, acov_est36, acov_est37, acov_est38, acov_est39,
                       acor_est, acor_est21, acor_est22, acor_est23, acor_est24, acor_est31, acor_est32, acor_est33, acor_est34, acor_est35, acor_est36, acor_est37, acor_est38, acor_est39)

    # estimation with bootstrap
    bootstrap <- boot(data = equantity, statistic = tojmomentest5, R = R)

  }


  # estimate
  estimate <- bootstrap$t0
  names(estimate) <- c("E(mean)", "E(acov)", "E(acor)",
                       "var(mean)", "var(acov)", "var(acor)",
                       "cor(mean, acov)", "cor(mean, acor)", "cor(acov, acor)")

  # standard errors
  se <- apply(bootstrap$t, MARGIN = 2, sd)
  names(se) <- c("E(mean)", "E(acov)", "E(acor)",
                 "var(mean)", "var(acov)", "var(acor)",
                 "cor(mean, acov)", "cor(mean, acor)", "cor(acov, acor)")

  # confidence intervals
  temp <- t(t(bootstrap$t) - bootstrap$t0)
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


#' computing TOJ bias-corrected estimate for T equivalent to 0 modulo 6
#'
#' @param quantity N * (3 * 6) matrix of the means, autocovariances, autocorrelations
#' @param indices indices for bootstrap replications
#'
tojmomentest0 <- function(quantity, indices) {

  # estimated quantities
  mean_est <- quantity[indices, 1]
  mean_est21 <- quantity[indices, 2]
  mean_est22 <- quantity[indices, 3]
  mean_est31 <- quantity[indices, 4]
  mean_est32 <- quantity[indices, 5]
  mean_est33 <- quantity[indices, 6]

  acov_est <- quantity[indices, 7]
  acov_est21 <- quantity[indices, 8]
  acov_est22 <- quantity[indices, 9]
  acov_est31 <- quantity[indices, 10]
  acov_est32 <- quantity[indices, 11]
  acov_est33 <- quantity[indices, 12]

  acor_est <- quantity[indices, 13]
  acor_est21 <- quantity[indices, 14]
  acor_est22 <- quantity[indices, 15]
  acor_est31 <- quantity[indices, 16]
  acor_est32 <- quantity[indices, 17]
  acor_est33 <- quantity[indices, 18]

  # the means
  mean_mean <- mean(mean_est)
  mean_mean21 <- mean(mean_est21)
  mean_mean22 <- mean(mean_est22)
  mean_mean31 <- mean(mean_est31)
  mean_mean32 <- mean(mean_est32)
  mean_mean33 <- mean(mean_est33)

  acov_mean <- mean(acov_est)
  acov_mean21 <- mean(acov_est21)
  acov_mean22 <- mean(acov_est22)
  acov_mean31 <- mean(acov_est31)
  acov_mean32 <- mean(acov_est32)
  acov_mean33 <- mean(acov_est33)

  acor_mean <- mean(acor_est)
  acor_mean21 <- mean(acor_est21)
  acor_mean22 <- mean(acor_est22)
  acor_mean31 <- mean(acor_est31)
  acor_mean32 <- mean(acor_est32)
  acor_mean33 <- mean(acor_est33)

  # the variances
  mean_var <- var(mean_est)
  mean_var21 <- var(mean_est21)
  mean_var22 <- var(mean_est22)
  mean_var31 <- var(mean_est31)
  mean_var32 <- var(mean_est32)
  mean_var33 <- var(mean_est33)

  acov_var <- var(acov_est)
  acov_var21 <- var(acov_est21)
  acov_var22 <- var(acov_est22)
  acov_var31 <- var(acov_est31)
  acov_var32 <- var(acov_est32)
  acov_var33 <- var(acov_est33)

  acor_var <- var(acor_est)
  acor_var21 <- var(acor_est21)
  acor_var22 <- var(acor_est22)
  acor_var31 <- var(acor_est31)
  acor_var32 <- var(acor_est32)
  acor_var33 <- var(acor_est33)

  # the correlations
  mean_acov_cor <- cor(mean_est, acov_est)
  mean_acov_cor21 <- cor(mean_est21, acov_est21)
  mean_acov_cor22 <- cor(mean_est22, acov_est22)
  mean_acov_cor31 <- cor(mean_est31, acov_est31)
  mean_acov_cor32 <- cor(mean_est32, acov_est32)
  mean_acov_cor33 <- cor(mean_est33, acov_est33)

  mean_acor_cor <- cor(mean_est, acor_est)
  mean_acor_cor21 <- cor(mean_est21, acor_est21)
  mean_acor_cor22 <- cor(mean_est22, acor_est22)
  mean_acor_cor31 <- cor(mean_est31, acor_est31)
  mean_acor_cor32 <- cor(mean_est32, acor_est32)
  mean_acor_cor33 <- cor(mean_est33, acor_est33)

  acov_acor_cor <- cor(acov_est, acor_est)
  acov_acor_cor21 <- cor(acov_est21, acor_est21)
  acov_acor_cor22 <- cor(acov_est22, acor_est22)
  acov_acor_cor31 <- cor(acov_est31, acor_est31)
  acov_acor_cor32 <- cor(acov_est32, acor_est32)
  acov_acor_cor33 <- cor(acov_est33, acor_est33)

  # estimates
  estimate <- c(mean_mean, acov_mean, acor_mean, mean_var, acov_var, acor_var, mean_acov_cor, mean_acor_cor, acov_acor_cor)
  estimate21 <- c(mean_mean21, acov_mean21, acor_mean21, mean_var21, acov_var21, acor_var21, mean_acov_cor21, mean_acor_cor21, acov_acor_cor21)
  estimate22 <- c(mean_mean22, acov_mean22, acor_mean22, mean_var22, acov_var22, acor_var22, mean_acov_cor22, mean_acor_cor22, acov_acor_cor22)
  estimate31 <- c(mean_mean31, acov_mean31, acor_mean31, mean_var31, acov_var31, acor_var31, mean_acov_cor31, mean_acor_cor31, acov_acor_cor31)
  estimate32 <- c(mean_mean32, acov_mean32, acor_mean32, mean_var32, acov_var32, acor_var32, mean_acov_cor32, mean_acor_cor32, acov_acor_cor32)
  estimate33 <- c(mean_mean33, acov_mean33, acor_mean33, mean_var33, acov_var33, acor_var33, mean_acov_cor33, mean_acor_cor33, acov_acor_cor33)

  # TOJ estimate
  tojestimate <- 3.536 * estimate - 4.072 * (estimate21 + estimate22) / 2 + 1.536 * (estimate31 + estimate32 + estimate33) / 3

  return(tojestimate)


}

#' computing TOJ bias-corrected estimate for T equivalent to 1 modulo 6
#'
#' @param quantity N * (3 * 14) matrix of the means, autocovariances, autocorrelations
#' @param indices indices for bootstrap replications
#'
tojmomentest1 <- function(quantity, indices) {

  # estimated quantities
  mean_est <- quantity[indices, 1]
  mean_est21 <- quantity[indices, 2]
  mean_est22 <- quantity[indices, 3]
  mean_est23 <- quantity[indices, 4]
  mean_est24 <- quantity[indices, 5]
  mean_est31 <- quantity[indices, 6]
  mean_est32 <- quantity[indices, 7]
  mean_est33 <- quantity[indices, 8]
  mean_est34 <- quantity[indices, 9]
  mean_est35 <- quantity[indices, 10]
  mean_est36 <- quantity[indices, 11]
  mean_est37 <- quantity[indices, 12]
  mean_est38 <- quantity[indices, 13]
  mean_est39 <- quantity[indices, 14]

  acov_est <- quantity[indices, 15]
  acov_est21 <- quantity[indices, 16]
  acov_est22 <- quantity[indices, 17]
  acov_est23 <- quantity[indices, 18]
  acov_est24 <- quantity[indices, 19]
  acov_est31 <- quantity[indices, 20]
  acov_est32 <- quantity[indices, 21]
  acov_est33 <- quantity[indices, 22]
  acov_est34 <- quantity[indices, 23]
  acov_est35 <- quantity[indices, 24]
  acov_est36 <- quantity[indices, 25]
  acov_est37 <- quantity[indices, 26]
  acov_est38 <- quantity[indices, 27]
  acov_est39 <- quantity[indices, 28]

  acor_est <- quantity[indices, 29]
  acor_est21 <- quantity[indices, 30]
  acor_est22 <- quantity[indices, 31]
  acor_est23 <- quantity[indices, 32]
  acor_est24 <- quantity[indices, 33]
  acor_est31 <- quantity[indices, 34]
  acor_est32 <- quantity[indices, 35]
  acor_est33 <- quantity[indices, 36]
  acor_est34 <- quantity[indices, 37]
  acor_est35 <- quantity[indices, 38]
  acor_est36 <- quantity[indices, 39]
  acor_est37 <- quantity[indices, 40]
  acor_est38 <- quantity[indices, 41]
  acor_est39 <- quantity[indices, 42]

  # the means
  mean_mean <- mean(mean_est)
  mean_mean21 <- mean(mean_est21)
  mean_mean22 <- mean(mean_est22)
  mean_mean23 <- mean(mean_est23)
  mean_mean24 <- mean(mean_est24)
  mean_mean31 <- mean(mean_est31)
  mean_mean32 <- mean(mean_est32)
  mean_mean33 <- mean(mean_est33)
  mean_mean34 <- mean(mean_est34)
  mean_mean35 <- mean(mean_est35)
  mean_mean36 <- mean(mean_est36)
  mean_mean37 <- mean(mean_est37)
  mean_mean38 <- mean(mean_est38)
  mean_mean39 <- mean(mean_est39)

  acov_mean <- mean(acov_est)
  acov_mean21 <- mean(acov_est21)
  acov_mean22 <- mean(acov_est22)
  acov_mean23 <- mean(acov_est23)
  acov_mean24 <- mean(acov_est24)
  acov_mean31 <- mean(acov_est31)
  acov_mean32 <- mean(acov_est32)
  acov_mean33 <- mean(acov_est33)
  acov_mean34 <- mean(acov_est34)
  acov_mean35 <- mean(acov_est35)
  acov_mean36 <- mean(acov_est36)
  acov_mean37 <- mean(acov_est37)
  acov_mean38 <- mean(acov_est38)
  acov_mean39 <- mean(acov_est39)

  acor_mean <- mean(acor_est)
  acor_mean21 <- mean(acor_est21)
  acor_mean22 <- mean(acor_est22)
  acor_mean23 <- mean(acor_est23)
  acor_mean24 <- mean(acor_est24)
  acor_mean31 <- mean(acor_est31)
  acor_mean32 <- mean(acor_est32)
  acor_mean33 <- mean(acor_est33)
  acor_mean34 <- mean(acor_est34)
  acor_mean35 <- mean(acor_est35)
  acor_mean36 <- mean(acor_est36)
  acor_mean37 <- mean(acor_est37)
  acor_mean38 <- mean(acor_est38)
  acor_mean39 <- mean(acor_est39)

  # the variances
  mean_var <- var(mean_est)
  mean_var21 <- var(mean_est21)
  mean_var22 <- var(mean_est22)
  mean_var23 <- var(mean_est23)
  mean_var24 <- var(mean_est24)
  mean_var31 <- var(mean_est31)
  mean_var32 <- var(mean_est32)
  mean_var33 <- var(mean_est33)
  mean_var34 <- var(mean_est34)
  mean_var35 <- var(mean_est35)
  mean_var36 <- var(mean_est36)
  mean_var37 <- var(mean_est37)
  mean_var38 <- var(mean_est38)
  mean_var39 <- var(mean_est39)

  acov_var <- var(acov_est)
  acov_var21 <- var(acov_est21)
  acov_var22 <- var(acov_est22)
  acov_var23 <- var(acov_est23)
  acov_var24 <- var(acov_est24)
  acov_var31 <- var(acov_est31)
  acov_var32 <- var(acov_est32)
  acov_var33 <- var(acov_est33)
  acov_var34 <- var(acov_est34)
  acov_var35 <- var(acov_est35)
  acov_var36 <- var(acov_est36)
  acov_var37 <- var(acov_est37)
  acov_var38 <- var(acov_est38)
  acov_var39 <- var(acov_est39)

  acor_var <- var(acor_est)
  acor_var21 <- var(acor_est21)
  acor_var22 <- var(acor_est22)
  acor_var23 <- var(acor_est23)
  acor_var24 <- var(acor_est24)
  acor_var31 <- var(acor_est31)
  acor_var32 <- var(acor_est32)
  acor_var33 <- var(acor_est33)
  acor_var34 <- var(acor_est34)
  acor_var35 <- var(acor_est35)
  acor_var36 <- var(acor_est36)
  acor_var37 <- var(acor_est37)
  acor_var38 <- var(acor_est38)
  acor_var39 <- var(acor_est39)

  # the correlations
  mean_acov_cor <- cor(mean_est, acov_est)
  mean_acov_cor21 <- cor(mean_est21, acov_est21)
  mean_acov_cor22 <- cor(mean_est22, acov_est22)
  mean_acov_cor23 <- cor(mean_est23, acov_est23)
  mean_acov_cor24 <- cor(mean_est24, acov_est24)
  mean_acov_cor31 <- cor(mean_est31, acov_est31)
  mean_acov_cor32 <- cor(mean_est32, acov_est32)
  mean_acov_cor33 <- cor(mean_est33, acov_est33)
  mean_acov_cor34 <- cor(mean_est34, acov_est34)
  mean_acov_cor35 <- cor(mean_est35, acov_est35)
  mean_acov_cor36 <- cor(mean_est36, acov_est36)
  mean_acov_cor37 <- cor(mean_est37, acov_est37)
  mean_acov_cor38 <- cor(mean_est38, acov_est38)
  mean_acov_cor39 <- cor(mean_est39, acov_est39)

  mean_acor_cor <- cor(mean_est, acor_est)
  mean_acor_cor21 <- cor(mean_est21, acor_est21)
  mean_acor_cor22 <- cor(mean_est22, acor_est22)
  mean_acor_cor23 <- cor(mean_est23, acor_est23)
  mean_acor_cor24 <- cor(mean_est24, acor_est24)
  mean_acor_cor31 <- cor(mean_est31, acor_est31)
  mean_acor_cor32 <- cor(mean_est32, acor_est32)
  mean_acor_cor33 <- cor(mean_est33, acor_est33)
  mean_acor_cor34 <- cor(mean_est34, acor_est34)
  mean_acor_cor35 <- cor(mean_est35, acor_est35)
  mean_acor_cor36 <- cor(mean_est36, acor_est36)
  mean_acor_cor37 <- cor(mean_est37, acor_est37)
  mean_acor_cor38 <- cor(mean_est38, acor_est38)
  mean_acor_cor39 <- cor(mean_est39, acor_est39)

  acov_acor_cor <- cor(acov_est, acor_est)
  acov_acor_cor21 <- cor(acov_est21, acor_est21)
  acov_acor_cor22 <- cor(acov_est22, acor_est22)
  acov_acor_cor23 <- cor(acov_est23, acor_est23)
  acov_acor_cor24 <- cor(acov_est24, acor_est24)
  acov_acor_cor31 <- cor(acov_est31, acor_est31)
  acov_acor_cor32 <- cor(acov_est32, acor_est32)
  acov_acor_cor33 <- cor(acov_est33, acor_est33)
  acov_acor_cor34 <- cor(acov_est34, acor_est34)
  acov_acor_cor35 <- cor(acov_est35, acor_est35)
  acov_acor_cor36 <- cor(acov_est36, acor_est36)
  acov_acor_cor37 <- cor(acov_est37, acor_est37)
  acov_acor_cor38 <- cor(acov_est38, acor_est38)
  acov_acor_cor39 <- cor(acov_est39, acor_est39)

  # estimates
  estimate <- c(mean_mean, acov_mean, acor_mean, mean_var, acov_var, acor_var, mean_acov_cor, mean_acor_cor, acov_acor_cor)
  estimate21 <- c(mean_mean21, acov_mean21, acor_mean21, mean_var21, acov_var21, acor_var21, mean_acov_cor21, mean_acor_cor21, acov_acor_cor21)
  estimate22 <- c(mean_mean22, acov_mean22, acor_mean22, mean_var22, acov_var22, acor_var22, mean_acov_cor22, mean_acor_cor22, acov_acor_cor22)
  estimate23 <- c(mean_mean23, acov_mean23, acor_mean23, mean_var23, acov_var23, acor_var23, mean_acov_cor23, mean_acor_cor23, acov_acor_cor23)
  estimate24 <- c(mean_mean24, acov_mean24, acor_mean24, mean_var24, acov_var24, acor_var24, mean_acov_cor24, mean_acor_cor24, acov_acor_cor24)
  estimate31 <- c(mean_mean31, acov_mean31, acor_mean31, mean_var31, acov_var31, acor_var31, mean_acov_cor31, mean_acor_cor31, acov_acor_cor31)
  estimate32 <- c(mean_mean32, acov_mean32, acor_mean32, mean_var32, acov_var32, acor_var32, mean_acov_cor32, mean_acor_cor32, acov_acor_cor32)
  estimate33 <- c(mean_mean33, acov_mean33, acor_mean33, mean_var33, acov_var33, acor_var33, mean_acov_cor33, mean_acor_cor33, acov_acor_cor33)
  estimate34 <- c(mean_mean34, acov_mean34, acor_mean34, mean_var34, acov_var34, acor_var34, mean_acov_cor34, mean_acor_cor34, acov_acor_cor34)
  estimate35 <- c(mean_mean35, acov_mean35, acor_mean35, mean_var35, acov_var35, acor_var35, mean_acov_cor35, mean_acor_cor35, acov_acor_cor35)
  estimate36 <- c(mean_mean36, acov_mean36, acor_mean36, mean_var36, acov_var36, acor_var36, mean_acov_cor36, mean_acor_cor36, acov_acor_cor36)
  estimate37 <- c(mean_mean37, acov_mean37, acor_mean37, mean_var37, acov_var37, acor_var37, mean_acov_cor37, mean_acor_cor37, acov_acor_cor37)
  estimate38 <- c(mean_mean38, acov_mean38, acor_mean38, mean_var38, acov_var38, acor_var38, mean_acov_cor38, mean_acor_cor38, acov_acor_cor38)
  estimate39 <- c(mean_mean39, acov_mean39, acor_mean39, mean_var39, acov_var39, acor_var39, mean_acov_cor39, mean_acor_cor39, acov_acor_cor39)

  # TOJ estimate
  tojestimate <- 3.536 * estimate - 4.072 * (estimate21 + estimate22 + estimate23 + estimate24) / 4 + 1.536 * (estimate31 + estimate32 + estimate33 + estimate34 + estimate35 + estimate36 + estimate37 + estimate38 + estimate39) / 9

  return(tojestimate)

}

#' computing TOJ bias-corrected estimate for T equivalent to 2 modulo 6
#'
#' @param quantity N * (3 * 12) matrix of the means, autocovariances, autocorrelations
#' @param indices indices for bootstrap replications
#'
tojmomentest2 <- function(quantity, indices) {

  # estimated quantities
  mean_est <- quantity[indices, 1]
  mean_est21 <- quantity[indices, 2]
  mean_est22 <- quantity[indices, 3]
  mean_est31 <- quantity[indices, 4]
  mean_est32 <- quantity[indices, 5]
  mean_est33 <- quantity[indices, 6]
  mean_est34 <- quantity[indices, 7]
  mean_est35 <- quantity[indices, 8]
  mean_est36 <- quantity[indices, 9]
  mean_est37 <- quantity[indices, 10]
  mean_est38 <- quantity[indices, 11]
  mean_est39 <- quantity[indices, 12]

  acov_est <- quantity[indices, 13]
  acov_est21 <- quantity[indices, 14]
  acov_est22 <- quantity[indices, 15]
  acov_est31 <- quantity[indices, 16]
  acov_est32 <- quantity[indices, 17]
  acov_est33 <- quantity[indices, 18]
  acov_est34 <- quantity[indices, 19]
  acov_est35 <- quantity[indices, 20]
  acov_est36 <- quantity[indices, 21]
  acov_est37 <- quantity[indices, 22]
  acov_est38 <- quantity[indices, 23]
  acov_est39 <- quantity[indices, 24]

  acor_est <- quantity[indices, 25]
  acor_est21 <- quantity[indices, 26]
  acor_est22 <- quantity[indices, 27]
  acor_est31 <- quantity[indices, 28]
  acor_est32 <- quantity[indices, 29]
  acor_est33 <- quantity[indices, 30]
  acor_est34 <- quantity[indices, 31]
  acor_est35 <- quantity[indices, 32]
  acor_est36 <- quantity[indices, 33]
  acor_est37 <- quantity[indices, 34]
  acor_est38 <- quantity[indices, 35]
  acor_est39 <- quantity[indices, 36]

  # the means
  mean_mean <- mean(mean_est)
  mean_mean21 <- mean(mean_est21)
  mean_mean22 <- mean(mean_est22)
  mean_mean31 <- mean(mean_est31)
  mean_mean32 <- mean(mean_est32)
  mean_mean33 <- mean(mean_est33)
  mean_mean34 <- mean(mean_est34)
  mean_mean35 <- mean(mean_est35)
  mean_mean36 <- mean(mean_est36)
  mean_mean37 <- mean(mean_est37)
  mean_mean38 <- mean(mean_est38)
  mean_mean39 <- mean(mean_est39)

  acov_mean <- mean(acov_est)
  acov_mean21 <- mean(acov_est21)
  acov_mean22 <- mean(acov_est22)
  acov_mean31 <- mean(acov_est31)
  acov_mean32 <- mean(acov_est32)
  acov_mean33 <- mean(acov_est33)
  acov_mean34 <- mean(acov_est34)
  acov_mean35 <- mean(acov_est35)
  acov_mean36 <- mean(acov_est36)
  acov_mean37 <- mean(acov_est37)
  acov_mean38 <- mean(acov_est38)
  acov_mean39 <- mean(acov_est39)

  acor_mean <- mean(acor_est)
  acor_mean21 <- mean(acor_est21)
  acor_mean22 <- mean(acor_est22)
  acor_mean31 <- mean(acor_est31)
  acor_mean32 <- mean(acor_est32)
  acor_mean33 <- mean(acor_est33)
  acor_mean34 <- mean(acor_est34)
  acor_mean35 <- mean(acor_est35)
  acor_mean36 <- mean(acor_est36)
  acor_mean37 <- mean(acor_est37)
  acor_mean38 <- mean(acor_est38)
  acor_mean39 <- mean(acor_est39)

  # the variances
  mean_var <- var(mean_est)
  mean_var21 <- var(mean_est21)
  mean_var22 <- var(mean_est22)
  mean_var31 <- var(mean_est31)
  mean_var32 <- var(mean_est32)
  mean_var33 <- var(mean_est33)
  mean_var34 <- var(mean_est34)
  mean_var35 <- var(mean_est35)
  mean_var36 <- var(mean_est36)
  mean_var37 <- var(mean_est37)
  mean_var38 <- var(mean_est38)
  mean_var39 <- var(mean_est39)

  acov_var <- var(acov_est)
  acov_var21 <- var(acov_est21)
  acov_var22 <- var(acov_est22)
  acov_var31 <- var(acov_est31)
  acov_var32 <- var(acov_est32)
  acov_var33 <- var(acov_est33)
  acov_var34 <- var(acov_est34)
  acov_var35 <- var(acov_est35)
  acov_var36 <- var(acov_est36)
  acov_var37 <- var(acov_est37)
  acov_var38 <- var(acov_est38)
  acov_var39 <- var(acov_est39)

  acor_var <- var(acor_est)
  acor_var21 <- var(acor_est21)
  acor_var22 <- var(acor_est22)
  acor_var31 <- var(acor_est31)
  acor_var32 <- var(acor_est32)
  acor_var33 <- var(acor_est33)
  acor_var34 <- var(acor_est34)
  acor_var35 <- var(acor_est35)
  acor_var36 <- var(acor_est36)
  acor_var37 <- var(acor_est37)
  acor_var38 <- var(acor_est38)
  acor_var39 <- var(acor_est39)

  # the correlations
  mean_acov_cor <- cor(mean_est, acov_est)
  mean_acov_cor21 <- cor(mean_est21, acov_est21)
  mean_acov_cor22 <- cor(mean_est22, acov_est22)
  mean_acov_cor31 <- cor(mean_est31, acov_est31)
  mean_acov_cor32 <- cor(mean_est32, acov_est32)
  mean_acov_cor33 <- cor(mean_est33, acov_est33)
  mean_acov_cor34 <- cor(mean_est34, acov_est34)
  mean_acov_cor35 <- cor(mean_est35, acov_est35)
  mean_acov_cor36 <- cor(mean_est36, acov_est36)
  mean_acov_cor37 <- cor(mean_est37, acov_est37)
  mean_acov_cor38 <- cor(mean_est38, acov_est38)
  mean_acov_cor39 <- cor(mean_est39, acov_est39)

  mean_acor_cor <- cor(mean_est, acor_est)
  mean_acor_cor21 <- cor(mean_est21, acor_est21)
  mean_acor_cor22 <- cor(mean_est22, acor_est22)
  mean_acor_cor31 <- cor(mean_est31, acor_est31)
  mean_acor_cor32 <- cor(mean_est32, acor_est32)
  mean_acor_cor33 <- cor(mean_est33, acor_est33)
  mean_acor_cor34 <- cor(mean_est34, acor_est34)
  mean_acor_cor35 <- cor(mean_est35, acor_est35)
  mean_acor_cor36 <- cor(mean_est36, acor_est36)
  mean_acor_cor37 <- cor(mean_est37, acor_est37)
  mean_acor_cor38 <- cor(mean_est38, acor_est38)
  mean_acor_cor39 <- cor(mean_est39, acor_est39)

  acov_acor_cor <- cor(acov_est, acor_est)
  acov_acor_cor21 <- cor(acov_est21, acor_est21)
  acov_acor_cor22 <- cor(acov_est22, acor_est22)
  acov_acor_cor31 <- cor(acov_est31, acor_est31)
  acov_acor_cor32 <- cor(acov_est32, acor_est32)
  acov_acor_cor33 <- cor(acov_est33, acor_est33)
  acov_acor_cor34 <- cor(acov_est34, acor_est34)
  acov_acor_cor35 <- cor(acov_est35, acor_est35)
  acov_acor_cor36 <- cor(acov_est36, acor_est36)
  acov_acor_cor37 <- cor(acov_est37, acor_est37)
  acov_acor_cor38 <- cor(acov_est38, acor_est38)
  acov_acor_cor39 <- cor(acov_est39, acor_est39)

  # estimates
  estimate <- c(mean_mean, acov_mean, acor_mean, mean_var, acov_var, acor_var, mean_acov_cor, mean_acor_cor, acov_acor_cor)
  estimate21 <- c(mean_mean21, acov_mean21, acor_mean21, mean_var21, acov_var21, acor_var21, mean_acov_cor21, mean_acor_cor21, acov_acor_cor21)
  estimate22 <- c(mean_mean22, acov_mean22, acor_mean22, mean_var22, acov_var22, acor_var22, mean_acov_cor22, mean_acor_cor22, acov_acor_cor22)
  estimate31 <- c(mean_mean31, acov_mean31, acor_mean31, mean_var31, acov_var31, acor_var31, mean_acov_cor31, mean_acor_cor31, acov_acor_cor31)
  estimate32 <- c(mean_mean32, acov_mean32, acor_mean32, mean_var32, acov_var32, acor_var32, mean_acov_cor32, mean_acor_cor32, acov_acor_cor32)
  estimate33 <- c(mean_mean33, acov_mean33, acor_mean33, mean_var33, acov_var33, acor_var33, mean_acov_cor33, mean_acor_cor33, acov_acor_cor33)
  estimate34 <- c(mean_mean34, acov_mean34, acor_mean34, mean_var34, acov_var34, acor_var34, mean_acov_cor34, mean_acor_cor34, acov_acor_cor34)
  estimate35 <- c(mean_mean35, acov_mean35, acor_mean35, mean_var35, acov_var35, acor_var35, mean_acov_cor35, mean_acor_cor35, acov_acor_cor35)
  estimate36 <- c(mean_mean36, acov_mean36, acor_mean36, mean_var36, acov_var36, acor_var36, mean_acov_cor36, mean_acor_cor36, acov_acor_cor36)
  estimate37 <- c(mean_mean37, acov_mean37, acor_mean37, mean_var37, acov_var37, acor_var37, mean_acov_cor37, mean_acor_cor37, acov_acor_cor37)
  estimate38 <- c(mean_mean38, acov_mean38, acor_mean38, mean_var38, acov_var38, acor_var38, mean_acov_cor38, mean_acor_cor38, acov_acor_cor38)
  estimate39 <- c(mean_mean39, acov_mean39, acor_mean39, mean_var39, acov_var39, acor_var39, mean_acov_cor39, mean_acor_cor39, acov_acor_cor39)

  # TOJ estimate
  tojestimate <- 3.536 * estimate - 4.072 * (estimate21 + estimate22) / 2 + 1.536 * (estimate31 + estimate32 + estimate33 + estimate34 + estimate35 + estimate36 + estimate37 + estimate38 + estimate39) / 9

  return(tojestimate)

}

#' computing TOJ bias-corrected estimate for T equivalent to 3 modulo 6
#'
#' @param quantity N * (3 * 8) matrix of the means, autocovariances, autocorrelations
#' @param indices indices for bootstrap replications
#'
tojmomentest3 <- function(quantity, indices) {

  # estimated quantities
  mean_est <- quantity[indices, 1]
  mean_est21 <- quantity[indices, 2]
  mean_est22 <- quantity[indices, 3]
  mean_est23 <- quantity[indices, 4]
  mean_est24 <- quantity[indices, 5]
  mean_est31 <- quantity[indices, 6]
  mean_est32 <- quantity[indices, 7]
  mean_est33 <- quantity[indices, 8]

  acov_est <- quantity[indices, 9]
  acov_est21 <- quantity[indices, 10]
  acov_est22 <- quantity[indices, 11]
  acov_est23 <- quantity[indices, 12]
  acov_est24 <- quantity[indices, 13]
  acov_est31 <- quantity[indices, 14]
  acov_est32 <- quantity[indices, 15]
  acov_est33 <- quantity[indices, 16]

  acor_est <- quantity[indices, 17]
  acor_est21 <- quantity[indices, 18]
  acor_est22 <- quantity[indices, 19]
  acor_est23 <- quantity[indices, 20]
  acor_est24 <- quantity[indices, 21]
  acor_est31 <- quantity[indices, 22]
  acor_est32 <- quantity[indices, 23]
  acor_est33 <- quantity[indices, 24]

  # the means
  mean_mean <- mean(mean_est)
  mean_mean21 <- mean(mean_est21)
  mean_mean22 <- mean(mean_est22)
  mean_mean23 <- mean(mean_est23)
  mean_mean24 <- mean(mean_est24)
  mean_mean31 <- mean(mean_est31)
  mean_mean32 <- mean(mean_est32)
  mean_mean33 <- mean(mean_est33)

  acov_mean <- mean(acov_est)
  acov_mean21 <- mean(acov_est21)
  acov_mean22 <- mean(acov_est22)
  acov_mean23 <- mean(acov_est23)
  acov_mean24 <- mean(acov_est24)
  acov_mean31 <- mean(acov_est31)
  acov_mean32 <- mean(acov_est32)
  acov_mean33 <- mean(acov_est33)

  acor_mean <- mean(acor_est)
  acor_mean21 <- mean(acor_est21)
  acor_mean22 <- mean(acor_est22)
  acor_mean23 <- mean(acor_est23)
  acor_mean24 <- mean(acor_est24)
  acor_mean31 <- mean(acor_est31)
  acor_mean32 <- mean(acor_est32)
  acor_mean33 <- mean(acor_est33)

  # the variances
  mean_var <- var(mean_est)
  mean_var21 <- var(mean_est21)
  mean_var22 <- var(mean_est22)
  mean_var23 <- var(mean_est23)
  mean_var24 <- var(mean_est24)
  mean_var31 <- var(mean_est31)
  mean_var32 <- var(mean_est32)
  mean_var33 <- var(mean_est33)

  acov_var <- var(acov_est)
  acov_var21 <- var(acov_est21)
  acov_var22 <- var(acov_est22)
  acov_var23 <- var(acov_est23)
  acov_var24 <- var(acov_est24)
  acov_var31 <- var(acov_est31)
  acov_var32 <- var(acov_est32)
  acov_var33 <- var(acov_est33)

  acor_var <- var(acor_est)
  acor_var21 <- var(acor_est21)
  acor_var22 <- var(acor_est22)
  acor_var23 <- var(acor_est23)
  acor_var24 <- var(acor_est24)
  acor_var31 <- var(acor_est31)
  acor_var32 <- var(acor_est32)
  acor_var33 <- var(acor_est33)

  # the correlations
  mean_acov_cor <- cor(mean_est, acov_est)
  mean_acov_cor21 <- cor(mean_est21, acov_est21)
  mean_acov_cor22 <- cor(mean_est22, acov_est22)
  mean_acov_cor23 <- cor(mean_est23, acov_est23)
  mean_acov_cor24 <- cor(mean_est24, acov_est24)
  mean_acov_cor31 <- cor(mean_est31, acov_est31)
  mean_acov_cor32 <- cor(mean_est32, acov_est32)
  mean_acov_cor33 <- cor(mean_est33, acov_est33)

  mean_acor_cor <- cor(mean_est, acor_est)
  mean_acor_cor21 <- cor(mean_est21, acor_est21)
  mean_acor_cor22 <- cor(mean_est22, acor_est22)
  mean_acor_cor23 <- cor(mean_est23, acor_est23)
  mean_acor_cor24 <- cor(mean_est24, acor_est24)
  mean_acor_cor31 <- cor(mean_est31, acor_est31)
  mean_acor_cor32 <- cor(mean_est32, acor_est32)
  mean_acor_cor33 <- cor(mean_est33, acor_est33)

  acov_acor_cor <- cor(acov_est, acor_est)
  acov_acor_cor21 <- cor(acov_est21, acor_est21)
  acov_acor_cor22 <- cor(acov_est22, acor_est22)
  acov_acor_cor23 <- cor(acov_est23, acor_est23)
  acov_acor_cor24 <- cor(acov_est24, acor_est24)
  acov_acor_cor31 <- cor(acov_est31, acor_est31)
  acov_acor_cor32 <- cor(acov_est32, acor_est32)
  acov_acor_cor33 <- cor(acov_est33, acor_est33)

  # estimates
  estimate <- c(mean_mean, acov_mean, acor_mean, mean_var, acov_var, acor_var, mean_acov_cor, mean_acor_cor, acov_acor_cor)
  estimate21 <- c(mean_mean21, acov_mean21, acor_mean21, mean_var21, acov_var21, acor_var21, mean_acov_cor21, mean_acor_cor21, acov_acor_cor21)
  estimate22 <- c(mean_mean22, acov_mean22, acor_mean22, mean_var22, acov_var22, acor_var22, mean_acov_cor22, mean_acor_cor22, acov_acor_cor22)
  estimate23 <- c(mean_mean23, acov_mean23, acor_mean23, mean_var23, acov_var23, acor_var23, mean_acov_cor23, mean_acor_cor23, acov_acor_cor23)
  estimate24 <- c(mean_mean24, acov_mean24, acor_mean24, mean_var24, acov_var24, acor_var24, mean_acov_cor24, mean_acor_cor24, acov_acor_cor24)
  estimate31 <- c(mean_mean31, acov_mean31, acor_mean31, mean_var31, acov_var31, acor_var31, mean_acov_cor31, mean_acor_cor31, acov_acor_cor31)
  estimate32 <- c(mean_mean32, acov_mean32, acor_mean32, mean_var32, acov_var32, acor_var32, mean_acov_cor32, mean_acor_cor32, acov_acor_cor32)
  estimate33 <- c(mean_mean33, acov_mean33, acor_mean33, mean_var33, acov_var33, acor_var33, mean_acov_cor33, mean_acor_cor33, acov_acor_cor33)

  # TOJ estimate
  tojestimate <- 3.536 * estimate - 4.072 * (estimate21 + estimate22 + estimate23 + estimate24) / 4 + 1.536 * (estimate31 + estimate32 + estimate33) / 3

  return(tojestimate)
}

#' computing TOJ bias-corrected estimate for T equivalent to 4 modulo 6
#'
#' @param quantity N * (3 * 12) matrix of the means, autocovariances, autocorrelations
#' @param indices indices for bootstrap replications
#'
tojmomentest4 <- function(quantity, indices) {

  # estimated quantities
  mean_est <- quantity[indices, 1]
  mean_est21 <- quantity[indices, 2]
  mean_est22 <- quantity[indices, 3]
  mean_est31 <- quantity[indices, 4]
  mean_est32 <- quantity[indices, 5]
  mean_est33 <- quantity[indices, 6]
  mean_est34 <- quantity[indices, 7]
  mean_est35 <- quantity[indices, 8]
  mean_est36 <- quantity[indices, 9]
  mean_est37 <- quantity[indices, 10]
  mean_est38 <- quantity[indices, 11]
  mean_est39 <- quantity[indices, 12]

  acov_est <- quantity[indices, 13]
  acov_est21 <- quantity[indices, 14]
  acov_est22 <- quantity[indices, 15]
  acov_est31 <- quantity[indices, 16]
  acov_est32 <- quantity[indices, 17]
  acov_est33 <- quantity[indices, 18]
  acov_est34 <- quantity[indices, 19]
  acov_est35 <- quantity[indices, 20]
  acov_est36 <- quantity[indices, 21]
  acov_est37 <- quantity[indices, 22]
  acov_est38 <- quantity[indices, 23]
  acov_est39 <- quantity[indices, 24]

  acor_est <- quantity[indices, 25]
  acor_est21 <- quantity[indices, 26]
  acor_est22 <- quantity[indices, 27]
  acor_est31 <- quantity[indices, 28]
  acor_est32 <- quantity[indices, 29]
  acor_est33 <- quantity[indices, 30]
  acor_est34 <- quantity[indices, 31]
  acor_est35 <- quantity[indices, 32]
  acor_est36 <- quantity[indices, 33]
  acor_est37 <- quantity[indices, 34]
  acor_est38 <- quantity[indices, 35]
  acor_est39 <- quantity[indices, 36]

  # the means
  mean_mean <- mean(mean_est)
  mean_mean21 <- mean(mean_est21)
  mean_mean22 <- mean(mean_est22)
  mean_mean31 <- mean(mean_est31)
  mean_mean32 <- mean(mean_est32)
  mean_mean33 <- mean(mean_est33)
  mean_mean34 <- mean(mean_est34)
  mean_mean35 <- mean(mean_est35)
  mean_mean36 <- mean(mean_est36)
  mean_mean37 <- mean(mean_est37)
  mean_mean38 <- mean(mean_est38)
  mean_mean39 <- mean(mean_est39)

  acov_mean <- mean(acov_est)
  acov_mean21 <- mean(acov_est21)
  acov_mean22 <- mean(acov_est22)
  acov_mean31 <- mean(acov_est31)
  acov_mean32 <- mean(acov_est32)
  acov_mean33 <- mean(acov_est33)
  acov_mean34 <- mean(acov_est34)
  acov_mean35 <- mean(acov_est35)
  acov_mean36 <- mean(acov_est36)
  acov_mean37 <- mean(acov_est37)
  acov_mean38 <- mean(acov_est38)
  acov_mean39 <- mean(acov_est39)

  acor_mean <- mean(acor_est)
  acor_mean21 <- mean(acor_est21)
  acor_mean22 <- mean(acor_est22)
  acor_mean31 <- mean(acor_est31)
  acor_mean32 <- mean(acor_est32)
  acor_mean33 <- mean(acor_est33)
  acor_mean34 <- mean(acor_est34)
  acor_mean35 <- mean(acor_est35)
  acor_mean36 <- mean(acor_est36)
  acor_mean37 <- mean(acor_est37)
  acor_mean38 <- mean(acor_est38)
  acor_mean39 <- mean(acor_est39)

  # the variances
  mean_var <- var(mean_est)
  mean_var21 <- var(mean_est21)
  mean_var22 <- var(mean_est22)
  mean_var31 <- var(mean_est31)
  mean_var32 <- var(mean_est32)
  mean_var33 <- var(mean_est33)
  mean_var34 <- var(mean_est34)
  mean_var35 <- var(mean_est35)
  mean_var36 <- var(mean_est36)
  mean_var37 <- var(mean_est37)
  mean_var38 <- var(mean_est38)
  mean_var39 <- var(mean_est39)

  acov_var <- var(acov_est)
  acov_var21 <- var(acov_est21)
  acov_var22 <- var(acov_est22)
  acov_var31 <- var(acov_est31)
  acov_var32 <- var(acov_est32)
  acov_var33 <- var(acov_est33)
  acov_var34 <- var(acov_est34)
  acov_var35 <- var(acov_est35)
  acov_var36 <- var(acov_est36)
  acov_var37 <- var(acov_est37)
  acov_var38 <- var(acov_est38)
  acov_var39 <- var(acov_est39)

  acor_var <- var(acor_est)
  acor_var21 <- var(acor_est21)
  acor_var22 <- var(acor_est22)
  acor_var31 <- var(acor_est31)
  acor_var32 <- var(acor_est32)
  acor_var33 <- var(acor_est33)
  acor_var34 <- var(acor_est34)
  acor_var35 <- var(acor_est35)
  acor_var36 <- var(acor_est36)
  acor_var37 <- var(acor_est37)
  acor_var38 <- var(acor_est38)
  acor_var39 <- var(acor_est39)

  # the correlations
  mean_acov_cor <- cor(mean_est, acov_est)
  mean_acov_cor21 <- cor(mean_est21, acov_est21)
  mean_acov_cor22 <- cor(mean_est22, acov_est22)
  mean_acov_cor31 <- cor(mean_est31, acov_est31)
  mean_acov_cor32 <- cor(mean_est32, acov_est32)
  mean_acov_cor33 <- cor(mean_est33, acov_est33)
  mean_acov_cor34 <- cor(mean_est34, acov_est34)
  mean_acov_cor35 <- cor(mean_est35, acov_est35)
  mean_acov_cor36 <- cor(mean_est36, acov_est36)
  mean_acov_cor37 <- cor(mean_est37, acov_est37)
  mean_acov_cor38 <- cor(mean_est38, acov_est38)
  mean_acov_cor39 <- cor(mean_est39, acov_est39)

  mean_acor_cor <- cor(mean_est, acor_est)
  mean_acor_cor21 <- cor(mean_est21, acor_est21)
  mean_acor_cor22 <- cor(mean_est22, acor_est22)
  mean_acor_cor31 <- cor(mean_est31, acor_est31)
  mean_acor_cor32 <- cor(mean_est32, acor_est32)
  mean_acor_cor33 <- cor(mean_est33, acor_est33)
  mean_acor_cor34 <- cor(mean_est34, acor_est34)
  mean_acor_cor35 <- cor(mean_est35, acor_est35)
  mean_acor_cor36 <- cor(mean_est36, acor_est36)
  mean_acor_cor37 <- cor(mean_est37, acor_est37)
  mean_acor_cor38 <- cor(mean_est38, acor_est38)
  mean_acor_cor39 <- cor(mean_est39, acor_est39)

  acov_acor_cor <- cor(acov_est, acor_est)
  acov_acor_cor21 <- cor(acov_est21, acor_est21)
  acov_acor_cor22 <- cor(acov_est22, acor_est22)
  acov_acor_cor31 <- cor(acov_est31, acor_est31)
  acov_acor_cor32 <- cor(acov_est32, acor_est32)
  acov_acor_cor33 <- cor(acov_est33, acor_est33)
  acov_acor_cor34 <- cor(acov_est34, acor_est34)
  acov_acor_cor35 <- cor(acov_est35, acor_est35)
  acov_acor_cor36 <- cor(acov_est36, acor_est36)
  acov_acor_cor37 <- cor(acov_est37, acor_est37)
  acov_acor_cor38 <- cor(acov_est38, acor_est38)
  acov_acor_cor39 <- cor(acov_est39, acor_est39)

  # estimates
  estimate <- c(mean_mean, acov_mean, acor_mean, mean_var, acov_var, acor_var, mean_acov_cor, mean_acor_cor, acov_acor_cor)
  estimate21 <- c(mean_mean21, acov_mean21, acor_mean21, mean_var21, acov_var21, acor_var21, mean_acov_cor21, mean_acor_cor21, acov_acor_cor21)
  estimate22 <- c(mean_mean22, acov_mean22, acor_mean22, mean_var22, acov_var22, acor_var22, mean_acov_cor22, mean_acor_cor22, acov_acor_cor22)
  estimate31 <- c(mean_mean31, acov_mean31, acor_mean31, mean_var31, acov_var31, acor_var31, mean_acov_cor31, mean_acor_cor31, acov_acor_cor31)
  estimate32 <- c(mean_mean32, acov_mean32, acor_mean32, mean_var32, acov_var32, acor_var32, mean_acov_cor32, mean_acor_cor32, acov_acor_cor32)
  estimate33 <- c(mean_mean33, acov_mean33, acor_mean33, mean_var33, acov_var33, acor_var33, mean_acov_cor33, mean_acor_cor33, acov_acor_cor33)
  estimate34 <- c(mean_mean34, acov_mean34, acor_mean34, mean_var34, acov_var34, acor_var34, mean_acov_cor34, mean_acor_cor34, acov_acor_cor34)
  estimate35 <- c(mean_mean35, acov_mean35, acor_mean35, mean_var35, acov_var35, acor_var35, mean_acov_cor35, mean_acor_cor35, acov_acor_cor35)
  estimate36 <- c(mean_mean36, acov_mean36, acor_mean36, mean_var36, acov_var36, acor_var36, mean_acov_cor36, mean_acor_cor36, acov_acor_cor36)
  estimate37 <- c(mean_mean37, acov_mean37, acor_mean37, mean_var37, acov_var37, acor_var37, mean_acov_cor37, mean_acor_cor37, acov_acor_cor37)
  estimate38 <- c(mean_mean38, acov_mean38, acor_mean38, mean_var38, acov_var38, acor_var38, mean_acov_cor38, mean_acor_cor38, acov_acor_cor38)
  estimate39 <- c(mean_mean39, acov_mean39, acor_mean39, mean_var39, acov_var39, acor_var39, mean_acov_cor39, mean_acor_cor39, acov_acor_cor39)

  # TOJ estimate
  tojestimate <- 3.536 * estimate - 4.072 * (estimate21 + estimate22) / 2 + 1.536 * (estimate31 + estimate32 + estimate33 + estimate34 + estimate35 + estimate36 + estimate37 + estimate38 + estimate39) / 9

  return(tojestimate)

}

#' computing TOJ bias-corrected estimate for T equivalent to 5 modulo 6
#'
#' @param quantity N * (3 * 14) matrix of the means, autocovariances, autocorrelations
#' @param indices indices for bootstrap replications
#'
tojmomentest5 <- function(quantity, indices) {

  # estimated quantities
  mean_est <- quantity[indices, 1]
  mean_est21 <- quantity[indices, 2]
  mean_est22 <- quantity[indices, 3]
  mean_est23 <- quantity[indices, 4]
  mean_est24 <- quantity[indices, 5]
  mean_est31 <- quantity[indices, 6]
  mean_est32 <- quantity[indices, 7]
  mean_est33 <- quantity[indices, 8]
  mean_est34 <- quantity[indices, 9]
  mean_est35 <- quantity[indices, 10]
  mean_est36 <- quantity[indices, 11]
  mean_est37 <- quantity[indices, 12]
  mean_est38 <- quantity[indices, 13]
  mean_est39 <- quantity[indices, 14]

  acov_est <- quantity[indices, 15]
  acov_est21 <- quantity[indices, 16]
  acov_est22 <- quantity[indices, 17]
  acov_est23 <- quantity[indices, 18]
  acov_est24 <- quantity[indices, 19]
  acov_est31 <- quantity[indices, 20]
  acov_est32 <- quantity[indices, 21]
  acov_est33 <- quantity[indices, 22]
  acov_est34 <- quantity[indices, 23]
  acov_est35 <- quantity[indices, 24]
  acov_est36 <- quantity[indices, 25]
  acov_est37 <- quantity[indices, 26]
  acov_est38 <- quantity[indices, 27]
  acov_est39 <- quantity[indices, 28]

  acor_est <- quantity[indices, 29]
  acor_est21 <- quantity[indices, 30]
  acor_est22 <- quantity[indices, 31]
  acor_est23 <- quantity[indices, 32]
  acor_est24 <- quantity[indices, 33]
  acor_est31 <- quantity[indices, 34]
  acor_est32 <- quantity[indices, 35]
  acor_est33 <- quantity[indices, 36]
  acor_est34 <- quantity[indices, 37]
  acor_est35 <- quantity[indices, 38]
  acor_est36 <- quantity[indices, 39]
  acor_est37 <- quantity[indices, 40]
  acor_est38 <- quantity[indices, 41]
  acor_est39 <- quantity[indices, 42]

  # the means
  mean_mean <- mean(mean_est)
  mean_mean21 <- mean(mean_est21)
  mean_mean22 <- mean(mean_est22)
  mean_mean23 <- mean(mean_est23)
  mean_mean24 <- mean(mean_est24)
  mean_mean31 <- mean(mean_est31)
  mean_mean32 <- mean(mean_est32)
  mean_mean33 <- mean(mean_est33)
  mean_mean34 <- mean(mean_est34)
  mean_mean35 <- mean(mean_est35)
  mean_mean36 <- mean(mean_est36)
  mean_mean37 <- mean(mean_est37)
  mean_mean38 <- mean(mean_est38)
  mean_mean39 <- mean(mean_est39)

  acov_mean <- mean(acov_est)
  acov_mean21 <- mean(acov_est21)
  acov_mean22 <- mean(acov_est22)
  acov_mean23 <- mean(acov_est23)
  acov_mean24 <- mean(acov_est24)
  acov_mean31 <- mean(acov_est31)
  acov_mean32 <- mean(acov_est32)
  acov_mean33 <- mean(acov_est33)
  acov_mean34 <- mean(acov_est34)
  acov_mean35 <- mean(acov_est35)
  acov_mean36 <- mean(acov_est36)
  acov_mean37 <- mean(acov_est37)
  acov_mean38 <- mean(acov_est38)
  acov_mean39 <- mean(acov_est39)

  acor_mean <- mean(acor_est)
  acor_mean21 <- mean(acor_est21)
  acor_mean22 <- mean(acor_est22)
  acor_mean23 <- mean(acor_est23)
  acor_mean24 <- mean(acor_est24)
  acor_mean31 <- mean(acor_est31)
  acor_mean32 <- mean(acor_est32)
  acor_mean33 <- mean(acor_est33)
  acor_mean34 <- mean(acor_est34)
  acor_mean35 <- mean(acor_est35)
  acor_mean36 <- mean(acor_est36)
  acor_mean37 <- mean(acor_est37)
  acor_mean38 <- mean(acor_est38)
  acor_mean39 <- mean(acor_est39)

  # the variances
  mean_var <- var(mean_est)
  mean_var21 <- var(mean_est21)
  mean_var22 <- var(mean_est22)
  mean_var23 <- var(mean_est23)
  mean_var24 <- var(mean_est24)
  mean_var31 <- var(mean_est31)
  mean_var32 <- var(mean_est32)
  mean_var33 <- var(mean_est33)
  mean_var34 <- var(mean_est34)
  mean_var35 <- var(mean_est35)
  mean_var36 <- var(mean_est36)
  mean_var37 <- var(mean_est37)
  mean_var38 <- var(mean_est38)
  mean_var39 <- var(mean_est39)

  acov_var <- var(acov_est)
  acov_var21 <- var(acov_est21)
  acov_var22 <- var(acov_est22)
  acov_var23 <- var(acov_est23)
  acov_var24 <- var(acov_est24)
  acov_var31 <- var(acov_est31)
  acov_var32 <- var(acov_est32)
  acov_var33 <- var(acov_est33)
  acov_var34 <- var(acov_est34)
  acov_var35 <- var(acov_est35)
  acov_var36 <- var(acov_est36)
  acov_var37 <- var(acov_est37)
  acov_var38 <- var(acov_est38)
  acov_var39 <- var(acov_est39)

  acor_var <- var(acor_est)
  acor_var21 <- var(acor_est21)
  acor_var22 <- var(acor_est22)
  acor_var23 <- var(acor_est23)
  acor_var24 <- var(acor_est24)
  acor_var31 <- var(acor_est31)
  acor_var32 <- var(acor_est32)
  acor_var33 <- var(acor_est33)
  acor_var34 <- var(acor_est34)
  acor_var35 <- var(acor_est35)
  acor_var36 <- var(acor_est36)
  acor_var37 <- var(acor_est37)
  acor_var38 <- var(acor_est38)
  acor_var39 <- var(acor_est39)

  # the correlations
  mean_acov_cor <- cor(mean_est, acov_est)
  mean_acov_cor21 <- cor(mean_est21, acov_est21)
  mean_acov_cor22 <- cor(mean_est22, acov_est22)
  mean_acov_cor23 <- cor(mean_est23, acov_est23)
  mean_acov_cor24 <- cor(mean_est24, acov_est24)
  mean_acov_cor31 <- cor(mean_est31, acov_est31)
  mean_acov_cor32 <- cor(mean_est32, acov_est32)
  mean_acov_cor33 <- cor(mean_est33, acov_est33)
  mean_acov_cor34 <- cor(mean_est34, acov_est34)
  mean_acov_cor35 <- cor(mean_est35, acov_est35)
  mean_acov_cor36 <- cor(mean_est36, acov_est36)
  mean_acov_cor37 <- cor(mean_est37, acov_est37)
  mean_acov_cor38 <- cor(mean_est38, acov_est38)
  mean_acov_cor39 <- cor(mean_est39, acov_est39)

  mean_acor_cor <- cor(mean_est, acor_est)
  mean_acor_cor21 <- cor(mean_est21, acor_est21)
  mean_acor_cor22 <- cor(mean_est22, acor_est22)
  mean_acor_cor23 <- cor(mean_est23, acor_est23)
  mean_acor_cor24 <- cor(mean_est24, acor_est24)
  mean_acor_cor31 <- cor(mean_est31, acor_est31)
  mean_acor_cor32 <- cor(mean_est32, acor_est32)
  mean_acor_cor33 <- cor(mean_est33, acor_est33)
  mean_acor_cor34 <- cor(mean_est34, acor_est34)
  mean_acor_cor35 <- cor(mean_est35, acor_est35)
  mean_acor_cor36 <- cor(mean_est36, acor_est36)
  mean_acor_cor37 <- cor(mean_est37, acor_est37)
  mean_acor_cor38 <- cor(mean_est38, acor_est38)
  mean_acor_cor39 <- cor(mean_est39, acor_est39)

  acov_acor_cor <- cor(acov_est, acor_est)
  acov_acor_cor21 <- cor(acov_est21, acor_est21)
  acov_acor_cor22 <- cor(acov_est22, acor_est22)
  acov_acor_cor23 <- cor(acov_est23, acor_est23)
  acov_acor_cor24 <- cor(acov_est24, acor_est24)
  acov_acor_cor31 <- cor(acov_est31, acor_est31)
  acov_acor_cor32 <- cor(acov_est32, acor_est32)
  acov_acor_cor33 <- cor(acov_est33, acor_est33)
  acov_acor_cor34 <- cor(acov_est34, acor_est34)
  acov_acor_cor35 <- cor(acov_est35, acor_est35)
  acov_acor_cor36 <- cor(acov_est36, acor_est36)
  acov_acor_cor37 <- cor(acov_est37, acor_est37)
  acov_acor_cor38 <- cor(acov_est38, acor_est38)
  acov_acor_cor39 <- cor(acov_est39, acor_est39)

  # estimates
  estimate <- c(mean_mean, acov_mean, acor_mean, mean_var, acov_var, acor_var, mean_acov_cor, mean_acor_cor, acov_acor_cor)
  estimate21 <- c(mean_mean21, acov_mean21, acor_mean21, mean_var21, acov_var21, acor_var21, mean_acov_cor21, mean_acor_cor21, acov_acor_cor21)
  estimate22 <- c(mean_mean22, acov_mean22, acor_mean22, mean_var22, acov_var22, acor_var22, mean_acov_cor22, mean_acor_cor22, acov_acor_cor22)
  estimate23 <- c(mean_mean23, acov_mean23, acor_mean23, mean_var23, acov_var23, acor_var23, mean_acov_cor23, mean_acor_cor23, acov_acor_cor23)
  estimate24 <- c(mean_mean24, acov_mean24, acor_mean24, mean_var24, acov_var24, acor_var24, mean_acov_cor24, mean_acor_cor24, acov_acor_cor24)
  estimate31 <- c(mean_mean31, acov_mean31, acor_mean31, mean_var31, acov_var31, acor_var31, mean_acov_cor31, mean_acor_cor31, acov_acor_cor31)
  estimate32 <- c(mean_mean32, acov_mean32, acor_mean32, mean_var32, acov_var32, acor_var32, mean_acov_cor32, mean_acor_cor32, acov_acor_cor32)
  estimate33 <- c(mean_mean33, acov_mean33, acor_mean33, mean_var33, acov_var33, acor_var33, mean_acov_cor33, mean_acor_cor33, acov_acor_cor33)
  estimate34 <- c(mean_mean34, acov_mean34, acor_mean34, mean_var34, acov_var34, acor_var34, mean_acov_cor34, mean_acor_cor34, acov_acor_cor34)
  estimate35 <- c(mean_mean35, acov_mean35, acor_mean35, mean_var35, acov_var35, acor_var35, mean_acov_cor35, mean_acor_cor35, acov_acor_cor35)
  estimate36 <- c(mean_mean36, acov_mean36, acor_mean36, mean_var36, acov_var36, acor_var36, mean_acov_cor36, mean_acor_cor36, acov_acor_cor36)
  estimate37 <- c(mean_mean37, acov_mean37, acor_mean37, mean_var37, acov_var37, acor_var37, mean_acov_cor37, mean_acor_cor37, acov_acor_cor37)
  estimate38 <- c(mean_mean38, acov_mean38, acor_mean38, mean_var38, acov_var38, acor_var38, mean_acov_cor38, mean_acor_cor38, acov_acor_cor38)
  estimate39 <- c(mean_mean39, acov_mean39, acor_mean39, mean_var39, acov_var39, acor_var39, mean_acov_cor39, mean_acor_cor39, acov_acor_cor39)

  # TOJ estimate
  tojestimate <- 3.536 * estimate - 4.072 * (estimate21 + estimate22 + estimate23 + estimate24) / 4 + 1.536 * (estimate31 + estimate32 + estimate33 + estimate34 + estimate35 + estimate36 + estimate37 + estimate38 + estimate39) / 9

  return(tojestimate)

}