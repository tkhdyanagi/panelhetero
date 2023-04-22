#' The TOJ bias-corrected empirical CDF estimation
#'
#' The `tojecdf()` function enables to implement the TOJ bias-corrected
#' estimation of the cumulative distribution function (CDF) of
#' the heterogeneous mean, the heterogeneous autocovariance, and
#' the heterogeneous autocorrelation.
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
#' panelhetero::tojecdf(data = data, R = 50)
#'
#' @references Okui, R. and Yanagi, T., 2019.
#' Panel data analysis with heterogeneous dynamics.
#' Journal of Econometrics, 212(2), pp.451-475.
#'
#' @export
#'
tojecdf <- function(data,
                    acov_order = 0,
                    acor_order = 1,
                    R = 1000,
                    ci = TRUE) {

  # Error handling -------------------------------------------------------------

  error1(data = data,
         acov_order = acov_order,
         acor_order = acor_order,
         R = R)

  # Variable definitions -------------------------------------------------------

  x <- y <- NULL

  # Omit NA
  data <- stats::na.omit(data)

  # Sample size
  N <- nrow(data)
  S <- ncol(data)

  # Estimated means, autocovariances, autocorrelations
  mean_est <- rowMeans(data)
  acov_est <- apply(data, MARGIN = 1, acov, acov_order = acov_order)
  acor_est <- apply(data, MARGIN = 1, acor, acor_order = acor_order)

  # TOJ bias-corrected estimation ----------------------------------------------

  if (S %% 6 == 0) {

    # Split panel data for T equivalent to 0 modulo 6
    data21 <- data[, 1:(S / 2)]
    data22 <- data[, (S / 2 + 1):S]
    data31 <- data[, 1:(S / 3)]
    data32 <- data[, (S / 3 + 1):(2*S / 3)]
    data33 <- data[, (2 * S / 3 + 1):S]

    # Estimated quantities for split panel data
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

    # Function for bootstrap confidence interval
    mean_ci_func <- Vectorize(FUN = function(x) {

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(mean_est,
                                           mean_est21,
                                           mean_est22,
                                           mean_est31,
                                           mean_est32,
                                           mean_est33) - x,
                              statistic = toj0_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acov_est,
                                           acov_est21,
                                           acov_est22,
                                           acov_est31,
                                           acov_est32,
                                           acov_est33) - x,
                              statistic = toj0_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acor_est,
                                           acor_est21,
                                           acor_est22,
                                           acor_est31,
                                           acor_est32,
                                           acor_est33) - x,
                              statistic = toj0_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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
      mean_x <- seq(min(mean_est),
                    max(mean_est),
                    length = 101)

      mean_y <- tojecdfest0(x   = mean_x,
                            X   = mean_est,
                            X21 = mean_est21,
                            X22 = mean_est22,
                            X31 = mean_est31,
                            X32 = mean_est32,
                            X33 = mean_est33)

      mean_toj <- Rearrangement::rearrangement(x = data.frame(x = mean_x),
                                               y = mean_y)

      mean_plot <- ggplot2::ggplot(
        data.frame(x = mean_x, y = mean_toj),
        ggplot2::aes(x = x, y = y)
      ) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(mean_est), max(mean_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::theme_bw()

      # Autocovariance
      acov_x <- seq(min(acov_est),
                    max(acov_est),
                    length = 101)

      acov_y <- tojecdfest0(x   = acov_x,
                            X   = acov_est,
                            X21 = acov_est21,
                            X22 = acov_est22,
                            X31 = acov_est31,
                            X32 = acov_est32,
                            X33 = acov_est33)

      acov_toj <- Rearrangement::rearrangement(x = data.frame(x = acov_x),
                                               y = acov_y)

      acov_plot <- ggplot2::ggplot(
        data.frame(x = acov_x, y = acov_toj),
        ggplot2::aes(x = x, y = y)
      ) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acov_est), max(acov_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_x <- seq(min(acor_est),
                    max(acor_est),
                    length = 101)

      acor_y <- tojecdfest0(x   = acor_x,
                            X   = acor_est,
                            X21 = acor_est21,
                            X22 = acor_est22,
                            X31 = acor_est31,
                            X32 = acor_est32,
                            X33 = acor_est33)

      acor_toj <- Rearrangement::rearrangement(x = data.frame(x = acor_x),
                                               y = acor_y)

      acor_plot <- ggplot2::ggplot(
        data.frame(x = acor_x, y = acor_toj),
        ggplot2::aes(x = x, y = y)
      ) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acor_est), max(acor_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::theme_bw()

    }

    if (ci) {

      # Mean
      mean_toj <- cbind(tojecdfest0(x   = mean_grid,
                                    X   = mean_est,
                                    X21 = mean_est21,
                                    X22 = mean_est22,
                                    X31 = mean_est31,
                                    X32 = mean_est32,
                                    X33 = mean_est33),
                        t(mean_ci))

      mean_plot <- ggplot2::ggplot(
        data = data.frame(x = mean_grid),
        ggplot2::aes(x = mean_grid)
      ) +
        ggplot2::geom_line(
          ggplot2::aes(
            x = mean_grid,
            y = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_toj[, 1]
            )
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = mean_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_toj[, 2]
            ),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_toj[, 3]
            )
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocovariance
      acov_toj <- cbind(tojecdfest0(x   = acov_grid,
                                    X   = acov_est,
                                    X21 = acov_est21,
                                    X22 = acov_est22,
                                    X31 = acov_est31,
                                    X32 = acov_est32,
                                    X33 = acov_est33),
                        t(acov_ci))

      acov_plot <- ggplot2::ggplot(
        data = data.frame(x = acov_grid),
        ggplot2::aes(x = acov_grid)
      ) +
        ggplot2::geom_line(
          ggplot2::aes(
            x = acov_grid,
            y = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_toj[, 1]
            )
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = acov_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_toj[, 2]
            ),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_toj[, 3]
            )
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_toj <- cbind(tojecdfest0(x   = acor_grid,
                                    X   = acor_est,
                                    X21 = acor_est21,
                                    X22 = acor_est22,
                                    X31 = acor_est31,
                                    X32 = acor_est32,
                                    X33 = acor_est33),
                        t(acor_ci))

      acor_plot <- ggplot2::ggplot(
        data = data.frame(x = acor_grid),
        ggplot2::aes(x = acor_grid)
      ) +
        ggplot2::geom_line(
          ggplot2::aes(
            x = acor_grid,
            y = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_toj[, 1]
            )
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = acor_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_toj[, 2]
            ),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_toj[, 3]
            )
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

    }

    # Functions without rearrangement
    mean_func <- function(x) {
      tojecdfest0(x   = x,
                  X   = mean_est,
                  X21 = mean_est21,
                  X22 = mean_est22,
                  X31 = mean_est31,
                  X32 = mean_est32,
                  X33 = mean_est33)
    }

    acov_func <- function(x) {
      tojecdfest0(x   = x,
                  X   = acov_est,
                  X21 = acov_est21,
                  X22 = acov_est22,
                  X31 = acov_est31,
                  X32 = acov_est32,
                  X33 = acov_est33)
    }

    acor_func <- function(x) {
      tojecdfest0(x   = x,
                  X   = acor_est,
                  X21 = acor_est21,
                  X22 = acor_est22,
                  X31 = acor_est31,
                  X32 = acor_est32,
                  X33 = acor_est33)
    }

  } else if (S %% 6 == 1) {

    # Split  panel data for T equivalent to 1 modulo 6
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

    # Function for bootstrap confidence interval
    mean_ci_func <- Vectorize(FUN = function(x) {

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(mean_est,
                                           mean_est21,
                                           mean_est22,
                                           mean_est23,
                                           mean_est24,
                                           mean_est31,
                                           mean_est32,
                                           mean_est33,
                                           mean_est34,
                                           mean_est35,
                                           mean_est36,
                                           mean_est37,
                                           mean_est38,
                                           mean_est39) - x,
                              statistic = toj1_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acov_est,
                                           acov_est21,
                                           acov_est22,
                                           acov_est23,
                                           acov_est24,
                                           acov_est31,
                                           acov_est32,
                                           acov_est33,
                                           acov_est34,
                                           acov_est35,
                                           acov_est36,
                                           acov_est37,
                                           acov_est38,
                                           acov_est39) - x,
                              statistic = toj1_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acor_est,
                                           acor_est21,
                                           acor_est22,
                                           acor_est23,
                                           acor_est24,
                                           acor_est31,
                                           acor_est32,
                                           acor_est33,
                                           acor_est34,
                                           acor_est35,
                                           acor_est36,
                                           acor_est37,
                                           acor_est38,
                                           acor_est39) - x,
                              statistic = toj1_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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
      mean_x <- seq(min(mean_est),
                    max(mean_est),
                    length = 1000)

      mean_y <- tojecdfest1(x   = mean_x,
                            X   = mean_est,
                            X21 = mean_est21,
                            X22 = mean_est22,
                            X23 = mean_est23,
                            X24 = mean_est24,
                            X31 = mean_est31,
                            X32 = mean_est32,
                            X33 = mean_est33,
                            X34 = mean_est34,
                            X35 = mean_est35,
                            X36 = mean_est36,
                            X37 = mean_est37,
                            X38 = mean_est38,
                            X39 = mean_est39)

      mean_toj <- Rearrangement::rearrangement(x = data.frame(x = mean_x),
                                               y = mean_y)

      mean_plot <- ggplot2::ggplot(data.frame(x = mean_x, y = mean_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(mean_est), max(mean_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::theme_bw()

      # Autocovariance
      acov_x <- seq(min(acov_est),
                    max(acov_est),
                    length = 1000)

      acov_y <- tojecdfest1(x   = acov_x,
                            X   = acov_est,
                            X21 = acov_est21,
                            X22 = acov_est22,
                            X23 = acov_est23,
                            X24 = acov_est24,
                            X31 = acov_est31,
                            X32 = acov_est32,
                            X33 = acov_est33,
                            X34 = acov_est34,
                            X35 = acov_est35,
                            X36 = acov_est36,
                            X37 = acov_est37,
                            X38 = acov_est38,
                            X39 = acov_est39)

      acov_toj <- Rearrangement::rearrangement(x = data.frame(x = acov_x),
                                               y = acov_y)

      acov_plot <- ggplot2::ggplot(data.frame(x = acov_x, y = acov_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acov_est), max(acov_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_x <- seq(min(acor_est),
                    max(acor_est),
                    length = 1000)

      acor_y <- tojecdfest1(x   = acor_x,
                            X   = acor_est,
                            X21 = acor_est21,
                            X22 = acor_est22,
                            X23 = acor_est23,
                            X24 = acor_est24,
                            X31 = acor_est31,
                            X32 = acor_est32,
                            X33 = acor_est33,
                            X34 = acor_est34,
                            X35 = acor_est35,
                            X36 = acor_est36,
                            X37 = acor_est37,
                            X38 = acor_est38,
                            X39 = acor_est39)

      acor_toj <- Rearrangement::rearrangement(x = data.frame(x = acor_x),
                                               y = acor_y)

      acor_plot <- ggplot2::ggplot(data.frame(x = acor_x, y = acor_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acor_est), max(acor_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::theme_bw()

    }

    if (ci) {

      # Mean
      mean_toj <- cbind(tojecdfest1(x   = mean_grid,
                                    X   = mean_est,
                                    X21 = mean_est21,
                                    X22 = mean_est22,
                                    X23 = mean_est23,
                                    X24 = mean_est24,
                                    X31 = mean_est31,
                                    X32 = mean_est32,
                                    X33 = mean_est33,
                                    X34 = mean_est34,
                                    X35 = mean_est35,
                                    X36 = mean_est36,
                                    X37 = mean_est37,
                                    X38 = mean_est38,
                                    X39 = mean_est39),
                        t(mean_ci))

      mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_grid),
                                   ggplot2::aes(x = mean_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(x = mean_grid,
                       y = Rearrangement::rearrangement(
                         x = data.frame(x = mean_grid),
                         y = mean_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(x = mean_grid,
                       ymin = Rearrangement::rearrangement(
                         x = data.frame(x = mean_grid),
                         y = mean_toj[, 2]),
                       ymax = Rearrangement::rearrangement(
                         x = data.frame(x = mean_grid),
                         y = mean_toj[, 3])
          ), alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocovariance
      acov_toj <- cbind(tojecdfest1(x   = acov_grid,
                                    X   = acov_est,
                                    X21 = acov_est21,
                                    X22 = acov_est22,
                                    X23 = acov_est23,
                                    X24 = acov_est24,
                                    X31 = acov_est31,
                                    X32 = acov_est32,
                                    X33 = acov_est33,
                                    X34 = acov_est34,
                                    X35 = acov_est35,
                                    X36 = acov_est36,
                                    X37 = acov_est37,
                                    X38 = acov_est38,
                                    X39 = acov_est39),
                        t(acov_ci))

      acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_grid),
                                   ggplot2::aes(x = acov_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(x = acov_grid,
                       y = Rearrangement::rearrangement(
                         x = data.frame(x = acov_grid),
                         y = acov_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(x = acov_grid,
                       ymin = Rearrangement::rearrangement(
                         x = data.frame(x = acov_grid),
                         y = acov_toj[, 2]),
                       ymax = Rearrangement::rearrangement(
                         x = data.frame(x = acov_grid),
                         y = acov_toj[, 3])
          ), alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_toj <- cbind(tojecdfest1(x   = acor_grid,
                                    X   = acor_est,
                                    X21 = acor_est21,
                                    X22 = acor_est22,
                                    X23 = acor_est23,
                                    X24 = acor_est24,
                                    X31 = acor_est31,
                                    X32 = acor_est32,
                                    X33 = acor_est33,
                                    X34 = acor_est34,
                                    X35 = acor_est35,
                                    X36 = acor_est36,
                                    X37 = acor_est37,
                                    X38 = acor_est38,
                                    X39 = acor_est39),
                        t(acor_ci))

      acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_grid),
                                   ggplot2::aes(x = acor_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(x = acor_grid,
                       y = Rearrangement::rearrangement(
                         x = data.frame(x = acor_grid),
                         y = acor_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(x = acor_grid,
                       ymin = Rearrangement::rearrangement(
                         x = data.frame(x = acor_grid),
                         y = acor_toj[, 2]),
                       ymax = Rearrangement::rearrangement(
                         x = data.frame(x = acor_grid),
                         y = acor_toj[, 3])
          ), alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

    }

    # Functions without rearrangement
    mean_func <- function(x) {
      tojecdfest1(x   = x,
                  X   = mean_est,
                  X21 = mean_est21,
                  X22 = mean_est22,
                  X23 = mean_est23,
                  X24 = mean_est24,
                  X31 = mean_est31,
                  X32 = mean_est32,
                  X33 = mean_est33,
                  X34 = mean_est34,
                  X35 = mean_est35,
                  X36 = mean_est36,
                  X37 = mean_est37,
                  X38 = mean_est38,
                  X39 = mean_est39)
    }

    acov_func <- function(x) {
      tojecdfest1(x = x,
                  X = acov_est,
                  X21 = acov_est21,
                  X22 = acov_est22,
                  X23 = acov_est23,
                  X24 = acov_est24,
                  X31 = acov_est31,
                  X32 = acov_est32,
                  X33 = acov_est33,
                  X34 = acov_est34,
                  X35 = acov_est35,
                  X36 = acov_est36,
                  X37 = acov_est37,
                  X38 = acov_est38,
                  X39 = acov_est39)
    }

    acor_func <- function(x) {
      tojecdfest1(x = x,
                  X = acor_est,
                  X21 = acor_est21,
                  X22 = acor_est22,
                  X23 = acor_est23,
                  X24 = acor_est24,
                  X31 = acor_est31,
                  X32 = acor_est32,
                  X33 = acor_est33,
                  X34 = acor_est34,
                  X35 = acor_est35,
                  X36 = acor_est36,
                  X37 = acor_est37,
                  X38 = acor_est38,
                  X39 = acor_est39)
    }


  } else if (S %% 6 == 2) {

    # Split  panel data for T equivalent to 2 modulo 6
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

    # Estimated quantities for split panel data
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

    # Function for bootstrap confidence interval
    mean_ci_func <- Vectorize(FUN = function(x) {

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(mean_est,
                                           mean_est21,
                                           mean_est22,
                                           mean_est31,
                                           mean_est32,
                                           mean_est33,
                                           mean_est34,
                                           mean_est35,
                                           mean_est36,
                                           mean_est37,
                                           mean_est38,
                                           mean_est39) - x,
                              statistic = toj2_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acov_est,
                                           acov_est21,
                                           acov_est22,
                                           acov_est31,
                                           acov_est32,
                                           acov_est33,
                                           acov_est34,
                                           acov_est35,
                                           acov_est36,
                                           acov_est37,
                                           acov_est38,
                                           acov_est39) - x,
                              statistic = toj2_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acor_est,
                                           acor_est21,
                                           acor_est22,
                                           acor_est31,
                                           acor_est32,
                                           acor_est33,
                                           acor_est34,
                                           acor_est35,
                                           acor_est36,
                                           acor_est37,
                                           acor_est38,
                                           acor_est39) - x,
                              statistic = toj2_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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
      mean_x <- seq(min(mean_est),
                    max(mean_est),
                    length = 1000)

      mean_y <- tojecdfest2(x   = mean_x,
                            X   = mean_est,
                            X21 = mean_est21,
                            X22 = mean_est22,
                            X31 = mean_est31,
                            X32 = mean_est32,
                            X33 = mean_est33,
                            X34 = mean_est34,
                            X35 = mean_est35,
                            X36 = mean_est36,
                            X37 = mean_est37,
                            X38 = mean_est38,
                            X39 = mean_est39)

      mean_toj <- Rearrangement::rearrangement(x = data.frame(x = mean_x),
                                               y = mean_y)

      mean_plot <- ggplot2::ggplot(data.frame(x = mean_x, y = mean_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(mean_est), max(mean_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::theme_bw()

      # Autocovariance
      acov_x <- seq(min(acov_est),
                    max(acov_est),
                    length = 1000)

      acov_y <- tojecdfest2(x   = acov_x,
                            X   = acov_est,
                            X21 = acov_est21,
                            X22 = acov_est22,
                            X31 = acov_est31,
                            X32 = acov_est32,
                            X33 = acov_est33,
                            X34 = acov_est34,
                            X35 = acov_est35,
                            X36 = acov_est36,
                            X37 = acov_est37,
                            X38 = acov_est38,
                            X39 = acov_est39)

      acov_toj <- Rearrangement::rearrangement(x = data.frame(x = acov_x),
                                               y = acov_y)

      acov_plot <- ggplot2::ggplot(data.frame(x = acov_x, y = acov_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acov_est), max(acov_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_x <- seq(min(acor_est),
                    max(acor_est),
                    length = 1000)

      acor_y <- tojecdfest2(x   = acor_x,
                            X   = acor_est,
                            X21 = acor_est21,
                            X22 = acor_est22,
                            X31 = acor_est31,
                            X32 = acor_est32,
                            X33 = acor_est33,
                            X34 = acor_est34,
                            X35 = acor_est35,
                            X36 = acor_est36,
                            X37 = acor_est37,
                            X38 = acor_est38,
                            X39 = acor_est39)

      acor_toj <- Rearrangement::rearrangement(x = data.frame(x = acor_x),
                                               y = acor_y)

      acor_plot <- ggplot2::ggplot(data.frame(x = acor_x, y = acor_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acor_est), max(acor_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::theme_bw()

    }

    if (ci) {

      # Mean
      mean_toj <- cbind(tojecdfest2(x   = mean_grid,
                                    X   = mean_est,
                                    X21 = mean_est21,
                                    X22 = mean_est22,
                                    X31 = mean_est31,
                                    X32 = mean_est32,
                                    X33 = mean_est33,
                                    X34 = mean_est34,
                                    X35 = mean_est35,
                                    X36 = mean_est36,
                                    X37 = mean_est37,
                                    X38 = mean_est38,
                                    X39 = mean_est39),
                        t(mean_ci))

      mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_grid),
                                   ggplot2::aes(x = mean_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(x = mean_grid,
                       y = Rearrangement::rearrangement(
                         x = data.frame(x = mean_grid),
                         y = mean_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(x = mean_grid,
                       ymin = Rearrangement::rearrangement(
                         x = data.frame(x = mean_grid),
                         y = mean_toj[, 2]),
                       ymax = Rearrangement::rearrangement(
                         x = data.frame(x = mean_grid),
                         y = mean_toj[, 3])
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocovariance
      acov_toj <- cbind(tojecdfest2(x   = acov_grid,
                                    X   = acov_est,
                                    X21 = acov_est21,
                                    X22 = acov_est22,
                                    X31 = acov_est31,
                                    X32 = acov_est32,
                                    X33 = acov_est33,
                                    X34 = acov_est34,
                                    X35 = acov_est35,
                                    X36 = acov_est36,
                                    X37 = acov_est37,
                                    X38 = acov_est38,
                                    X39 = acov_est39),
                        t(acov_ci))

      acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_grid),
                                   ggplot2::aes(x = acov_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(x = acov_grid,
                       y = Rearrangement::rearrangement(
                         x = data.frame(x = acov_grid),
                         y = acov_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(x = acov_grid,
                       ymin = Rearrangement::rearrangement(
                         x = data.frame(x = acov_grid),
                         y = acov_toj[, 2]),
                       ymax = Rearrangement::rearrangement(
                         x = data.frame(x = acov_grid),
                         y = acov_toj[, 3])
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_toj <- cbind(tojecdfest2(x   = acor_grid,
                                    X   = acor_est,
                                    X21 = acor_est21,
                                    X22 = acor_est22,
                                    X31 = acor_est31,
                                    X32 = acor_est32,
                                    X33 = acor_est33,
                                    X34 = acor_est34,
                                    X35 = acor_est35,
                                    X36 = acor_est36,
                                    X37 = acor_est37,
                                    X38 = acor_est38,
                                    X39 = acor_est39),
                        t(acor_ci))

      acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_grid),
                                   ggplot2::aes(x = acor_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(x = acor_grid,
                       y = Rearrangement::rearrangement(
                         x = data.frame(x = acor_grid),
                         y = acor_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(x = acor_grid,
                       ymin = Rearrangement::rearrangement(
                         x = data.frame(x = acor_grid),
                         y = acor_toj[, 2]),
                       ymax = Rearrangement::rearrangement(
                         x = data.frame(x = acor_grid),
                         y = acor_toj[, 3])
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

    }

    # Functions without rearrangement
    mean_func <- function(x) {
      tojecdfest2(x   = x,
                  X   = mean_est,
                  X21 = mean_est21,
                  X22 = mean_est22,
                  X31 = mean_est31,
                  X32 = mean_est32,
                  X33 = mean_est33,
                  X34 = mean_est34,
                  X35 = mean_est35,
                  X36 = mean_est36,
                  X37 = mean_est37,
                  X38 = mean_est38,
                  X39 = mean_est39)
    }

    acov_func <- function(x) {
      tojecdfest2(x   = x,
                  X   = acov_est,
                  X21 = acov_est21,
                  X22 = acov_est22,
                  X31 = acov_est31,
                  X32 = acov_est32,
                  X33 = acov_est33,
                  X34 = acov_est34,
                  X35 = acov_est35,
                  X36 = acov_est36,
                  X37 = acov_est37,
                  X38 = acov_est38,
                  X39 = acov_est39)
    }

    acor_func <- function(x) {
      tojecdfest2(x   = x,
                  X   = acor_est,
                  X21 = acor_est21,
                  X22 = acor_est22,
                  X31 = acor_est31,
                  X32 = acor_est32,
                  X33 = acor_est33,
                  X34 = acor_est34,
                  X35 = acor_est35,
                  X36 = acor_est36,
                  X37 = acor_est37,
                  X38 = acor_est38,
                  X39 = acor_est39)
    }

  } else if (S %% 6 == 3) {

    # Split  panel data for T equivalent to 3 modulo 6
    data21 <- data[, 1:floor(S / 2)]
    data22 <- data[, (floor(S / 2) + 1):S]
    data23 <- data[, 1:ceiling(S / 2)]
    data24 <- data[, (ceiling(S / 2) + 1):S]
    data31 <- data[, 1:(S / 3)]
    data32 <- data[, (S / 3 + 1):(2*S / 3)]
    data33 <- data[, (2 * S / 3 + 1):S]

    # Estimated quantities for split panel data
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

    # Function for bootstrap confidence interval
    mean_ci_func <- Vectorize(FUN = function(x) {

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(mean_est,
                                           mean_est21,
                                           mean_est22,
                                           mean_est23,
                                           mean_est24,
                                           mean_est31,
                                           mean_est32,
                                           mean_est33) - x,
                              statistic = toj3_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acov_est,
                                           acov_est21,
                                           acov_est22,
                                           acov_est23,
                                           acov_est24,
                                           acov_est31,
                                           acov_est32,
                                           acov_est33) - x,
                              statistic = toj3_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acor_est,
                                           acor_est21,
                                           acor_est22,
                                           acor_est23,
                                           acor_est24,
                                           acor_est31,
                                           acor_est32,
                                           acor_est33) - x,
                              statistic = toj3_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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
      mean_x <- seq(min(mean_est),
                    max(mean_est),
                    length = 1000)

      mean_y <- tojecdfest3(x   = mean_x,
                            X   = mean_est,
                            X21 = mean_est21,
                            X22 = mean_est22,
                            X23 = mean_est23,
                            X24 = mean_est24,
                            X31 = mean_est31,
                            X32 = mean_est32,
                            X33 = mean_est33)

      mean_toj <- Rearrangement::rearrangement(x = data.frame(x = mean_x),
                                               y = mean_y)

      mean_plot <- ggplot2::ggplot(data.frame(x = mean_x, y = mean_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(mean_est), max(mean_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::theme_bw()

      # Autocovariance
      acov_x <- seq(min(acov_est),
                    max(acov_est),
                    length = 1000)

      acov_y <- tojecdfest3(x   = acov_x,
                            X   = acov_est,
                            X21 = acov_est21,
                            X22 = acov_est22,
                            X23 = acov_est23,
                            X24 = acov_est24,
                            X31 = acov_est31,
                            X32 = acov_est32,
                            X33 = acov_est33)

      acov_toj <- Rearrangement::rearrangement(x = data.frame(x = acov_x),
                                               y = acov_y)

      acov_plot <- ggplot2::ggplot(data.frame(x = acov_x, y = acov_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acov_est), max(acov_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_x <- seq(min(acor_est),
                    max(acor_est),
                    length = 1000)

      acor_y <- tojecdfest3(x   = acor_x,
                            X   = acor_est,
                            X21 = acor_est21,
                            X22 = acor_est22,
                            X23 = acor_est23,
                            X24 = acor_est24,
                            X31 = acor_est31,
                            X32 = acor_est32,
                            X33 = acor_est33)

      acor_toj <- Rearrangement::rearrangement(x = data.frame(x = acor_x),
                                               y = acor_y)

      acor_plot <- ggplot2::ggplot(data.frame(x = acor_x, y = acor_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acor_est), max(acor_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::theme_bw()

    }

    if (ci) {

      # Mean
      mean_toj <- cbind(tojecdfest3(x   = mean_grid,
                                    X   = mean_est,
                                    X21 = mean_est21,
                                    X22 = mean_est22,
                                    X23 = mean_est23,
                                    X24 = mean_est24,
                                    X31 = mean_est31,
                                    X32 = mean_est32,
                                    X33 = mean_est33),
                        t(mean_ci))

      mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_grid),
                                   ggplot2::aes(x = mean_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(x = mean_grid,
                       y = Rearrangement::rearrangement(
                         x = data.frame(x = mean_grid),
                         y = mean_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = mean_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_toj[, 2]),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_toj[, 3])
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocovariance
      acov_toj <- cbind(tojecdfest3(x   = acov_grid,
                                    X   = acov_est,
                                    X21 = acov_est21,
                                    X22 = acov_est22,
                                    X23 = acov_est23,
                                    X24 = acov_est24,
                                    X31 = acov_est31,
                                    X32 = acov_est32,
                                    X33 = acov_est33),
                        t(acov_ci))

      acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_grid),
                                   ggplot2::aes(x = acov_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(x = acov_grid,
                       y = Rearrangement::rearrangement(
                         x = data.frame(x = acov_grid),
                         y = acov_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = acov_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_toj[, 2]),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_toj[, 3])
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_toj <- cbind(tojecdfest3(x   = acor_grid,
                                    X   = acor_est,
                                    X21 = acor_est21,
                                    X22 = acor_est22,
                                    X23 = acor_est23,
                                    X24 = acor_est24,
                                    X31 = acor_est31,
                                    X32 = acor_est32,
                                    X33 = acor_est33),
                        t(acor_ci))

      acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_grid),
                                   ggplot2::aes(x = acor_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(x = acor_grid,
                       y = Rearrangement::rearrangement(
                         x = data.frame(x = acor_grid),
                         y = acor_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = acor_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_toj[, 2]),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_toj[, 3])
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

    }

    # Functions without rearrangement
    mean_func <- function(x) {
      tojecdfest3(x   = x,
                  X   = mean_est,
                  X21 = mean_est21,
                  X22 = mean_est22,
                  X23 = mean_est23,
                  X24 = mean_est24,
                  X31 = mean_est31,
                  X32 = mean_est32,
                  X33 = mean_est33)
    }

    acov_func <- function(x) {
      tojecdfest3(x   = x,
                  X   = acov_est,
                  X21 = acov_est21,
                  X22 = acov_est22,
                  X23 = acov_est23,
                  X24 = acov_est24,
                  X31 = acov_est31,
                  X32 = acov_est32,
                  X33 = acov_est33)
    }

    acor_func <- function(x) {
      tojecdfest3(x   = x,
                  X   = acor_est,
                  X21 = acor_est21,
                  X22 = acor_est22,
                  X23 = acor_est23,
                  X24 = acor_est24,
                  X31 = acor_est31,
                  X32 = acor_est32,
                  X33 = acor_est33)
    }

  } else if (S %% 6 == 4) {

    # Split panel data for T equivalent to 4 modulo 6
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

    # Estimated quantities for split panel data
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

    # Function for bootstrap confidence interval
    mean_ci_func <- Vectorize(FUN = function(x) {

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(mean_est,
                                           mean_est21,
                                           mean_est22,
                                           mean_est31,
                                           mean_est32,
                                           mean_est33,
                                           mean_est34,
                                           mean_est35,
                                           mean_est36,
                                           mean_est37,
                                           mean_est38,
                                           mean_est39) - x,
                              statistic = toj4_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acov_est,
                                           acov_est21,
                                           acov_est22,
                                           acov_est31,
                                           acov_est32,
                                           acov_est33,
                                           acov_est34,
                                           acov_est35,
                                           acov_est36,
                                           acov_est37,
                                           acov_est38,
                                           acov_est39) - x,
                              statistic = toj4_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acor_est,
                                           acor_est21,
                                           acor_est22,
                                           acor_est31,
                                           acor_est32,
                                           acor_est33,
                                           acor_est34,
                                           acor_est35,
                                           acor_est36,
                                           acor_est37,
                                           acor_est38,
                                           acor_est39) - x,
                              statistic = toj4_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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
      mean_x <- seq(min(mean_est),
                    max(mean_est),
                    length = 1000)

      mean_y <- tojecdfest4(x   = mean_x,
                            X   = mean_est,
                            X21 = mean_est21,
                            X22 = mean_est22,
                            X31 = mean_est31,
                            X32 = mean_est32,
                            X33 = mean_est33,
                            X34 = mean_est34,
                            X35 = mean_est35,
                            X36 = mean_est36,
                            X37 = mean_est37,
                            X38 = mean_est38,
                            X39 = mean_est39)

      mean_toj <- Rearrangement::rearrangement(x = data.frame(x = mean_x),
                                               y = mean_y)

      mean_plot <- ggplot2::ggplot(data.frame(x = mean_x, y = mean_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(mean_est), max(mean_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::theme_bw()

      # Autocovariance
      acov_x <- seq(min(acov_est),
                    max(acov_est),
                    length = 1000)

      acov_y <- tojecdfest4(x   = acov_x,
                            X   = acov_est,
                            X21 = acov_est21,
                            X22 = acov_est22,
                            X31 = acov_est31,
                            X32 = acov_est32,
                            X33 = acov_est33,
                            X34 = acov_est34,
                            X35 = acov_est35,
                            X36 = acov_est36,
                            X37 = acov_est37,
                            X38 = acov_est38,
                            X39 = acov_est39)

      acov_toj <- Rearrangement::rearrangement(x = data.frame(x = acov_x),
                                               y = acov_y)

      acov_plot <- ggplot2::ggplot(data.frame(x = acov_x, y = acov_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acov_est), max(acov_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_x <- seq(min(acor_est),
                    max(acor_est),
                    length = 1000)

      acor_y <- tojecdfest4(x   = acor_x,
                            X   = acor_est,
                            X21 = acor_est21,
                            X22 = acor_est22,
                            X31 = acor_est31,
                            X32 = acor_est32,
                            X33 = acor_est33,
                            X34 = acor_est34,
                            X35 = acor_est35,
                            X36 = acor_est36,
                            X37 = acor_est37,
                            X38 = acor_est38,
                            X39 = acor_est39)

      acor_toj <- Rearrangement::rearrangement(x = data.frame(x = acor_x),
                                               y = acor_y)

      acor_plot <- ggplot2::ggplot(data.frame(x = acor_x, y = acor_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acor_est), max(acor_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::theme_bw()

    }

    if (ci) {

      # Mean
      mean_toj <- cbind(tojecdfest4(x   = mean_grid,
                                    X   = mean_est,
                                    X21 = mean_est21,
                                    X22 = mean_est22,
                                    X31 = mean_est31,
                                    X32 = mean_est32,
                                    X33 = mean_est33,
                                    X34 = mean_est34,
                                    X35 = mean_est35,
                                    X36 = mean_est36,
                                    X37 = mean_est37,
                                    X38 = mean_est38,
                                    X39 = mean_est39),
                        t(mean_ci))

      mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_grid),
                                   ggplot2::aes(x = mean_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(
            x = mean_grid,
            y = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = mean_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_toj[, 2]),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_toj[, 3])
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocovariance
      acov_toj <- cbind(tojecdfest4(x   = acov_grid,
                                    X   = acov_est,
                                    X21 = acov_est21,
                                    X22 = acov_est22,
                                    X31 = acov_est31,
                                    X32 = acov_est32,
                                    X33 = acov_est33,
                                    X34 = acov_est34,
                                    X35 = acov_est35,
                                    X36 = acov_est36,
                                    X37 = acov_est37,
                                    X38 = acov_est38,
                                    X39 = acov_est39),
                        t(acov_ci))

      acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_grid),
                                   ggplot2::aes(x = acov_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(
            x = acov_grid,
            y = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = acov_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_toj[, 2]),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_toj[, 3])
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_toj <- cbind(tojecdfest4(x   = acor_grid,
                                    X   = acor_est,
                                    X21 = acor_est21,
                                    X22 = acor_est22,
                                    X31 = acor_est31,
                                    X32 = acor_est32,
                                    X33 = acor_est33,
                                    X34 = acor_est34,
                                    X35 = acor_est35,
                                    X36 = acor_est36,
                                    X37 = acor_est37,
                                    X38 = acor_est38,
                                    X39 = acor_est39),
                        t(acor_ci))

      acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_grid),
                                   ggplot2::aes(x = acor_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(
            x = acor_grid,
            y = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = acor_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_toj[, 2]),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_toj[, 3])
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

    }

    # Functions without rearrangement
    mean_func <- function(x) {
      tojecdfest4(x   = x,
                  X   = mean_est,
                  X21 = mean_est21,
                  X22 = mean_est22,
                  X31 = mean_est31,
                  X32 = mean_est32,
                  X33 = mean_est33,
                  X34 = mean_est34,
                  X35 = mean_est35,
                  X36 = mean_est36,
                  X37 = mean_est37,
                  X38 = mean_est38,
                  X39 = mean_est39)
    }

    acov_func <- function(x) {
      tojecdfest4(x   = x,
                  X   = acov_est,
                  X21 = acov_est21,
                  X22 = acov_est22,
                  X31 = acov_est31,
                  X32 = acov_est32,
                  X33 = acov_est33,
                  X34 = acov_est34,
                  X35 = acov_est35,
                  X36 = acov_est36,
                  X37 = acov_est37,
                  X38 = acov_est38,
                  X39 = acov_est39)
    }

    acor_func <- function(x) {
      tojecdfest4(x   = x,
                  X   = acor_est,
                  X21 = acor_est21,
                  X22 = acor_est22,
                  X31 = acor_est31,
                  X32 = acor_est32,
                  X33 = acor_est33,
                  X34 = acor_est34,
                  X35 = acor_est35,
                  X36 = acor_est36,
                  X37 = acor_est37,
                  X38 = acor_est38,
                  X39 = acor_est39)
    }

  } else {

    # Split  panel data for T equivalent to 5 modulo 6
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

    # Estimated quantities for split panel data
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

    # Function for bootstrap confidence interval
    mean_ci_func <- Vectorize(FUN = function(x) {

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(mean_est,
                                           mean_est21,
                                           mean_est22,
                                           mean_est23,
                                           mean_est24,
                                           mean_est31,
                                           mean_est32,
                                           mean_est33,
                                           mean_est34,
                                           mean_est35,
                                           mean_est36,
                                           mean_est37,
                                           mean_est38,
                                           mean_est39) - x,
                              statistic = toj5_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acov_est,
                                           acov_est21,
                                           acov_est22,
                                           acov_est23,
                                           acov_est24,
                                           acov_est31,
                                           acov_est32,
                                           acov_est33,
                                           acov_est34,
                                           acov_est35,
                                           acov_est36,
                                           acov_est37,
                                           acov_est38,
                                           acov_est39) - x,
                              statistic = toj5_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

      # Confidence interval
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

      # TOJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acor_est,
                                           acor_est21,
                                           acor_est22,
                                           acor_est23,
                                           acor_est24,
                                           acor_est31,
                                           acor_est32,
                                           acor_est33,
                                           acor_est34,
                                           acor_est35,
                                           acor_est36,
                                           acor_est37,
                                           acor_est38,
                                           acor_est39) - x,
                              statistic = toj5_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t - bootstrap$t0

      # Confidence interval
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

    # Making figures using ggplot2
    if (!ci) {

      # Mean
      mean_x <- seq(min(mean_est),
                    max(mean_est),
                    length = 1000)

      mean_y <- tojecdfest5(x   = mean_x,
                            X   = mean_est,
                            X21 = mean_est21,
                            X22 = mean_est22,
                            X23 = mean_est23,
                            X24 = mean_est24,
                            X31 = mean_est31,
                            X32 = mean_est32,
                            X33 = mean_est33,
                            X34 = mean_est34,
                            X35 = mean_est35,
                            X36 = mean_est36,
                            X37 = mean_est37,
                            X38 = mean_est38,
                            X39 = mean_est39)

      mean_toj <- Rearrangement::rearrangement(x = data.frame(x = mean_x),
                                               y = mean_y)

      mean_plot <- ggplot2::ggplot(data.frame(x = mean_x, y = mean_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(mean_est), max(mean_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::theme_bw()

      # Autocovariance
      acov_x <- seq(min(acov_est),
                    max(acov_est),
                    length = 1000)

      acov_y <- tojecdfest5(x   = acov_x,
                            X   = acov_est,
                            X21 = acov_est21,
                            X22 = acov_est22,
                            X23 = acov_est23,
                            X24 = acov_est24,
                            X31 = acov_est31,
                            X32 = acov_est32,
                            X33 = acov_est33,
                            X34 = acov_est34,
                            X35 = acov_est35,
                            X36 = acov_est36,
                            X37 = acov_est37,
                            X38 = acov_est38,
                            X39 = acov_est39)

      acov_toj <- Rearrangement::rearrangement(x = data.frame(x = acov_x),
                                               y = acov_y)

      acov_plot <- ggplot2::ggplot(data.frame(x = acov_x, y = acov_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acov_est), max(acov_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_x <- seq(min(acor_est),
                    max(acor_est),
                    length = 1000)

      acor_y <- tojecdfest5(x   = acor_x,
                            X   = acor_est,
                            X21 = acor_est21,
                            X22 = acor_est22,
                            X23 = acor_est23,
                            X24 = acor_est24,
                            X31 = acor_est31,
                            X32 = acor_est32,
                            X33 = acor_est33,
                            X34 = acor_est34,
                            X35 = acor_est35,
                            X36 = acor_est36,
                            X37 = acor_est37,
                            X38 = acor_est38,
                            X39 = acor_est39)

      acor_toj <- Rearrangement::rearrangement(x = data.frame(x = acor_x),
                                               y = acor_y)

      acor_plot <- ggplot2::ggplot(data.frame(x = acor_x, y = acor_toj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acor_est), max(acor_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::theme_bw()

    }

    if (ci) {

      # Mean
      mean_toj <- cbind(tojecdfest5(x   = mean_grid,
                                    X   = mean_est,
                                    X21 = mean_est21,
                                    X22 = mean_est22,
                                    X23 = mean_est23,
                                    X24 = mean_est24,
                                    X31 = mean_est31,
                                    X32 = mean_est32,
                                    X33 = mean_est33,
                                    X34 = mean_est34,
                                    X35 = mean_est35,
                                    X36 = mean_est36,
                                    X37 = mean_est37,
                                    X38 = mean_est38,
                                    X39 = mean_est39),
                        t(mean_ci))

      mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_grid),
                                   ggplot2::aes(x = mean_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(x = mean_grid,
                       y = Rearrangement::rearrangement(
                         x = data.frame(x = mean_grid),
                         y = mean_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = mean_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_toj[, 2]),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_toj[, 3])
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocovariance
      acov_toj <- cbind(tojecdfest5(x   = acov_grid,
                                    X   = acov_est,
                                    X21 = acov_est21,
                                    X22 = acov_est22,
                                    X23 = acov_est23,
                                    X24 = acov_est24,
                                    X31 = acov_est31,
                                    X32 = acov_est32,
                                    X33 = acov_est33,
                                    X34 = acov_est34,
                                    X35 = acov_est35,
                                    X36 = acov_est36,
                                    X37 = acov_est37,
                                    X38 = acov_est38,
                                    X39 = acov_est39),
                        t(acov_ci))

      acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_grid),
                                   ggplot2::aes(x = acov_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(x = acov_grid,
                       y = Rearrangement::rearrangement(
                         x = data.frame(x = acov_grid),
                         y = acov_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = acov_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_toj[, 2]),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_toj[, 3])
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_toj <- cbind(tojecdfest5(x   = acor_grid,
                                    X   = acor_est,
                                    X21 = acor_est21,
                                    X22 = acor_est22,
                                    X23 = acor_est23,
                                    X24 = acor_est24,
                                    X31 = acor_est31,
                                    X32 = acor_est32,
                                    X33 = acor_est33,
                                    X34 = acor_est34,
                                    X35 = acor_est35,
                                    X36 = acor_est36,
                                    X37 = acor_est37,
                                    X38 = acor_est38,
                                    X39 = acor_est39),
                        t(acor_ci))

      acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_grid),
                                   ggplot2::aes(x = acor_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(x = acor_grid,
                       y = Rearrangement::rearrangement(
                         x = data.frame(x = acor_grid),
                         y = acor_toj[, 1])
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = acor_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_toj[, 2]),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_toj[, 3])
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

    }

    # Functions without rearrangement
    mean_func <- function(x) {
      tojecdfest5(x   = x,
                  X   = mean_est,
                  X21 = mean_est21,
                  X22 = mean_est22,
                  X23 = mean_est23,
                  X24 = mean_est24,
                  X31 = mean_est31,
                  X32 = mean_est32,
                  X33 = mean_est33,
                  X34 = mean_est34,
                  X35 = mean_est35,
                  X36 = mean_est36,
                  X37 = mean_est37,
                  X38 = mean_est38,
                  X39 = mean_est39)
    }

    acov_func <- function(x) {
      tojecdfest5(x   = x,
                  X   = acov_est,
                  X21 = acov_est21,
                  X22 = acov_est22,
                  X23 = acov_est23,
                  X24 = acov_est24,
                  X31 = acov_est31,
                  X32 = acov_est32,
                  X33 = acov_est33,
                  X34 = acov_est34,
                  X35 = acov_est35,
                  X36 = acov_est36,
                  X37 = acov_est37,
                  X38 = acov_est38,
                  X39 = acov_est39)
    }

    acor_func <- function(x) {
      tojecdfest5(x   = x,
                  X   = acor_est,
                  X21 = acor_est21,
                  X22 = acor_est22,
                  X23 = acor_est23,
                  X24 = acor_est24,
                  X31 = acor_est31,
                  X32 = acor_est32,
                  X33 = acor_est33,
                  X34 = acor_est34,
                  X35 = acor_est35,
                  X36 = acor_est36,
                  X37 = acor_est37,
                  X38 = acor_est38,
                  X39 = acor_est39)
    }

  }

  # Results
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

#' Compute TOJ empirical CDF estimate for T equivalent to 0 modulo 6
#'
#' @param x An evaluation point
#' @param X A vector of original cross-sectional data
#' @param X21 A vector of half-panel cross-sectional data 1
#' @param X22 A vector of half-panel cross-sectional data 2
#' @param X31 A vector of third-panel cross-sectional data 1
#' @param X32 A vector of third-panel cross-sectional data 2
#' @param X33 A vector of third-panel cross-sectional data 3
#'
#' @noRd
#'
tojecdfest0 <- Vectorize(FUN = function(x, X, X21, X22, X31, X32, X33) {

  # Estimates
  est   <- mean(X   <= x)
  est21 <- mean(X21 <= x)
  est22 <- mean(X22 <= x)
  est31 <- mean(X31 <= x)
  est32 <- mean(X32 <= x)
  est33 <- mean(X33 <= x)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22) / 2 +
    1.536 * (est31 + est32 + est33) / 3

  # Ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}, vectorize.args = "x")

#' Compute TOJ empirical CDF estimate for T equivalent to 1 modulo 6
#'
#' @param x An evaluation point
#' @param X A vector of original cross-sectional data
#' @param X21 A vector of half-panel cross-sectional data 1
#' @param X22 A vector of half-panel cross-sectional data 2
#' @param X23 A vector of half-panel cross-sectional data 3
#' @param X24 A vector of half-panel cross-sectional data 4
#' @param X31 A vector of third-panel cross-sectional data 1
#' @param X32 A vector of third-panel cross-sectional data 2
#' @param X33 A vector of third-panel cross-sectional data 3
#' @param X34 A vector of third-panel cross-sectional data 4
#' @param X35 A vector of third-panel cross-sectional data 5
#' @param X36 A vector of third-panel cross-sectional data 6
#' @param X37 A vector of third-panel cross-sectional data 7
#' @param X38 A vector of third-panel cross-sectional data 8
#' @param X39 A vector of third-panel cross-sectional data 9
#'
#' @noRd
#'
tojecdfest1 <- Vectorize(FUN = function(x,
                                        X,
                                        X21,
                                        X22,
                                        X23,
                                        X24,
                                        X31,
                                        X32,
                                        X33,
                                        X34,
                                        X35,
                                        X36,
                                        X37,
                                        X38,
                                        X39) {

  # Estimates
  est   <- mean(X   <= x)
  est21 <- mean(X21 <= x)
  est22 <- mean(X22 <= x)
  est23 <- mean(X23 <= x)
  est24 <- mean(X24 <= x)
  est31 <- mean(X31 <= x)
  est32 <- mean(X32 <= x)
  est33 <- mean(X33 <= x)
  est34 <- mean(X34 <= x)
  est35 <- mean(X35 <= x)
  est36 <- mean(X36 <= x)
  est37 <- mean(X37 <= x)
  est38 <- mean(X38 <= x)
  est39 <- mean(X39 <= x)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22 + est23 + est24) / 4 +
    1.536 * (est31 + est32 + est33 + est34 +
               est35 + est36 + est37 + est38 + est39) / 9

  # Ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}, vectorize.args = "x")

#' Compute TOJ empirical CDF estimate for T equivalent to 2 modulo 6
#'
#' @param x An evaluation point
#' @param X A vector of original cross-sectional data
#' @param X21 A vector of half-panel cross-sectional data 1
#' @param X22 A vector of half-panel cross-sectional data 2
#' @param X31 A vector of third-panel cross-sectional data 1
#' @param X32 A vector of third-panel cross-sectional data 2
#' @param X33 A vector of third-panel cross-sectional data 3
#' @param X34 A vector of third-panel cross-sectional data 4
#' @param X35 A vector of third-panel cross-sectional data 5
#' @param X36 A vector of third-panel cross-sectional data 6
#' @param X37 A vector of third-panel cross-sectional data 7
#' @param X38 A vector of third-panel cross-sectional data 8
#' @param X39 A vector of third-panel cross-sectional data 9
#'
#' @noRd
#'
tojecdfest2 <- Vectorize(FUN = function(x,
                                        X,
                                        X21,
                                        X22,
                                        X31,
                                        X32,
                                        X33,
                                        X34,
                                        X35,
                                        X36,
                                        X37,
                                        X38,
                                        X39) {

  # Estimates
  est   <- mean(X   <= x)
  est21 <- mean(X21 <= x)
  est22 <- mean(X22 <= x)
  est31 <- mean(X31 <= x)
  est32 <- mean(X32 <= x)
  est33 <- mean(X33 <= x)
  est34 <- mean(X34 <= x)
  est35 <- mean(X35 <= x)
  est36 <- mean(X36 <= x)
  est37 <- mean(X37 <= x)
  est38 <- mean(X38 <= x)
  est39 <- mean(X39 <= x)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22) / 2 +
    1.536 * (est31 + est32 + est33 + est34 +
               est35 + est36 + est37 + est38 + est39) / 9

  # Ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}, vectorize.args = "x")

#' Compute TOJ empirical CDF estimate for T equivalent to 3 modulo 6
#'
#' @param x An evaluation point
#' @param X A vector of original cross-sectional data
#' @param X21 A vector of half-panel cross-sectional data 1
#' @param X22 A vector of half-panel cross-sectional data 2
#' @param X23 A vector of half-panel cross-sectional data 3
#' @param X24 A vector of half-panel cross-sectional data 4
#' @param X31 A vector of third-panel cross-sectional data 1
#' @param X32 A vector of third-panel cross-sectional data 2
#' @param X33 A vector of third-panel cross-sectional data 3
#'
#' @noRd
#'
tojecdfest3 <- Vectorize(FUN = function(x,
                                        X,
                                        X21,
                                        X22,
                                        X23,
                                        X24,
                                        X31,
                                        X32,
                                        X33) {

  # Estimates
  est   <- mean(X   <= x)
  est21 <- mean(X21 <= x)
  est22 <- mean(X22 <= x)
  est23 <- mean(X23 <= x)
  est24 <- mean(X24 <= x)
  est31 <- mean(X31 <= x)
  est32 <- mean(X32 <= x)
  est33 <- mean(X33 <= x)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22 + est23 + est24) / 4 +
    1.536 * (est31 + est32 + est33) / 3

  # Ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}, vectorize.args = "x")

#' Compute TOJ empirical CDF estimate for T equivalent to 4 modulo 6
#'
#' @param x An evaluation point
#' @param X A vector of original cross-sectional data
#' @param X21 A vector of half-panel cross-sectional data 1
#' @param X22 A vector of half-panel cross-sectional data 2
#' @param X31 A vector of third-panel cross-sectional data 1
#' @param X32 A vector of third-panel cross-sectional data 2
#' @param X33 A vector of third-panel cross-sectional data 3
#' @param X34 A vector of third-panel cross-sectional data 4
#' @param X35 A vector of third-panel cross-sectional data 5
#' @param X36 A vector of third-panel cross-sectional data 6
#' @param X37 A vector of third-panel cross-sectional data 7
#' @param X38 A vector of third-panel cross-sectional data 8
#' @param X39 A vector of third-panel cross-sectional data 9
#'
#' @noRd
#'
tojecdfest4 <- Vectorize(FUN = function(x,
                                        X,
                                        X21,
                                        X22,
                                        X31,
                                        X32,
                                        X33,
                                        X34,
                                        X35,
                                        X36,
                                        X37,
                                        X38,
                                        X39) {

  # Estimates
  est   <- mean(X   <= x)
  est21 <- mean(X21 <= x)
  est22 <- mean(X22 <= x)
  est31 <- mean(X31 <= x)
  est32 <- mean(X32 <= x)
  est33 <- mean(X33 <= x)
  est34 <- mean(X34 <= x)
  est35 <- mean(X35 <= x)
  est36 <- mean(X36 <= x)
  est37 <- mean(X37 <= x)
  est38 <- mean(X38 <= x)
  est39 <- mean(X39 <= x)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22) / 2 +
    1.536 * (est31 + est32 + est33 + est34 +
               est35 + est36 + est37 + est38 + est39) / 9

  # Ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}, vectorize.args = "x")

#' Compute TOJ empirical CDF estimate for T equivalent to 5 modulo 6
#'
#' @param x An evaluation point
#' @param X A vector of original cross-sectional data
#' @param X21 A vector of half-panel cross-sectional data 1
#' @param X22 A vector of half-panel cross-sectional data 2
#' @param X23 A vector of half-panel cross-sectional data 3
#' @param X24 A vector of half-panel cross-sectional data 4
#' @param X31 A vector of third-panel cross-sectional data 1
#' @param X32 A vector of third-panel cross-sectional data 2
#' @param X33 A vector of third-panel cross-sectional data 3
#' @param X34 A vector of third-panel cross-sectional data 4
#' @param X35 A vector of third-panel cross-sectional data 5
#' @param X36 A vector of third-panel cross-sectional data 6
#' @param X37 A vector of third-panel cross-sectional data 7
#' @param X38 A vector of third-panel cross-sectional data 8
#' @param X39 A vector of third-panel cross-sectional data 9
#'
#' @noRd
#'
tojecdfest5 <- Vectorize(FUN = function(x,
                                        X,
                                        X21,
                                        X22,
                                        X23,
                                        X24,
                                        X31,
                                        X32,
                                        X33,
                                        X34,
                                        X35,
                                        X36,
                                        X37,
                                        X38,
                                        X39) {

  # Estimates
  est   <- mean(X   <= x)
  est21 <- mean(X21 <= x)
  est22 <- mean(X22 <= x)
  est23 <- mean(X23 <= x)
  est24 <- mean(X24 <= x)
  est31 <- mean(X31 <= x)
  est32 <- mean(X32 <= x)
  est33 <- mean(X33 <= x)
  est34 <- mean(X34 <= x)
  est35 <- mean(X35 <= x)
  est36 <- mean(X36 <= x)
  est37 <- mean(X37 <= x)
  est38 <- mean(X38 <= x)
  est39 <- mean(X39 <= x)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22 + est23 + est24) / 4 +
    1.536 * (est31 + est32 + est33 + est34 +
               est35 + est36 + est37 + est38 + est39) / 9

  # Ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}, vectorize.args = "x")

#' Compute bootstrap TOJ empirical CDF estimate for T equivalent to 0 modulo 6
#'
#' @param quantity An N * 6 matrix of estimates
#' @param indices A vector of indices for bootstrap repetitions
#'
#' @noRd
#'
toj0_boot <- function(quantity, indices) {

  # Estimates
  est   <- mean(quantity[indices, 1] <= 0)
  est21 <- mean(quantity[indices, 2] <= 0)
  est22 <- mean(quantity[indices, 3] <= 0)
  est31 <- mean(quantity[indices, 4] <= 0)
  est32 <- mean(quantity[indices, 5] <= 0)
  est33 <- mean(quantity[indices, 6] <= 0)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22) / 2 +
    1.536 * (est31 + est32 + est33) / 3

  # Ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)
}

#' Compute bootstrap TOJ empirical CDF estimate for T equivalent to 1 modulo 6
#'
#' @param quantity An N * 14 matrix of estimates
#' @param indices A vector of indices for bootstrap repetitions
#'
#' @noRd
#'
toj1_boot <- function(quantity, indices) {

  # Estimates
  est   <- mean(quantity[indices, 1]  <= 0)
  est21 <- mean(quantity[indices, 2]  <= 0)
  est22 <- mean(quantity[indices, 3]  <= 0)
  est23 <- mean(quantity[indices, 4]  <= 0)
  est24 <- mean(quantity[indices, 5]  <= 0)
  est31 <- mean(quantity[indices, 6]  <= 0)
  est32 <- mean(quantity[indices, 7]  <= 0)
  est33 <- mean(quantity[indices, 8]  <= 0)
  est34 <- mean(quantity[indices, 9]  <= 0)
  est35 <- mean(quantity[indices, 10] <= 0)
  est36 <- mean(quantity[indices, 11] <= 0)
  est37 <- mean(quantity[indices, 12] <= 0)
  est38 <- mean(quantity[indices, 13] <= 0)
  est39 <- mean(quantity[indices, 14] <= 0)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22 + est23 + est24) / 4 +
    1.536 * (est31 + est32 + est33 + est34 +
               est35 + est36 + est37 + est38 + est39) / 9

  # Ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)
}

#' Compute bootstrap TOJ empirical CDF estimate for T equivalent to 2 modulo 6
#'
#' @param quantity An N * 12 matrix of estimates
#' @param indices A vector of indices for bootstrap repetitions
#'
#' @noRd
#'
toj2_boot <- function(quantity, indices) {

  # Estimates
  est   <- mean(quantity[indices, 1]  <= 0)
  est21 <- mean(quantity[indices, 2]  <= 0)
  est22 <- mean(quantity[indices, 3]  <= 0)
  est31 <- mean(quantity[indices, 4]  <= 0)
  est32 <- mean(quantity[indices, 5]  <= 0)
  est33 <- mean(quantity[indices, 6]  <= 0)
  est34 <- mean(quantity[indices, 7]  <= 0)
  est35 <- mean(quantity[indices, 8]  <= 0)
  est36 <- mean(quantity[indices, 9]  <= 0)
  est37 <- mean(quantity[indices, 10] <= 0)
  est38 <- mean(quantity[indices, 11] <= 0)
  est39 <- mean(quantity[indices, 12] <= 0)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22) / 2 +
    1.536 * (est31 + est32 + est33 + est34 +
               est35 + est36 + est37 + est38 + est39) / 9

  # Ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}

#' Compute bootstrap TOJ empirical CDF estimate for T equivalent to 3 modulo 6
#'
#' @param quantity An N * 8 matrix of estimates
#' @param indices A vector of indices for bootstrap repetitions
#'
#' @noRd
#'
toj3_boot <- function(quantity, indices) {

  # Estimates
  est   <- mean(quantity[indices, 1] <= 0)
  est21 <- mean(quantity[indices, 2] <= 0)
  est22 <- mean(quantity[indices, 3] <= 0)
  est23 <- mean(quantity[indices, 4] <= 0)
  est24 <- mean(quantity[indices, 5] <= 0)
  est31 <- mean(quantity[indices, 6] <= 0)
  est32 <- mean(quantity[indices, 7] <= 0)
  est33 <- mean(quantity[indices, 8] <= 0)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22 + est23 + est24) / 4 +
    1.536 * (est31 + est32 + est33) / 3

  # Ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}

#' Compute bootstrap TOJ empirical CDF estimate for T equivalent to 4 modulo 6
#'
#' @param quantity An N * 12 matrix of estimates
#' @param indices A vector of indices for bootstrap repetitions
#'
#' @noRd
#'
toj4_boot <- function(quantity, indices) {

  # Estimates
  est   <- mean(quantity[indices, 1]  <= 0)
  est21 <- mean(quantity[indices, 2]  <= 0)
  est22 <- mean(quantity[indices, 3]  <= 0)
  est31 <- mean(quantity[indices, 4]  <= 0)
  est32 <- mean(quantity[indices, 5]  <= 0)
  est33 <- mean(quantity[indices, 6]  <= 0)
  est34 <- mean(quantity[indices, 7]  <= 0)
  est35 <- mean(quantity[indices, 8]  <= 0)
  est36 <- mean(quantity[indices, 9]  <= 0)
  est37 <- mean(quantity[indices, 10] <= 0)
  est38 <- mean(quantity[indices, 11] <= 0)
  est39 <- mean(quantity[indices, 12] <= 0)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22) / 2 +
    1.536 * (est31 + est32 + est33 + est34 +
               est35 + est36 + est37 + est38 + est39) / 9

  # Ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}

#' Compute bootstrap TOJ empirical CDF estimate for T equivalent to 5 modulo 6
#'
#' @param quantity N * 14 matrix of estimates
#' @param indices A vector of indices for bootstrap repetitions
#'
#' @noRd
#'
toj5_boot <- function(quantity, indices) {

  # Estimates
  est   <- mean(quantity[indices, 1]  <= 0)
  est21 <- mean(quantity[indices, 2]  <= 0)
  est22 <- mean(quantity[indices, 3]  <= 0)
  est23 <- mean(quantity[indices, 4]  <= 0)
  est24 <- mean(quantity[indices, 5]  <= 0)
  est31 <- mean(quantity[indices, 6]  <= 0)
  est32 <- mean(quantity[indices, 7]  <= 0)
  est33 <- mean(quantity[indices, 8]  <= 0)
  est34 <- mean(quantity[indices, 9]  <= 0)
  est35 <- mean(quantity[indices, 10] <= 0)
  est36 <- mean(quantity[indices, 11] <= 0)
  est37 <- mean(quantity[indices, 12] <= 0)
  est38 <- mean(quantity[indices, 13] <= 0)
  est39 <- mean(quantity[indices, 14] <= 0)

  # TOJ estimate
  tojest <- 3.536 * est -
    4.072 * (est21 + est22 + est23 + est24) / 4 +
    1.536 * (est31 + est32 + est33 + est34 +
               est35 + est36 + est37 + est38 + est39) / 9

  # Ensure valid estimates
  tojest <- ifelse(tojest >= 0, tojest, 0)
  tojest <- ifelse(tojest <= 1, tojest, 1)

  return(tojest)

}