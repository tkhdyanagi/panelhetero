#' The HPJ bias-corrected empirical CDF estimation
#'
#' The `hpjecdf()` function enables to implement the HPJ bias-corrected
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
#' panelhetero::hpjecdf(data = data, R = 50)
#'
#' @references Okui, R. and Yanagi, T., 2019.
#' Panel data analysis with heterogeneous dynamics.
#' Journal of Econometrics, 212(2), pp.451-475.
#'
#' @export
#'
hpjecdf <- function(data,
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

  # HPJ bias-corrected estimation ----------------------------------------------

  if (S %% 2 == 0) {

    # Half panel data for even T
    data1 <- data[, 1:(S / 2)]
    data2 <- data[, (S / 2 + 1):S]

    # Estimated quantities for half panel data
    mean_est1 <- rowMeans(data1)
    mean_est2 <- rowMeans(data2)

    acov_est1 <- apply(data1, MARGIN = 1, acov, acov_order = acov_order)
    acov_est2 <- apply(data2, MARGIN = 1, acov, acov_order = acov_order)

    acor_est1 <- apply(data1, MARGIN = 1, acor, acor_order = acor_order)
    acor_est2 <- apply(data2, MARGIN = 1, acor, acor_order = acor_order)

    # Function for bootstrap confidence interval
    mean_ci_func <- Vectorize(FUN = function(x) {

      # HPJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(mean_est,
                                           mean_est1,
                                           mean_est2) - x,
                              statistic = hpj1_boot,
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

      # HPJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acov_est,
                                           acov_est1,
                                           acov_est2) - x,
                              statistic = hpj1_boot,
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

      # HPJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acor_est,
                                           acor_est1,
                                           acor_est2) - x,
                              statistic = hpj1_boot,
                              R = R)
      estimate <- bootstrap$t0
      temp <- bootstrap$t- bootstrap$t0

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
      mean_x <- seq(min(mean_est),
                    max(mean_est),
                    length = 1000)

      mean_y <- hpjecdfest1(x = mean_x,
                            X = mean_est,
                            X1 = mean_est1,
                            X2 = mean_est2)

      mean_hpj <- Rearrangement::rearrangement(x = data.frame(x = mean_x),
                                               y = mean_y)

      mean_plot <- ggplot2::ggplot(data.frame(x = mean_x,
                                              y = mean_hpj),
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

      acov_y <- hpjecdfest1(x = acov_x,
                            X = acov_est,
                            X1 = acov_est1,
                            X2 = acov_est2)

      acov_hpj <- Rearrangement::rearrangement(x = data.frame(x = acov_x),
                                               y = acov_y)

      acov_plot <- ggplot2::ggplot(data.frame(x = acov_x,
                                              y = acov_hpj),
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

      acor_y <- hpjecdfest1(x = acor_x,
                            X = acor_est,
                            X1 = acor_est1,
                            X2 = acor_est2)

      acor_hpj <- Rearrangement::rearrangement(x = data.frame(x = acor_x),
                                               y = acor_y)

      acor_plot <- ggplot2::ggplot(data.frame(x = acor_x,
                                              y = acor_hpj),
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
      mean_hpj <- cbind(hpjecdfest1(x  = mean_grid,
                                    X  = mean_est,
                                    X1 = mean_est1,
                                    X2 = mean_est2),
                        t(mean_ci))

      mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_grid),
                                   ggplot2::aes(x = mean_grid)) +
        ggplot2::geom_line(ggplot2::aes(
          x = mean_grid,
          y = Rearrangement::rearrangement(
            x = data.frame(x = mean_grid),
            y = mean_hpj[, 1])
        )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = mean_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_hpj[, 2]
            ),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_hpj[, 3]
            )
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocovariance
      acov_hpj <- cbind(hpjecdfest1(x  = acov_grid,
                                    X  = acov_est,
                                    X1 = acov_est1,
                                    X2 = acov_est2),
                        t(acov_ci))

      acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_grid),
                                   ggplot2::aes(x = acov_grid)) +
        ggplot2::geom_line(ggplot2::aes(
          x = acov_grid,
          y = Rearrangement::rearrangement(
            x = data.frame(x = acov_grid),
            y = acov_hpj[, 1])
        )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = acov_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_hpj[, 2]
            ),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_hpj[, 3]
            )
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_hpj <- cbind(hpjecdfest1(x  = acor_grid,
                                    X  = acor_est,
                                    X1 = acor_est1,
                                    X2 = acor_est2),
                        t(acor_ci))

      acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_grid),
                                   ggplot2::aes(x = acor_grid)) +
        ggplot2::geom_line(ggplot2::aes(
          x = acor_grid,
          y = Rearrangement::rearrangement(
            x = data.frame(x = acor_grid),
            y = acor_hpj[, 1])
        )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = acor_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_hpj[, 2]
            ),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_hpj[, 3]
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
      hpjecdfest1(x = x, X = mean_est, X1 = mean_est1, X2 = mean_est2)
    }

    acov_func <- function(x) {
      hpjecdfest1(x = x, X = acov_est, X1 = acov_est1, X2 = acov_est2)
    }

    acor_func <- function(x) {
      hpjecdfest1(x = x, X = acor_est, X1 = acor_est1, X2 = acor_est2)
    }

  } else {

    # Half panel data for odd T
    data1 <- data[, 1:floor(S / 2)]
    data2 <- data[, (floor(S / 2) + 1):S]
    data3 <- data[, 1:ceiling(S / 2)]
    data4 <- data[, (ceiling(S / 2) + 1):S]

    # Estimated quantities for half panel data
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

    # Function for bootstrap confidence interval
    mean_ci_func <- Vectorize(FUN = function(x) {

      # HPJ bias-corretced estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(mean_est,
                                           mean_est1,
                                           mean_est2,
                                           mean_est3,
                                           mean_est4) - x,
                              statistic = hpj2_boot,
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

      # HPJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acov_est,
                                           acov_est1,
                                           acov_est2,
                                           acov_est3,
                                           acov_est4) - x,
                              statistic = hpj2_boot,
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

      # HPJ bias-corrected estimation with bootstrap
      bootstrap <- boot::boot(data = cbind(acor_est,
                                           acor_est1,
                                           acor_est2,
                                           acor_est3,
                                           acor_est4) - x,
                              statistic = hpj2_boot,
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
      mean_x <- seq(min(mean_est), max(mean_est), length = 1000)

      mean_y <- hpjecdfest2(x  = mean_x,
                            X  = mean_est,
                            X1 = mean_est1,
                            X2 = mean_est2,
                            X3 = mean_est3,
                            X4 = mean_est4)

      mean_hpj <- Rearrangement::rearrangement(
        x = data.frame(x = mean_x),
        y = mean_y
      )

      mean_plot <- ggplot2::ggplot(data.frame(x = mean_x, y = mean_hpj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(mean_est), max(mean_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::theme_bw()

      # Autocovariance
      acov_x <- seq(min(acov_est), max(acov_est), length = 1000)

      acov_y <- hpjecdfest2(x  = acov_x,
                            X  = acov_est,
                            X1 = acov_est1,
                            X2 = acov_est2,
                            X3 = acov_est3,
                            X4 = acov_est4)

      acov_hpj <- Rearrangement::rearrangement(
        x = data.frame(x = acov_x),
        y = acov_y
      )

      acov_plot <- ggplot2::ggplot(data.frame(x = acov_x, y = acov_hpj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acov_est), max(acov_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_x <- seq(min(acor_est), max(acor_est), length = 1000)

      acor_y <- hpjecdfest2(x  = acor_x,
                            X  = acor_est,
                            X1 = acor_est1,
                            X2 = acor_est2,
                            X3 = acor_est3,
                            X4 = acor_est4)

      acor_hpj <- Rearrangement::rearrangement(
        x = data.frame(x = acor_x),
        y = acor_y
      )

      acor_plot <- ggplot2::ggplot(data.frame(x = acor_x, y = acor_hpj),
                                   ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line() +
        ggplot2::xlim(min(acor_est), max(acor_est)) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocorrelation") +
        ggplot2::theme_bw()

    }

    if (ci){

      # Mean
      mean_hpj <- cbind(hpjecdfest2(x  = mean_grid,
                                    X  = mean_est,
                                    X1 = mean_est1,
                                    X2 = mean_est2,
                                    X3 = mean_est3,
                                    X4 = mean_est4),
                        t(mean_ci))

      mean_plot <- ggplot2::ggplot(data = data.frame(x = mean_grid),
                                   ggplot2::aes(x = mean_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(
            x = mean_grid,
            y = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_hpj[, 1]
            )
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = mean_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_hpj[, 2]
            ),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = mean_grid),
              y = mean_hpj[, 3]
            )
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous mean") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocovariance
      acov_hpj <- cbind(hpjecdfest2(x  = acov_grid,
                                    X  = acov_est,
                                    X1 = acov_est1,
                                    X2 = acov_est2,
                                    X3 = acov_est3,
                                    X4 = acov_est4),
                        t(acov_ci))

      acov_plot <- ggplot2::ggplot(data = data.frame(x = acov_grid),
                                   ggplot2::aes(x = acov_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(
            x = acov_grid,
            y = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_hpj[, 1]
            )
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = acov_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_hpj[, 2]
            ),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = acov_grid),
              y = acov_hpj[, 3]
            )
          ),
          alpha = 0.1) +
        ggplot2::labs(x = "x", y = "") +
        ggplot2::ggtitle("The heterogeneous autocovariance") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw()

      # Autocorrelation
      acor_hpj <- cbind(hpjecdfest2(x  = acor_grid,
                                    X  = acor_est,
                                    X1 = acor_est1,
                                    X2 = acor_est2,
                                    X3 = acor_est3,
                                    X4 = acor_est4),
                        t(acor_ci))

      acor_plot <- ggplot2::ggplot(data = data.frame(x = acor_grid),
                                   ggplot2::aes(x = acor_grid)) +
        ggplot2::geom_line(
          ggplot2::aes(
            x = acor_grid,
            y = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_hpj[, 1]
            )
          )
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = acor_grid,
            ymin = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_hpj[, 2]
            ),
            ymax = Rearrangement::rearrangement(
              x = data.frame(x = acor_grid),
              y = acor_hpj[, 3]
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
      hpjecdfest2(x  = x,
                  X  = mean_est,
                  X1 = mean_est1,
                  X2 = mean_est2,
                  X3 = mean_est3,
                  X4 = mean_est4)
    }

    acov_func <- function(x) {
      hpjecdfest2(x  = x,
                  X  = acov_est,
                  X1 = acov_est1,
                  X2 = acov_est2,
                  X3 = acov_est3,
                  X4 = acov_est4)
    }

    acor_func <- function(x) {
      hpjecdfest2(x  = x,
                  X  = acor_est,
                  X1 = acor_est1,
                  X2 = acor_est2,
                  X3 = acor_est3,
                  X4 = acor_est4)
    }

  }

  # Result
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


#' Compute HPJ bias-corrected empirical CDF estimate for even T
#'
#' @param x An evaluation point
#' @param X A vector of original cross-sectional data
#' @param X1 A vector of half-panel cross-sectional data 1
#' @param X2 A vector of half-panel cross-sectional data 2
#'
#' @noRd
#'
hpjecdfest1 <- Vectorize(FUN = function(x, X, X1, X2) {

  # Estimates
  est  <- mean(X  <= x)
  est1 <- mean(X1 <= x)
  est2 <- mean(X2 <= x)

  # HPJ bias-corrected estimates
  hpjest <- 2 * est - (est1 + est2) / 2

  # Ensure valid estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)
  hpjest <- ifelse(hpjest <= 1, hpjest, 1)

  return(hpjest)

}, vectorize.args = "x")


#' Compute HPJ bias-corrected empirical CDF estimate for odd T
#'
#' @param x An evaluation point
#' @param X A vector of original cross-sectional data
#' @param X1 A vector of half-panel cross-sectional data 1
#' @param X2 A vector of half-panel cross-sectional data 2
#' @param X3 A vector of half-panel cross-sectional data 3
#' @param X4 A vector of half-panel cross-sectional data 4
#'
#' @noRd
#'
hpjecdfest2 <- Vectorize(FUN = function(x, X, X1, X2, X3, X4) {

  # Estimates
  est  <- mean(X  <= x)
  est1 <- mean(X1 <= x)
  est2 <- mean(X2 <= x)
  est3 <- mean(X3 <= x)
  est4 <- mean(X4 <= x)

  # HPJ bias-corrected estimates
  hpjest <- 2 * est - (est1 + est2 + est3 + est4) / 4

  # Ensure valid estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)
  hpjest <- ifelse(hpjest <= 1, hpjest, 1)

  return(hpjest)

}, vectorize.args = "x")

#' Compute bootstrap HPJ empirical CDF estimate for odd T
#'
#' @param quantity An N * 3 matrix of estimates
#' @param indices A vector of indices for bootstrap repetitions
#'
#' @noRd
#'
hpj1_boot <- function(quantity, indices) {

  # Estimates
  est1 <- mean(quantity[indices, 1] <= 0)
  est2 <- mean(quantity[indices, 2] <= 0)
  est3 <- mean(quantity[indices, 3] <= 0)

  # HPJ bias-corrected estimate
  hpjest <- 2 * est1 - (est2 + est3) / 2

  # Ensure valid estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)
  hpjest <- ifelse(hpjest <= 1, hpjest, 1)

  return(hpjest)

}

#' Compute bootstrap HPJ empirical CDF estimate for even T
#'
#' @param quantity An N * 5 matrix of estimates
#' @param indices A vector of indices for bootstrap repetitions
#'
#' @noRd
#'
hpj2_boot <- function(quantity, indices) {

  # Estimates
  est1 <- mean(quantity[indices, 1] <= 0)
  est2 <- mean(quantity[indices, 2] <= 0)
  est3 <- mean(quantity[indices, 3] <= 0)
  est4 <- mean(quantity[indices, 4] <= 0)
  est5 <- mean(quantity[indices, 5] <= 0)

  # HPJ bias-corrected estimate
  hpjest <- 2 * est1 - (est2 + est3 + est4 + est5) / 4

  # Ensure valid estimates
  hpjest <- ifelse(hpjest >= 0, hpjest, 0)
  hpjest <- ifelse(hpjest <= 1, hpjest, 1)

  return(hpjest)

}