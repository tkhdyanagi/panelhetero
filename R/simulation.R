#' Generate artificial data
#'
#' The `simulation()` function enables to generate artificial data from an
#' AR(1) model with random coefficients.
#' The function is used in the package vignette.
#'
#' @param N The number of cross-sectional units
#' @param S The length of time series
#'
#' @returns An N times S matrix of panel data
#'
#' @examples
#' panelhetero::simulation(N = 300, S = 50)
#'
#' @export
#'
simulation <- function(N, S) {

  # Heterogeneous random coefficients
  c     <- stats::rnorm(N, mean = 0, sd = 1)
  phi   <- stats::rbeta(N, shape1 = 2, shape2 = 5)
  sigma <- stats::rbeta(N, shape1 = 2, shape2 = 5) + 0.5

  # Construct the N times (S + 1) matrix of panel data
  y <- matrix(0, nrow = N, ncol = S + 1)

  # Generate panel data from AR(1) model with random coefficients
  for (i in 1:N) {

    # Unobserved initial value
    y[i, 1] <- stats::rnorm(1, mean = c[i], sd = sigma[i])

    for (t in 2:(S + 1)) {
      y[i, t] <- (1 - phi[i]) * c[i] +
        phi[i] * y[i, t-1] +
        sqrt((1 - phi[i])^2 * sigma[i]^2) * stats::rnorm(1, mean = 0, sd = 1)
    }

  }

  # The N * S matrix of panel data
  y <- y[, -1]

  return(y)

}
