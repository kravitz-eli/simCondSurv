#' Draws samples from a Lognormal survival model conditional on a timepoint and adjusted for a hazard ratio, \eqn{T \sim f(t|t_0; \gamma)}
#'
#' @inheritParams cond_sample
#' @param params a named list containing \code{mu}, mean of lognormal distribution
#'  on the log-scale, and \code{sigma}, standard deviation of lognormal distribution on the log-scale
#'
#' @return vector of survival times following the Lognormal survival model parameterized by \eqn{mu} and \eqn{sigma}
#' @export
#'
#' @examples
#' log_normal_inverse_conditional(0.5, t0 = 1, mu = 0.2, sigma = 1, gamma = 1)
#' log_normal_inverse_conditional(seq(0, 1, by = 0.2))
cond_sample.lognormal  <- function(
  u,
  t0,
  HR,
  params
) {

  with(
    params,
    exp(
      sigma * qnorm(
        1 - (1 - u)^(1 / HR) * pnorm(-1 * ( log(t0) - mu ) / sigma )
      ) + mu
    ) - t0
  )


}
