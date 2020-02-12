#' Draws samples from a Lognormal survival model conditional on a timepoint and adjusted for a hazard ratio, \eqn{T \sim f(t|t_0; \gamma)}
#'
#' @param u vector of numbers between 0 and 1. Generally a vector of uniform [0,1] random variable
#' @param t0  time which the survival function is conditioned on, \eqn{S(t|t0)}. When \code{t0 = 0} there is no conditioning. Default to 0
#' @param mu  mean of lognormal distribution on the log-scale, default value of 0
#' @param sigma standard deviation of lognormal distribution on the log-scale, default value of 1
#' @param gamma  hazard ratio, default value of 1
#'
#' @return vector of survival times following the Lognormal survival model parameterized by \eqn{mu} and \eqn{sigma}
#' @export
#'
#' @examples
#' log_normal_inverse_conditional(0.5, t0 = 1, mu = 0.2, sigma = 1, gamma = 1)
#' log_normal_inverse_conditional(seq(0, 1, by = 0.2))
log_normal_inverse_conditional  <- function(
  u,
  t0 = 0,
  mu = 0,
  sigma = 1,
  gamma = 1
) {

  exp(
    sigma * qnorm(
      1 - (1 - u)^(1 / gamma) * pnorm(-1 * ( log(t0) - mu ) / sigma )
    ) + mu
  ) - t0


}
