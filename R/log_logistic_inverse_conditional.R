#'  Simulate random variables from a Weibull distribution. The
#'  function is parameterized with hazard: \eqn{\nu \lambda t^{nu-1}}
#'
#' @inheritParams cond_sample
#' @param params a named list containing \code{alpha}, the scale
#'  parameter (and media), and \code{beta}, the shape parameter
#'
#' @return a vector of failure times conditional on \code{t0}
#' @export
#'
#' @examples
inverse_log_logistic = function(
  u,
  t0,
  HR = 1,
  params
){

  with(
    params,
    alpha * (
      (1 + (t0/alpha) ^ beta) /
        ((1 - u)^ (1 / HR)) - 1) ^ (1 / beta) - t0
  )
}
