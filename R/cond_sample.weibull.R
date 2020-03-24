#'  Simulate random variables from a Weibull distribution. The
#'  function is parameterized with hazard: \eqn{\nu \lambda t^{nu-1}}
#'
#' @inheritParams cond_sample
#' @param params a named list containing \code{lambda}, the Weibull scale
#'  parameter, and \code{nu}, the Weibull shape parameter
#'
#' @return a vector of failure times conditional on \code{t0}
#' @export
#'
#' @examples
cond_sample.weibull <- function(
  u,
  t0,
  HR,
  params
){

  with(
    params,
    (-log(1 - u) / (lambda * HR) + t0^nu) ^ (1/nu) - t0
  )
}
