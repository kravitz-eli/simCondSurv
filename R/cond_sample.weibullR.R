#'  Simulate random variables from a Weibull distribution. The
#'  function is parameterized with hazard: \eqn{\nu \lambda t^{nu-1}}
#'
#' @inheritParams cond_sample
#' @param lambda scale parameter
#' @param nu shape parameter
#'
#' @return
#' @export
#'
#' @examples
cond_sample.weibull = function(
  u,
  t0,
  nu,
  lambda,
  HR = 1
){

  (-log(1 - u) / (lambda * HR) + t0^nu) ^ (1/nu) - t0

}
