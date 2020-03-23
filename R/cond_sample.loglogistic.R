#'  Simulate random variables from a log logistic distribution. The
#'  function is parameterized with hazard: \eqn{\frac{(\beta / \alpha)(t / \alpha)^{\beta-1}}{1+(t / \alpha)^{\beta}}}
#'
#' @inheritParams cond_sample
#' @param params a named list containing \code{alpha}, the scale
#'  parameter (and media), and \code{beta}, the shape parameter
#'
#' @return a vector of failure times conditional on \code{t0}
#' @export
#'
#' @examples
cond_sample.loglogistic <- function(
  u,
  t0,
  HR,
  params
){

  with(
    params,
    alpha * (
      (1 + (t0/alpha) ^ beta) /
        ((1 - u)^ (1 / HR)) - 1) ^ (1 / beta) - t0
  )
}
