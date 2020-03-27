#'  Simulate random variables from a Weibull distribution. The
#'  function is parameterized with hazard: \eqn{\nu \lambda t^{nu-1}}
#'
#' @inheritParams cond_sample
#' @param params a named list containing \code{lambda}, the Weibull scale
#'  parameter, and \code{nu}, the Weibull shape parameter, and \code{log_HR},
#'  log hazard ratio
#'
#' @return a vector of failure times conditional on \code{t0}
#' @export
#'
#' @examples
cond_sample.weibull <- function(
  u,
  t0,
  trt = 0,
  params
){

  with(
    params,
    (-log(1 - u) / (lambda * exp(log_HR * trt)) + t0^nu) ^ (1/nu) - t0
  )
}


eval_cdf.weibull = function(t, trt = 0, params){
  with(
    params,
    pweibull(
      t,
      shape = nu,
      scale = (exp(log_HR * trt) * lambda) ^ (-1 / nu)
    )
  )
}

eval_pdf.weibull = function(t, trt = 0, lambda, nu, log_HR = 0){
  dweibull(
    t,
    shape = nu,
    scale = (exp(log_HR * trt) * lambda) ^ (-1 / nu)
  )
}
