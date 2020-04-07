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
  log_HR = 0,
  params
){

  with(
    params,
    (-log(1 - u) / (lambda * exp(log_HR * trt)) + t0^nu) ^ (1/nu) - t0
  )
}


#' Alternate parameterization of Weibull CDF with hazard: \eqn{\nu \lambda t^{nu-1}}
#'
#' Probability density function of Weibull distribution with the proportional
#' hazards parameterzation. Wrapper to pweibull with arguments that match
#' the survival parameterrization
#'
#' @param x vector of quantiles
#' @param lambda scale parameter. This is multiplied by hazard ratio in
#' proportional hazards model
#' @param nu shape parameter
#'
#' @return Weibull CDF evaluated at \code{x}
#' @export
#'
#' @examples
pweibull_surv = function(x, lambda, nu) {
  pweibull(
    t,
    shape = nu,
    scale = lambda ^ (-1 / nu)
  )
}

#' Alternate parameterization of Weibull PDF with hazard: \eqn{\nu \lambda t^{nu-1}}
#'
#' Probability density function of Weibull distribution with the proportional
#' hazards parameterzation. Wrapper to pweibull with arguments that match
#' the survival parameterrization
#'
#' @param x vector of quantiles
#' @param lambda scale parameter. This is multiplied by hazard ratio in
#' proportional hazards model
#' @param nu shape parameter
#'
#' @return Weibull PDF evaluated at \code{x}
#' @export
#'
#' @examples
dweibull_surv = function(x, lambda, nu) {
  dweibull(
    t,
    shape = nu,
    scale = lambda ^ (-1 / nu)
  )
}


eval_cdf.weibull = function(t, trt = 0, log_HR = 0, params){
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
