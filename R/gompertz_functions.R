#'  Simulate random variables from a Gompertz distribution. The
#'  function is parameterized with hazard: \eqn{\lambda(x) =  a \exp(bx)}
#'
#' @inheritParams cond_sample
#' @param params list with elements a and b, the parameters of a gompertz distribution.
#'
#' @return
#' @export
cond_sample.gompertz <- function(
  u,
  t0,
  trt = 0,
  params
){

  with(
    params,
    1 / b * log(1 - b/(a * exp(log_HR * trt)) * log(1 - u) * exp(-b * t0))
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
