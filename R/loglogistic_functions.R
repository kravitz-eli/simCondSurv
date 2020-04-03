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
  trt,
  params
){

  with(
    params,
    alpha * (
      (1 + (t0/alpha) ^ beta) /
        ((1 - u)^ (1 / exp(log_HR * trt))) - 1) ^ (1 / beta) - t0
  )
}

eval_cdf.loglogistic = function(t, trt = 0, params){

  with(
    params,
    plog_logistic(t, alpha = alpha, beta = beta)^exp(log_HR * trt)
  )
}

eval_pdf.loglogistic = function(t, trt = 0, params){
  dlnorm(
    t,
    meanlog = ifelse(trt == 0, log_HR, 0) + mu,
    sdlog = sigma
  )
}

plog_logistic = function(x, alpha, beta){
   x^beta / (alpha^beta + x^beta)

}

dlog_logistic = function(x, alpha, beta){

  ((beta / alpha) * (x / alpha) ^ (beta - 1)) / (alpha + (x / alpha) ^ beta)^2
}
