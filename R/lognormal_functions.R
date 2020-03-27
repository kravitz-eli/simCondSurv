#' Draws samples from a Lognormal survival model conditional on a timepoint and adjusted for a hazard ratio, \eqn{T \sim f(t|t_0; \gamma)}
#'
#' @inheritParams cond_sample
#' @param params a named list containing \code{mu}, mean of lognormal distribution
#'  on the log-scale, and \code{sigma}, standard deviation of lognormal distribution on the log-scale
#'
#' @return vector of survival times following the Lognormal survival model parameterized by \eqn{mu} and \eqn{sigma}
#' @export
#'
cond_sample.lognormal  <- function(
  u,
  t0,
  trt,
  params
) {

  with(
    params,
    exp(
      sigma * qnorm(
        1 - (1 - u)^(1 / exp(log_HR * trt)) * pnorm(-1 * ( log(t0) - mu ) / sigma )
      ) + mu
    ) - t0
  )


}


get_paramater_names.lognormal = function() c("mu", "sigma")


eval_cdf.lognormal = function(t, trt = 0, params){
  with(
    params,
    plnorm(
      t,
      meanlog = ifelse(trt == 0, log_HR, 0) + mu,
      sdlog = sigma
    )
  )
}

eval_pdf.lognormal = function(t, trt = 0, lambda, nu, log_HR = 0){
  plnorm(
    t,
    meanlog = ifelse(trt == 0, log_HR, 0) + mu,
    sdlog = sigma
  )
}
