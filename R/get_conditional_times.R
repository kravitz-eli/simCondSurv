
#' Generate ONE sample from the posterior predictive of conditional survival times.
#'
#' Get the posterior predictive distribution of conditional survival times for a
#' SINGLE posterior sample. This function should be interated
#' over all posterior samples.
#'
#' @param time positive numeric. Survival times
#' @param event logical or 0/1. \code{event = 1} indicates failure,
#' \code{event = 0} indicates censoring
#' @param trt 0/1 vector. Indicates assignment to treatment or control
#' @param log_HR log hazard ratio, usually posterior sample
#' @param n_events numeric or NULL. How many events should be observed before
#' censoring. If NA, there will be no censoring time.
#'
#' @return
#' @export
#'
#' @examples
get_single_conditional_times = function(
  time,
  event,
  trt,
  distribution,
  n_events = NA,
  params
){

  u <- structure(runif(length(time), 0 , 1), class = c(distribution, "numeric"))
  additional_time <- (1 - event) * cond_sample(u,
                                               t0 = time,
                                               trt = trt,
                                               params = params)


  # Simulated data is censored
  if(!is.na(n_events)) {
    return(censor_survival_time(time, additional_time, event, n_events, trt))
  }
  else{   # Simulated data is NOT censored
    return(data.frame("time" = time + additional_time,
                      "event" = 1,
                      "trt" = trt,
                      "original_time" = time,
                      "additional_time" = additional_time,
                      "censor_time" = censor_time))
  }

}

#' Generate samples from the posterior predictive of conditional survival times.
#'
#'
#` @inheritParams get_single_conditional_times
#'
#' @return
#' @export
#'
#' @examples
get_conditional_times = function(
  time,
  event,
  trt,
  distribution,
  n_events = NA,
  params
){

  post_pred_times = vector("list", length = nrow(params))

  for (i in 1:nrow(params)) {
    post_pred_times[[i]] = get_single_conditional_times(
      time = time,
      event = event,
      trt = trt,
      distribution = distribution,
      n_events = 250,
      params = params[i, ]
    )
  }

  return(post_pred_times)

}
