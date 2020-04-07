# !!!!!! Turn this into methods for PH and NPH.

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
  prop_haz = TRUE,
  n_events = NA,
  params
){

  post_pred_times = vector("list", length = nrow(params))

  if (prop_haz) {
    for (i in 1:nrow(params)) {
      post_pred_times[[i]] = get_single_conditional_times_PH(
        time,
        event = event,
        trt = trt,
        distribution = distribution,
        n_events = 250,
        params = params[i, ]
      )
    }
  } else {
    for (i in 1:nrow(params)) {

      post_pred_times[[i]] = get_single_conditional_times_NPH(
        time,
        event = event,
        trt = trt,
        distribution = distribution,
        n_events = 250,
        params = params[i, ]
      )
    }
  }

  return(post_pred_times)

}


# get_single_conditional_times = function(x, ...){
#   browser()
#   UseMethod("get_single_conditional_times", x)
# }

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
get_single_conditional_times_PH = function(
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



  return(censor_survival_time(time, additional_time, event, n_events, trt))


}


get_single_conditional_times_NPH = function(
  time,
  event,
  trt,
  distribution,
  n_events = NA,
  params
){

  ctrl_patients <- which(trt == 0)
  trt_patients <- which(trt == 1)

  ctrl_param_index = which(endsWith(names(params), "[1]"))
  trt_param_index = which(endsWith(names(params), "[2]"))

  param_names = unique(remove_bracket(names(params)))

  # Seperate Vectors of random U[0,1] to simulate from both arms seperately
  u_ctrl <- structure(runif(length(ctrl_patients), 0 , 1),
                      class = c(distribution, "numeric"))
  u_trt <- structure(runif(length(trt_patients), 0 , 1),
                     class = c(distribution, "numeric"))

  additional_time_ctrl <- (1 - event[ctrl_patients]) * cond_sample(
    u_ctrl,
    t0 = time[ctrl_patients],
    trt = trt[ctrl_patients],
    params = setNames(params[, ctrl_param_index], param_names)
  )

  additional_time_trt <- (1 - event[trt_patients]) * cond_sample(
    u_ctrl,
    t0 = time[trt_patients],
    trt = trt[trt_patients],
    params = setNames(params[, trt_param_index], param_names)
  )

  additional_time = vector(length = length(trt))
  additional_time[ctrl_patients] = additional_time_ctrl
  additional_time[trt_patients] = additional_time_trt


  return(censor_survival_time(time, additional_time, event, n_events, trt))


}
