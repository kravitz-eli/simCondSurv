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
  data_and_model,
  n_events = NA,
){

  post_pred_times = vector("list", length = nrow(params))


  for (i in 1:nrow(params)) {
    post_pred_times[[i]] = get_single_conditional_times_PH(
      data_and_model,
      post_params_i = post_params,
      n_events = 250
    )
  }

  # if (prop_haz) {
  #   for (i in 1:nrow(params)) {
  #     post_pred_times[[i]] = get_single_conditional_times_PH(
  #       time,
  #       event = event,
  #       trt = trt,
  #       distribution = distribution,
  #       n_events = 250,
  #       params = params[i, ]
  #     )
  #   }
  # } else {
  #   for (i in 1:nrow(params)) {
  #
  #     post_pred_times[[i]] = get_single_conditional_times_NPH(
  #       time,
  #       event = event,
  #       trt = trt,
  #       distribution = distribution,
  #       n_events = 250,
  #       params = params[i, ]
  #     )
  #   }
  # }

  return(post_pred_times)

}


get_single_conditional_times = function(x, ...){
  UseMethod("get_single_conditional_times", x)
}

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
get_single_conditional_times.PH = function(
  data_and_model,
  post_params_i, # a single draw from the posterior distribution of each parameter. 1xp data.frame
  n_events = NA
){

  u <- structure(runif(length(time), 0 , 1), class = c(data_and_model$distribution, "numeric"))
  additional_time <- (1 - event) * cond_sample(u,
                                               t0 = data_and_model$time,
                                               trt = data_and_model$trt,
                                               params = post_params_i)



  return(censor_survival_time(time, additional_time, event, n_events, trt))


}


get_single_conditional_times.NPH = function(
  data_and_model,
  post_params_i, # a single draw from the posterior distribution of each parameter. 1xp data.frame
  n_events = NA
){



  # Seperately simulate survival times for two arms. Can't simulate at the same
  # time when there is no hazard ratio               ---------------------------

  # control arm: ----------------------------------
  ctrl_patients <- which(data_and_model$trt == 0)
  u_ctrl <- structure(runif(length(ctrl_patients), 0 , 1),
                      class = c(distribution, "numeric"))

  additional_time_ctrl <- (1 - data$event[ctrl_patients]) * cond_sample(
    u_ctrl,
    t0 = data$time[ctrl_patients],
    trt = data$trt[ctrl_patients],
    params = setNames(post_params_i[, data_and_model$ctrl_param_index],
                      data_and_model$param_names)
  )

  # Treatment arm --------------------------------------------
  u_trt <- structure(runif(length(trt_patients), 0 , 1),
                     class = c(distribution, "numeric"))
  trt_patients <- which(data_and_model$trt == 1)

  additional_time_trt <- (1 - data$event[trt_patients]) * cond_sample(
    u_trt,
    t0 = data$time[trt_patients],
    trt = data$trt[trt_patients],
    params = setNames(data_and_model$post_params[, data_and_model$trt_param_index],
                      data_and_model$param_names)
  )

  # Combine simulated survival times from both arms ------------------------
  additional_time = vector(length = length(data$trt))
  additional_time[ctrl_patients] = additional_time_ctrl
  additional_time[trt_patients] = additional_time_trt


  return(censor_survival_time(time, additional_time, event, n_events, trt))


}
