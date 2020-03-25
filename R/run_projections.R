get_projections = function(simulated_times, crit_values) {

  time_to_events = sapply(simulated_times, function(x) x$censor_time[[1]])

  HR_at_events = vector("numeric", length(simulated_times))
  post_pred_median_ctrl = vector("numeric", length(simulated_times))
  post_pred_median_trt = vector("numeric", length(simulated_times))

  for(i in seq_along(simulated_times)) {

    coxPH_model = survival::coxph(survival::Surv(time, event) ~ trt, data = simulated_times[[i]])
    HR_at_events[[i]] = exp(unname(coxPH_model$coefficients))

    KM_fit = survival::survfit(survival::Surv(time, event) ~trt, data = simulated_times[[i]])

    # Get the median survival times for the control and experiment arm
    post_pred_median_ctrl[[i]] = survival:::quantile.survfit(KM_fit, probs = 0.50)$quantile[1]
    post_pred_median_trt[[i]] = survival:::quantile.survfit(KM_fit, probs = 0.50)$quantile[2]
  }


  prob_reject = sapply(crit_values, function(x) mean(HR_at_event < x))


}
