censor_survival_time = function(time, additional_time, event, n_events, trt){
  
  censored_subjects <- which(event == 0)
  num_observed_events <- sum(event)
  n_additional_events <- n_events - num_observed_events
  
  censor_time <- sort(additional_time[censored_subjects])[n_additional_events]
  
  event[censored_subjects] <- as.numeric(additional_time[censored_subjects] <= censor_time)
  # total_time[censored_subjects] <- pmin(additional_time[censored_subjects], censor_time) 
  additional_time <- pmin(additional_time, censor_time) 
  
  
  return(data.frame("time" = time + additional_time, 
                    "event" = event, 
                    "trt" = trt,
                    "original_time" = time,
                    "additional_time" = additional_time))
  
}