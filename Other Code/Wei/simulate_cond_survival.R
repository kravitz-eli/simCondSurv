## Get the posterior predictive distribution of conditional survival times for a 
## SINGLE posterior sample. This function should be interated (for-loop or apply or map)
## over all posterior samples.
simulate_cond_survival = function(
  beta, # a single log hazard posterior sample
  lambda, # a single scale posterior sample
  nu, # a shape log hazard posterior sample
  time, # vector of survival times 
  event, # indicator if survival event has happened, 0 = no event (censored), 1 = event (uncensored)
  trt, #0/1 indicator of control (trt = 0) or treatment (trt = 1)
  add_censor = FALSE, # T/F should simulated conditional survival times be censored?
  n_events = NA # if simulated times will be censored, how many events should be observed before censoring starts
){
  
  # browser()
  
  # Get conditional survival times for one posterior sample
  # Dont get additional times for patients who had an event (event = 1)
  additional_time = (1 - event) * InvCondCDF(nsim = length(time), 
                                             C = time, 
                                             nu = nu,
                                             lambda = lambda, 
                                             beta = exp(beta * trt))
  
  # Append posterior conditional survival times to censored time
  
  # Should simulated data be censored?
  if(add_censor) {
    return(censor_survival_time(time, additional_time, event, n_events, trt))
  }
  else{
    return(data.frame("time" = time + additional_time, 
                      "event" = rep(1, length(time)), # everyone has the event
                      "trt" = trt,
                      "original_time" = time,
                      "additional_time" = additional_time)) 
  }
  
}


