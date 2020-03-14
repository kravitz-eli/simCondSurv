# Expected Time required to observe 300 events  -----------------------
new_events_at_interim = 116
# time in months from beginning of study that previous interim was conducted 
time_at_last_interim = 33.5 

# Subset the data to only patients censored at last interim
censored_patients = filter(data, event == 0)
simulate_censored_times = pmap(.l = mcmc_out, 
                               .f = simulate_cond_survival, 
                               time = censored_patients$time, 
                               event = censored_patients$event, 
                               trt = censored_patients$trt,
                               add_censor = TRUE,
                               n_events = 116)

# Time after interim analysis to reach 300 total events
# First calculate all the additional 
additional_time = map_dbl(
  simulate_censored_times, 
  ~{ .x %>% 
      mutate(
        "original_time" = censored_patients$time, # get the survival times in the original data
        "added_time" = time - original_time # see how much time was added to each subject
      ) %>% 
      filter(event == 1) %>%  # Only keep the patients who had an event
      summarise(max(time)) %>%  # Find the last timepoint an event occured
      unlist()
  })

median(additional_time)
plot(density(additional_time, adjust = 2), 
     main = "Time After Interim to Reach Endpoint",
     xlab = "Months",
     lwd = 2,
     col = "navy")
abline(v = median(additional_time), lty = 2, lwd = 2, col = "grey")



time_to_300 = time_at_last_interim + additional_time

median(time_to_300)
plot(density(time_to_300, adjust = 2), 
     main = "Time from Beginning of Study to Reach Endpoint",
     xlab = "Time",
     lwd = 2,
     col = "navy")
abline(v = median(time_to_300), lty = 2, lwd = 2, col = "grey")



# Posterior Predictive Median Survial Time -------------------------------------
# Method 1: Calculate from parametric weibull model 
# Distribution of median
control_median_surv = (log(2)/(mcmc_out$lambda))^(1/mcmc_out$nu) #formula for median
trt_median_surv = (log(2)/(exp(mcmc_out$beta)*mcmc_out$lambda))^(1/mcmc_out$nu) # formula for median

# Point estimate: median of distribution of median survival times
median(control_median_surv)
median(trt_median_surv)

# Method 2: Nonparametric using Kaplan Meier Estimate. 
# Median OS will not be mature for both arms at 300 events, so we'll generate event times
# without censoring and then estimate the median.

# Generate conditional survival times with no censoring
time_no_censor = pmap(.l = mcmc_out, 
                      .f = simulate_cond_survival, 
                      time = data$time, 
                      event = data$event, 
                      trt = data$trt,
                      add_censor = FALSE)

# Fit a Kaplan Meier estimate to each curve 
kaplan_meier_no_cnsr = map(
  time_no_censor,
  ~survfit(Surv(time, event) ~trt, data = .x)
)

# Get the median survival time from each survival curve in control arm
ctrl_post_pred_median = map_dbl(
  kaplan_meier_no_cnsr,
  ~quantile(.x,.5)$quantile[1]
)

median(ctrl_post_pred_median)

# Get the median survival time from each survival curve in control arm
trt_post_pred_median = map_dbl(
  kaplan_meier_no_cnsr,
  ~quantile(.x,.5)$quantile[2]
)

median(trt_post_pred_median)
