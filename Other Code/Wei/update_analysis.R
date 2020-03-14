# Median survival times projections  -------------------------------
# Calculate median survival time from Weibull model
# Formula for median:  (log(2) / lambda) ^ (1 / nu)

# #Posterior distribution for true median OS on placebo + letrozole
control_median_surv = (log(2)/(mcmc_out$lambda))^(1/mcmc_out$nu)
median(control_median_surv)
# 95% credible interval
quantile(control_median_surv, c(.025,.975))

plot(density(control_median_surv, adjust = 1.5), lwd = 3)

# #Posterior distribution for true median OS on Ribo arm
trt_median_surv = (log(2)/(exp(mcmc_out$beta)* mcmc_out$lambda))^(1/mcmc_out$nu)
plot(density(trt_median_surv, adjust = 1.5), lwd = 3)
# Posterior median 
median(trt_median_surv)
quantile(trt_median_surv,c(.025,.975))

# Probability median survival time hits some critical value
critical_value = 50
mean(trt_median_surv > 50)


# Probability median OS data is mature at interim analysis -------------------
kaplan_meier_fits = map(
  time_and_event,
  ~survfit(Surv(time, event) ~ trt, data = .x)
)

# Calculate median OS for control arm.
# Will be NA if >50% patients are projected to be alive
ctrl_median = map_dbl(
  kaplan_meier_fits,
  ~quantile(.x,.5)$quantile[1]
)

prob_ctrl_mature = mean(ifelse(is.na(ctrl_median), yes = 0, no = 1))

# Median if mature
ctrl_median %>% 
  median(na.rm = TRUE)

## Calculate median OS for Ribo arm
# Will be NA if >50% patients are projected to be alive
trt_median = map_dbl(
  kaplan_meier_fits,
  ~quantile(.x,.5)$quantile[2]
)

prob_ctrl_mature = mean(ifelse(is.na(trt_median), yes = 0, no = 1))

# Median if mature, I don't think this is a useful estimate of median survival
trt_median %>% 
  median(na.rm = TRUE)


# Expected Time required to observe 300 events  -------------------------------

# Get the time to 300 events for posterior OS simulation This be the censoring time
# for each simulation

# addit_time = time_and_event %>% 
#   map_dbl(~max(.x$time))
# 
# summary(addit_time)
# 
# censored_patients = data %>%
#   filter(event == 0)
# 
# cens_time_1 = data %>% 
#   filter(event == 1) %>%
#   summarise(max(time)) %>% 
#   unlist()
# 
# 
# cens_time_2 = time_and_event %>% 
#   map_dbl(~max(.x$time))
# 
# time_at_last_interim = 33.5 
# 
# time_from_interim = cens_time_2 - cens_time_1
# 
# interim_time = time_at_last_interim + time_from_interim

# time_censored = pmap(.l = mcmc_out,
#      .f = simulate_cond_survival,
#      time = censored_patients$time,
#      event = censored_patients$event,
#      trt = censored_patients$trt,
#      add_censor = TRUE,
#      n_event = 184)

# additional_time = map_dbl(time_censored, ~max(.x$additional_time))

# additional_time = map(time_and_event, ~filter(.x, event == 1)) %>%
  # map_dbl(~max(.x$time))

additional_time = map_dbl(time_and_event, ~max(.x$additional_time))
time_at_first_censor = map_dbl(time_and_event, ~max(.x$original_time))  

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
