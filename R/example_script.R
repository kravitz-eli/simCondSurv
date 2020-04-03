library(rjags)
library(magrittr)
library(survival)
library(survminer)

source("R/utils-pipe.R")
source("R/cond_sample.R")
source("R/cond_sample.weibull.R")

source("R/run_jags.R")
source("R/get_conditional_times.R")
source("R/censor_survival_time.R")
source("R/run_projections.R")

data("monaleesa2")

time = monaleesa2$time
event = monaleesa2$event
trt = monaleesa2$trt



distribution = "loglogistic"

# Fit the JAGS model ---------------------------------------
n_iter = 1e3
n_burn = 1e3
n_adapt = 1e3

model = run_jags(
  distribution = distribution,
  prop_haz = TRUE,
  time = time,
  event = event,
  trt = trt,
  chains = 1,
  n_iter = n_iter,
  n_burn = n_burn,
  n_adapt = n_adapt
)

post_params = tibble::as_tibble(model[[1]])
# names(post_params)[1] = "log_HR"

colMeans(post_params)


surv_times = get_conditional_times(
  time = time,
  event = event,
  trt = trt,
  distribution = distribution,
  n_events = 250,
  params = post_params
)

time_projections = get_projections(
  simulated_times = surv_times,
  crit_values = 0.75
)

plots = make_plots(
  params = post_params,
  data = monaleesa2,
  post_pred_HR = time_projections$HR_at_events,
  axis_min_time = 0,
  axis_max_time = 100
)

plots$survival_plot


