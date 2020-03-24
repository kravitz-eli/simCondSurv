library(rjags)

source("R/cond_sample.R")
source("R/cond_sample.weibull.R")

source("R/simulate_cond_survival.R")
source("R/censor_survival_time.R")

monaL = data("monaleesa2")

time = monaleesa2$time
event = monaleesa2$event
trt = monaleesa2$trt



distribution = "weibull"
params = list("lambda" = 1.5, "nu" = 0.20)

n_events = 250
log_HR = 0


model = run_jags(
  distribution = distribution,
  prop_haz = TRUE,
  time = time,
  event = event,
  trt = trt,
  chains = 1,
  iter = 1000,
  burn = 1000,
  n_adapt = 1000,
  track_variable_names = c("beta", "nu", "lambda")
)

post_params = as_tibble(model[[1]])

surv_times = get_conditional_times(
  time = time,
  event = event,
  trt = trt,
  distribution = distribution,
  n_events = 250,
  log_HR = tibble::deframe(post_params[1, "beta"]),
  params = as.list(post_params[1, c("lambda", "nu")])
)


