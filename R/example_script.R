library(rjags)
library(magrittr)

source("R/utils-pipe.R")
source("R/cond_sample.R")
source("R/cond_sample.weibull.R")

source("R/run_jags.R")
source("R/get_conditional_times.R")
source("R/censor_survival_time.R")

monaL = data("monaleesa2")

time = monaleesa2$time
event = monaleesa2$event
trt = monaleesa2$trt



distribution = "weibull"
# params = list("lambda" = 1.5, "nu" = 0.20)

# n_events = 250
# log_HR = 0

# Fit the JAGS model ----
n_iter = 5e3
n_burn = 1e4
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
  n_adapt = n_adapt,
  track_variable_names = c("beta", "nu", "lambda")
)

post_params = tibble::as_tibble(model[[1]]) %>%
  setNames(c("log_HR", "lambda", "nu"))
colMeans(post_params)


surv_times = get_conditional_times(
  time = time,
  event = event,
  trt = trt,
  distribution = distribution,
  n_events = 250,
  params = post_params
)
