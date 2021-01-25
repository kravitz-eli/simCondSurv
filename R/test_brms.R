## Compare BRMS parametric survival modeling to JAGS model
## Start with Weibull regression and see how other things work.

# Load libraries ----------------------------------
library(tidyverse)
library(brms)
library(survival)
library(flexsurv)

# Functions -------------------------------------------------------------
source(here::here("R", "simulate_weibull.R"))

# Simulate data to test brms -----------------------------------------------
n_pts <- 1e3
n_events = n_pts * 0.95
# Proportional hazards Weibull model (doesn't work)
set.seed(3)

pts <- simulate_weibull(N = n_pts, lambda = 0.31, rho = -0.4, beta = log(0.80), n_events = n_events)
hist(pts$time)
pts %>%
  group_by(x) %>%
  summarise(sum(status))

# Can't simulate good survival times. Just use ovarian dataset
freq_model <- survreg(Surv(futime, fustat) ~ rx, ovarian, dist='weibull')



# scale = 1 / shape
1 / freq_model$scale
# Scale = exp(intercept)
exp(freq_model$coefficients[["(Intercept)"]])



# Can't simulate good survival times. Just use ovarian dataset
freq_model <- survreg(Surv(time, status) ~ x, data = pts, dist='weibull')

freq_model <- flexsurvreg(Surv(time, status) ~ x, data = pts, dist = 'weibull')
freq_model$coefficients
# scale = 1 / shape
1 / freq_model$coefficients[["scale"]]
# Scale = exp(intercept)
freq_model$coefficients
# HR = exp(x)
freq_model$coefficients[["x"]] %>% exp()

# Try with BRMS -------------------------------------------------------------
bayes_model <- brm(time | cens(status)  ~ x,
                   data = pts,
                   family = weibull(),
                   save_model = "weibull_ph")
