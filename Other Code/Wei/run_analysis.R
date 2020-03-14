# Load libraries ------------------------------------------------
library(purrr) # replaces apply() functions with map()
library(survival)
library(rjags)

# Load functions -----------------------------------------------
source("simulate_cond_survival.R") 
source("censor_survival_time.R") 
source("inverse_conditional_weibull.R")
source("plot_surv.R")

# Load example data -------------
data = readRDS("example_weibull.rds")

## Weibull jags model ----------------------------------------------------------
model = '
model{

# Likelihood 
  for(i in 1:N){
    # Bernoulli random variable to add censoring
    ones[i] ~ dbern(p[i] / 100000)

    ## ith contribution to the joint likelihood, changes if there is censoring
    p[i] <- exp(-H[i] )* pow(h[i], event[i])

    ## hazard and cumulative hazard functions
    h[i] <- exp(beta* TRT[i] )*nu*lambda*pow(t[i],nu-1)
    H[i] <- exp(beta* TRT[i] )*lambda*pow(t[i],nu)
  }

#Priors 
  beta ~ dnorm(0,.001) #hazard ratio
  nu~dnorm(0,.1)T(0,) #shape parameter
  lambda~dgamma(.1,.1) #scale parameter

}'

# Run JAGS ------------------------------------------------------
input <- list(
  t = data$time, #survival times
  event = data$event, # event indicator
  N = nrow(data), #total number of patients
  TRT =  data$trt, # 0/1 indicator of control or treatment
  ones = rep(1, times = nrow(data))  # don't remove, used to add censored datapoints to likelihood
)

# Declare model
mod <- jags.model(textConnection(model),data = input, n.chains = 1)
# Burn in
update(mod, 1000)
# Posterior Samples
mcmc_out <- coda.samples(
  mod,
  variable.names = c('lambda','nu','beta'),
  n.iter = 1000
)

mcmc_out = as.data.frame(as.matrix(mcmc_out))

# saveRDS(mcmc_out, "example_MCMC.RDS")

# Plot -----------------------------------------------------------------------
# Plot Kaplan Meier curve with estimate survival function overlayed
plot_surv(
  beta = mcmc_out$beta,
  lambda = mcmc_out$lambda,
  nu = mcmc_out$nu, 
  data = data,
  upper_bound = max(data$time) * 2 # How you want to plot the estimated survival function
)

# Posterior Median survival times for each arm  -------------------------------
## median PFS or OS on control arm
control_median_surv = (log(2)/(mcmc_out$lambda))^(1/mcmc_out$nu)
plot(density(control_median_surv, adjust = 2), lwd = 3)
# posterior median of median survival times
median(control_median_surv) 
# 95% credible interval
quantile(control_median_surv, c(.025,.975))

##Posterior distribution for true median PFS on Ram arm
trt_median_surv = (log(2)/ (exp(mcmc_out$beta)*mcmc_out$lambda) )^(1/mcmc_out$nu)
plot(density(trt_median_surv, adjust = 2), lwd = 3)
# Posterior median 
median(trt_median_surv)
quantile(trt_median_surv,c(.025,.975))

# Probability median survival time hits some critical value
critical_value = 6
mean(trt_median_surv > 6)



# Calculate additional times for patients who are censored 
# from posterior predictive distribution ---------------------------------------

# Get a list with data.frame as every element with survival times and event indicato
# Simulated survival times are censored after a certain number of events
# number of events in controled by argument: n_event
time_and_event = pmap(.l = mcmc_out, 
                      .f = simulate_cond_survival, 
                      time = data$time, 
                      event = data$event, 
                      trt = data$trt,
                      add_censor = TRUE,
                      n_event = 294)


# Landmark survival -----------------------------------------------------------
# Subset to patients who meet the landmark, then estimate the hazard ratio from
# the subset
landmark_time = 13
# Remove projected survival times that don't meet landmark
landmark_subsets = map(time_and_event, ~subset(.x, time > landmark_time))

# Run your analysis. It could be something like:
# (1) number patients who will hit landmark in each arm
control_patients_landmark = map(landmark_subsets, ~ subset(.x, trt == 0)) %>% 
  map_dbl(nrow) %>% 
  mean() 

trt_patients_landmark = map(landmark_subsets, ~ subset(.x, trt == 1)) %>% 
  map_dbl(nrow) %>% 
  mean()

# Cox regression with subset of patients who hit landmark
cox_landmark = map(landmark_subsets, 
                   ~coxph(Surv(time, event) ~ trt, data =.x)) %>% 
  map_dbl("coefficients") %>% 
  exp()

mean(cox_landmark)

plot(density(cox_landmark, adjust = 2), lwd = 3)

# Posterior predictive distribution  ------------------------------------------
# Conditional on observed hazard ratio
post_pred_HR = map(
  time_and_event, 
  ~coxph(Surv(time, event) ~ trt, data = .x)
) %>% 
  map_dbl("coefficients") %>% 
  map_dbl(exp)

# Posterior predictive median survival time 
kaplan_meier_fits = map(
  time_and_event,
  ~survfit(Surv(time, event) ~trt, data = .x)
)

ctrl_post_pred_median = map_dbl(kaplan_meier_fits,
                                ~quantile(.x,.5)$quantile[1])

trt_post_pred_median = map_dbl(kaplan_meier_fits,
                               ~quantile(.x,.5)$quantile[2])


# Conditional on a HR of interest (like 0.737 in the slides) ------------------
HR_of_interest = 0.45

# Get the posterior hazard ratios and remove any that are less than HR of interest
posterior_subset = subset(mcmc_out, beta < log(HR_of_interest))

time_and_event_conditional = pmap(.l = posterior_subset, 
                                  .f = simulate_cond_survival, 
                                  time = data$time, 
                                  event = data$event, 
                                  trt = data$trt,
                                  add_censor = TRUE,
                                  n_event = 190)

post_pred_HR_conditional = map(
  time_and_event_conditional, 
  ~coxph(Surv(time, event) ~ trt, data = .x)
) %>% 
  map_dbl("coefficients") %>% 
  map_dbl(exp)


plot(
  density(post_pred_HR_conditional, adjust = 2),
  lwd = 3,
  col = "blue",
  xlab = "Overall Survival Hazard",
  ylab = "Posterior",
  xlim = c(0, 1)
)
par(new=TRUE) # make next plot overlay
plot(
  density(post_pred_HR, adjust = 2),
  lwd = 3,
  xlim = c(0, 1),
  ann=FALSE, 
  axes=FALSE
)
par(new=FALSE) #turn plots back to normal


# Get average time to required number of event, 300 in your example ------------
num_total_events = 150
time_at_interim = 3 #put the time the interim was conducted 
  
# Generate conditional survival times with no censoring
time_no_censor = pmap(.l = mcmc_out, 
                      .f = simulate_cond_survival, 
                      time = data$time, 
                      event = data$event, 
                      trt = data$trt,
                      add_censor = FALSE)

# Sort the survival times and find the time when the required number of events happen
time_to_n_events = time_no_censor %>% 
  map_dbl(~sort(.x$time)[num_total_events])

time_from_interim = time_to_n_events - time_at_interim

median(time_from_interim)
plot(density(time_from_interim, adjust = 1.5), 
     main = "Time from Interim to Reach Endpoint",
     xlab = "Time",
     lwd = 2)

