# Load libraries ------------------------------------------------
library(purrr) # replaces apply() functions with map()
library(survival)
library(rjags)
library(bayestestR)
library(dplyr)
library(ggplot2)

# Load functions -----------------------------------------------
# source("H:/CDK4_6/OS Projection/References/Reference Code/simulate_cond_survival.R") 
# source("H:/CDK4_6/OS Projection/References/Reference Code/censor_survival_time.R") 
# source("H:/CDK4_6/OS Projection/References/Reference Code/inverse_conditional_weibull.R")
# source("H:/CDK4_6/OS Projection/MONALEESA_2/plot_surv_monaleesa2.R")

source(here::here("simulate_cond_survival.R"))
source(here::here("censor_survival_time.R"))
source(here::here("inverse_conditional_weibull.R"))
source(here::here("plot_surv.R"))

# Load example data -------------
# data = read.csv('H:/CDK4_6/OS Projection/MONALEESA_2/Digitization/Derived/monaleesa2_os.csv')
data = read.csv(here::here("monaleesa2_os.csv"))

# KM curves of digitized data
km_fit=survfit(Surv(time, event) ~ trt, data=data)
#par(mar=c(1,1,1,1))  # To resolve figure margins too large issue
plot(km_fit,main="Kaplan-Meier Curve of Digitized Data",xlab="Time (months)",ylab="Probability of Overall Survival",col=c("red","blue"))
legend("bottomleft", legend=c("Ribo","Placebo"),col=c("blue","red"),horiz=FALSE, lty=1, cex=0.8)

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
mod <- jags.model(textConnection(model), data = input, n.chains = 2)
# Burn in
update(mod,1000)
# Posterior Samples
mcmc_out <- coda.samples(
  mod,
  variable.names = c('lambda','nu','beta'),
  n.iter= 1000
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
  upper_bound = 2 * max(data$time) #* 1.5 # How you want to plot the estimated survival function
)

# Posterior Median survival times for each arm  -------------------------------
## median PFS or OS on control arm
control_median_surv = (log(2)/(mcmc_out$lambda))^(1/mcmc_out$nu)
plot(density(control_median_surv, adjust = 1.5), lwd = 3)
# posterior median of median survival times
median(control_median_surv) 
# 95% credible interval
quantile(control_median_surv, c(.025,.975))

##Posterior distribution for true median PFS on Ribo arm
trt_median_surv = (log(2)/(exp(mcmc_out$beta)*mcmc_out$lambda))^(1/mcmc_out$nu)
plot(density(trt_median_surv, adjust = 1.5), lwd = 3)
# Posterior median 
median(trt_median_surv)
quantile(trt_median_surv,c(.025,.975))

# Probability median survival time hits some critical value
critical_value = 50
mean(trt_median_surv > 50)



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
                      n_event = 300)

# You can extract survival times or event indicators from 
# the time_and_event list if you need to
post_survival_times = map_dfc(time_and_event, "time")
post_event = map_dfc(time_and_event, "event")

# Landmark survival -----------------------------------------------------------
# Subset to patients who meet the landmark, then estimate the hazard ratio from
# the subset
landmark_time = 36
# Remove projected survival times that don't meet landmark
landmark_subsets = map(time_and_event, ~subset(.x, time > landmark_time))

# Run your analysis. It could be something like:
# (1) number patients who will hit landmark in each arm
control_patients_landmark = map(landmark_subsets, ~ subset(.x, trt == 0)) %>% 
  map_dbl(nrow)

trt_patients_landmark = map(landmark_subsets, ~ subset(.x, trt == 1)) %>% 
  map_dbl(nrow) 

# Cox regression with subset of patients who hit landmark
cox_landmark = map(landmark_subsets, 
                   ~coxph(Surv(time, event) ~ trt, data =.x)) %>% 
  map_dbl("coefficients")

plot(density(cox_landmark, adjust = 1.5), lwd = 3)

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
#HR_of_interest = 0.65

# Get the posterior hazard ratios and remove any that are less than HR of interest
#posterior_subset = subset(mcmc_out, beta < log(HR_of_interest))

#time_and_event = pmap(.l = posterior_subset, 
#                      .f = simulate_cond_survival, 
#                      time = data$time, 
#                      event = data$event, 
#                      trt = data$trt,
#                      add_censor = TRUE,
#                      n_event = 300)

#post_pred_HR_conditional = map(
# time_and_event, 
#  ~coxph(Surv(time, event) ~ trt, data = .x)
#) %>% 
#  map_dbl("coefficients") %>% 
#  map_dbl(exp)


#plot(
#  density(post_pred_HR_conditional, adjust = 2),
#  lwd = 3,
#  col = "blue",
#  xlab = "Overall Survival Hazard",
#  ylab = "Posterior",
#  xlim = c(-2, 2)
#  )
#par(new=TRUE) # make next plot overlay
plot(
  density(post_pred_HR, adjust = 2),
  lwd = 3,
  col = "black",
  main=c(paste('Posterior Predictive Distribution'),paste('Overall Survival Hazard Ratio at 300 Events')),
  xlab = "Overall Survival Hazard",
  ylab = "Posterior",
  xlim = c(-2, 2)
  #  ann=FALSE, 
  #  axes=FALSE
)

#par(new=FALSE)

# Expected HR and credible intervals
mean(post_pred_HR)
quantile(post_pred_HR, c(.025,.975))

# Expected median survival based on posterior predictive distribution
mean(ctrl_post_pred_median,na.rm = TRUE)
quantile(ctrl_post_pred_median, c(.025,.975), na.rm = TRUE)
mean(trt_post_pred_median, na.rm = TRUE)
quantile(trt_post_pred_median, c(.025,.975), na.rm = TRUE)

# Expected landmark rates
mean(control_patients_landmark/334)
quantile(control_patients_landmark/334, c(.025,.975))
mean(trt_patients_landmark/334)
quantile(trt_patients_landmark/334, c(.025,.975))


# Expected Time required to observe  300 in your example ------------
# 35.3 months at interim 1 with 116 events
num_total_events = 300
time_at_interim = 30 #put the time the interim was conducted 

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

# Posterior Predictive Median Survial Time
control_median_surv = (log(2)/(mcmc_out$lambda))^(1/mcmc_out$nu)
plot(density(control_median_surv, adjust = 1.5), lwd = 3)
# posterior median of median survival times
median(control_median_surv) 
# 95% credible interval
quantile(control_median_surv, c(.025,.975))

## Posterior predictive distribution for median survival time
# Posterior predictive median survival time 
kaplan_meier_no_cnsr = map(
  time_no_censor,
  ~survfit(Surv(time, event) ~trt, data = .x)
)

ctrl_post_pred_median = map_dbl(
  kaplan_meier_no_cnsr,
  ~quantile(.x,.5)$quantile[1]
)

median(ctrl_post_pred_median)

trt_post_pred_median = map_dbl(kaplan_meier_no_cnsr,
                               ~quantile(.x,.5)$quantile[2])

median(trt_post_pred_median)
