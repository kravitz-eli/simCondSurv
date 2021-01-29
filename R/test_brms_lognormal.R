# Libraries ------------------------------------------------
library(tidyverse)
library(survival)
library(flexsurv)
library(brms)

# Function to simulation lognormal survival times -----------------------------
sim_lnorm_surv <- function(n_pts, meanlog, sdlog, hr, n_events) {

  tibble(
    id = seq(1, n_pts),
    arm = sample(c(0, 1), n_pts, replace = TRUE),
    true_time = rlnorm(n_pts, meanlog = meanlog + log(hr) * arm, sdlog = sdlog)
  ) %>%
    mutate(
      cens_time = nth(true_time, n_events, order_by = true_time),
      event = as.numeric(true_time < cens_time),
      obs_time = ifelse(event, true_time, cens_time)
    )


}

# Try again with lognormal --------------------------------------------
n_pts <- 1e3
meanlog = log(20) # approx 3
sdlog <- 0.25
hr <- 0.65
n_events <- n_pts * 4/5

pts_df <- sim_lnorm_surv(n_pts, meanlog, sdlog, hr, n_events)
freq_model <- flexsurvreg(Surv(obs_time, event) ~ arm, dist = "lnorm", data = pts_df)


bayes_model <- brm(obs_time | cens(event) ~ arm,
                   family = "lognormal",
                   data = pts_df)
