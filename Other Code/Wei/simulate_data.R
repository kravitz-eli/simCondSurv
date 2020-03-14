simulWeib <- function(N, lambda, rho, beta, rateC, n_event)
{
  # covariate --> N Bernoulli trials
  x <- sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5))
  
  # Weibull latent event times
  v <- runif(n=N)
  surv_times <- (- log(v) / (lambda * exp(x * beta)))^(1 / rho)
  
  # censoring time
  C = surv_times[order(surv_times)][n_event]
  
  # follow-up times and event indicators
  time <- pmin(surv_times, C)
  status <- as.numeric(surv_times <= C)
  
  # data set
  data.frame(id = 1:N,
             time = time,
             event = status,
             trt = x)
}

set.seed(3)
data = simulWeib(N= 200, lambda = log(2)/6, rho = 1.2, beta=-0.6, n_event = 140)
saveRDS(data, "example_weibull.Rds")
