log_normal_inverse_conditional  <- function(
  n_sim, # number of random variables to draw
  t0, # conditioning time, NOT on log scale
  mu, # mean
  sigma, # std deviation
  gamma = 1 # hazard ratio
) {

  u = runif(n_sim, 0 , 1)

  exp(
    sigma * qnorm(
      1 - (1 - u)^(1 / gamma) * pnorm(-1 * ( log(t0) - mu ) / sigma )
    ) + mu
  ) - t0


}
