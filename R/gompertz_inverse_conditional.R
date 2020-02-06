# Gompertz -------------------------------------------
sample_gompertz <- function(
  n_sim, # number of random variables to draw
  t0, # conditioning time
  a,
  b,
  gamma = 1) {

  # Inverse Cumulative Hazard sampling
  u <- runif(n_sim, 0, 1)

  1 / b * log( exp(b * t0) - b / (gamma * a) * log(1 - u)) - t0

}
