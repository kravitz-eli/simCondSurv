#' Simulate right administratively censored survival times with Weibull baseline hazard.
#' Assignment to treatment happens with equal probability
#'
#' @param N number of survival times to simulate
#' @param lambda rate paramter
#' @param rho scale parameter
#' @param beta log hazard ratio
#' @param n_events number of events, after this happens survival times are censored
#'
#' @return a dataframe with ID, time, censoring status, and arm
#' @export
#'
#' @examples
simulate_weibull <- function(N, lambda, rho, beta, n_events) {

  # covariate --> N Bernoulli trials
  x <- sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5))

  # Weibull latent event times
  v <- runif(n=N)
  Tlat <- (- log(v) / (lambda * exp(x * beta)))^(1 / rho)

  # censoring time
  C <- dplyr::nth(Tlat, n_events, order_by = Tlat)

  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)

  # data set
  tibble::tibble(id=1:N,
             time=time,
             status=status,
             x=x)
}
