#' Draws samples from a user specified conditional survival function
#'
#' @param n_samples Positive integer. The number of samples to draw from
#' conditional survival function
#' @param distribution Character vector. The name of the conditional distribution
#' you want to sample from (FILL in options later)
#' @param t0 Nonnegative real number. The time the survival function is
#' conditioned on. \code{t0 = 0} will return samples from the unconditional distribution
#' @param HR hazard ratio between treatment arm and control arm
#' @param params named list of parameters for distribution
#'
#' @return a vector of conditional survival times from the specified distribution.
#' @export
#'
#' @examples
draw_samples = function(
  n_samples,
  distribution,
  t0,
  HR,
  params
){

  u <- runif(n_samples, 0, 1)

  class(u) <- distribution


}
