#' Generic function to call samplers
#'
#' @param u numeric vector, \eqn{0 \leq u  \leq 1} for all elements in u.
#' @param t0 Nonnegative real number. The time the survival function is
#' conditioned on. \code{t0 = 0} will return samples from the unconditional distribution
#' @param HR hazard ratio between treatment arm and control arm
#' @param params named list of parameters for distribution
#'
#' @return
#' @export
#'
#' @examples
cond_sample <- function(u, t0, HR, params) {
  UseMethod("cond_sample", u, t0, HR, params)
}
