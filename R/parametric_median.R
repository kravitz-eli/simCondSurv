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
parametric_median <- function(x, ...) {

  UseMethod("parametric_median", x)

}


#' Method for when the user supplied distribution doesn't match one of the
#' options in the code.
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
parametric_median.default <- function(x, ...){
  warning(paste("No method for median of distribution:",
                class(x)))

}




