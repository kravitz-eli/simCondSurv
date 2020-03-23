#'  Simulate random variables from a Gompertz distribution. The
#'  function is parameterized with hazard: \eqn{\lambda(x) =  a \exp(bx)}
#'
#' @inheritParams cond_sample
#' @param params list with elements a and b, the parameters of a gompertz distribution.
#'
#' @return
#' @export
#'
#' @examples
cond_sample.gompertz <- function(
  u,
  t0,
  params,
  HR = 1
)
{

  1 / b * log( exp(b * t0) - b / (beta * a) * log(1 - u)) - t0


}
