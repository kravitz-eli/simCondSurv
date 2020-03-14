# simulate conditional Weibull conditional on survival > T ---------------

# survival function is exp{-(T+t/b)^a} / exp{-(T/b)^a} = 1-F(t)
# n = number of points to return
# shape = shape parm of weibull
# scale = scale parm of weibull (default 1)
# t is minimum (default is 0, which makes the function act like rweibull)
my_rweibull <- function(n,shape,scale=1,t=0) {
  if (length(t)!=1 && length(t)!=n) {
    stop("length(t) is not 1 or n")
  }
  return(shape*(-log(runif(n))+(t/shape)^scale)^(1/scale))
}
