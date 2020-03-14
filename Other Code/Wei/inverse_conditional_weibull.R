# Inverse CDF for weibull conditional on a timepoint 
InvCondCDF = function(
  nsim = 1, # vector of numbers between 0 and 1. A vector of uniform [0,1] random variables will give you draws from the pdf
  C, # time which the survival function is conditioned on, i.e. S(t | C). This will
  nu,  # shape parameter
  lambda, # scale parameter
  beta = 1 # hazard ratio between control and treatment. 1 by default
){
  
  # browser()
  
  u = runif(nsim, 0, 1)
  t  = (-log((1-u)^(1 /beta)) / lambda + C^nu )^(1/nu) - C
  # t = (-log(1-u)/(lambda*beta)+C^nu)^(1/eta)-nu
  
  return(t)
}
