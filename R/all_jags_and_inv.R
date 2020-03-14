# Lognormal --------------

##JAGS Model
model = '
model{

  for(i in 1:N){
    ones[i] ~ dbern(p[i]/10000000000)
    p[i] <- exp(-H[i])*pow(h[i],event[i])

    h[i] <- exp((TRT[i]-1)*beta)*(1/(sigma*t[i]*sqrt(2*pi)))*exp(-(1/2)*pow((log(t[i])-mu)/sigma,2))/(1-phi((log(t[i])-mu)/sigma))
    H[i] <- exp((TRT[i]-1)*beta)*(-log(1-phi((log(t[i])-mu)/sigma)))
  }

  mu~dnorm(0,.001)
  sigma~dunif(0,10)
  beta~dnorm(0,.001)
  P = prod(p)


}'

InvCondCDF = function(u, t0, mu, sigma, beta){
  exp(
    sigma * qnorm(1 - (1 - u)^(1 / beta) * pnorm(-1 * ( log(t0) - mu ) / sigma )) + mu) - t0
}

# Weibull -----------
##JAGS Model
model = '
model{

  for(i in 1:N){
    ##Bernoulli is used to deal with custom joint likelihood with right censoring
    ones[i]~dbern(p[i]/10000000)

    ## ith contribution to the joint likelihood
    p[i] <- exp(-H[i])* pow(h[i],event[i])

    ## hazard and cumulative hazard functions
    h[i] <- nu[TRT[i]]*lambda[TRT[i]]*pow(t[i],nu[TRT[i]]-1)
    H[i] <- lambda[TRT[i]]*pow(t[i],nu[TRT[i]])
  }

  ##Priors
  nu[1]~dnorm(0,.1)T(0,)
  lambda[1]~dgamma(.1,.1)

  nu[2]~dnorm(0,.1)T(0,)
  lambda[2]~dgamma(.1,.1)

  P = prod(p)

}'


##Inverse CDF for conditional Weibull failure times (given survival to time C)
InvCondCDF = function(u,C,eta,lambda,beta){
  ##u uniform sample per use of PIT thm
  ##C is current censor time
  ##eta is a sample of Weibull shape (called 'nu' in JAGS output)
  ##lambda is a sample of Weibull scale
  t  = (-log(1-u)/(lambda*beta)+C^eta)^(1/eta)-C
  return(t)
}


# Log logistic ----------------------------

##JAGS Model
model = '
model{

  for(i in 1:N){
    ones[i] ~ dbern(p[i]/10000000000)
    p[i] <- S[i]*pow(h[i],event[i])

    h[i] <- (beta[TRT[i]]/alpha[TRT[i]])*pow(t[i]/alpha[TRT[i]],beta[TRT[i]]-1)/(1+pow(t[i]/alpha[TRT[i]],beta[TRT[i]]))
    S[i] <- 1/(1+pow(t[i]/alpha[TRT[i]],beta[TRT[i]]))
  }

  alpha[1] ~ dnorm(0,.001)T(0,)
  alpha[2] ~ dnorm(0,.001)T(0,)
  beta[1] ~ dnorm(0,.001)T(0,)
  beta[2] ~ dnorm(0,.001)T(0,)

  P = prod(p)


}'

##Inverse CDF for conditional Weibull failure times (given survival to time C)
InvCondCDF = function(u, t0, alpha, beta, lambda){
  out=alpha*((1+(t0/alpha)^beta)*(1-u)^(-1/lambda)-1)^(1/beta)-t0
  return(out)
}

# Gompertz -------------------
##JAGS Model
model = '
model{

  for(i in 1:N){

    ones[i] ~ dbern(p[i]/10000000000)
    p[i] <- exp(-H[i])*pow(h[i],event[i])

    h[i] <- a[TRT[i]]*exp(b[TRT[i]]*t[i])
    H[i] <- (a[TRT[i]]/b[TRT[i]])*(exp(b[TRT[i]]*t[i])-1)

  }

  a[1] ~ dnorm(0,1)T(0,)
  a[2] ~ dnorm(0,1)T(0,)
  b[1] ~ dnorm(0,1)T(0,)
  b[2] ~ dnorm(0,1)T(0,)

  P = prod(p)


}'


############################################################################################

InvCondCDF = function(u, t0, a, b, beta){
  # 1 / b * log( exp(b * t0) - b / (beta * a) * log(1 - u)) - t0

  # Wolfram alpha's crazy function:
  1/b * (log((b * ((a * exp(b * t0))/b - log((1 - u)^ (1 / beta))))/ a) - b * t0)

  # Simple function that matches Wolfram alpha
  # 1 / b * log(1 - b/(a * beta) * log(1 - u) * exp(-b * t0))


}
