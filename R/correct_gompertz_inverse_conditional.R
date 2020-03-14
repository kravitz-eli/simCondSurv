InvCondCDF = function(u, t0, a, b, beta){

  1/b * (log((b * ((a * exp(b * t0))/b - log((1 - u)^ (1 / beta))))/ a) - b * t0)
}


test1 = function(u, t0, a, b, beta){
1 / b * log( exp(b * t0) - b / (beta * a) * log((1 - u)^(1/beta))) - t0
}



a = 0.2
b = 0.45
beta = 1.5

samp = InvCondCDF(
  u = seq(0.01, 0.99, by = 0.01),
  t0 = 0,
  a = a,
  b = b,
  beta = beta
)

samp1 = test1(
  u = seq(0.01, 0.99, by = 0.01),
  t0 = 0,
  a = a,
  b = b,
  beta = beta
)

samp2 = qgompertz(
  p = seq(0.01, 0.99, by = 0.01),
  alpha = b,
  theta = a
)

par(mfrow = c(1,2))
plot(density(samp))
plot(density(samp2))

library(reliaR)
curve(qgompertz(x, alpha = b, theta = a), from = 0, to = 1)

dgompertz(, alpha = b, theta = a)



samp = InvCondCDF(
  u = seq(0.01, 0.99, by = 0.01),
  t0 = 0,
  a = a,
  b = b,
  beta = 1
)
