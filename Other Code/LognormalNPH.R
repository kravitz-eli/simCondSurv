#########################################################################################################
#########################################################################################################
##Title: Lognormal Projection of Overall Survival in Support of RELAY Submission
##Author: Zachary Thomas
##Last Modified: 6 Feb 2020
#########################################################################################################

##Install packages
library(survival)
library(rjags)
library(here)
library(tictoc) ; tic()
##Read in dataset
dat = read.csv(here("Data", "modified_data.csv"), header = T)

# Name for the output and txt file
outname = "lognormal_NPH"

# Add paths for different file types
pdf_path = here("Plots", paste0(outname, ".pdf"))
text_path = here("Text", paste0(outname, ".txt"))
mcmc_path = here("MCMC", paste0(outname, ".rds"))

##Take appropriate rows
sub = which(dat$PARAMCD=='OS' & dat$ITTFL=='Y')

##Current OS data
t = dat$AVAL[sub]
event = 1-dat$CNSR[sub]
TRT = ifelse(dat$TRT01P[sub]=='Ramucirumab+Erlotinib',2,1)

##Median OS follow-up
median(t)

##Number enrolled patients
N = length(t)

##JAGS Model
model = '
model{

  for(i in 1:N){
    ones[i] ~ dbern(p[i]/10000000000)
    p[i] <- exp(-H[i])*pow(h[i],event[i])

    h[i] <- (1/(sigma[TRT[i]]*t[i]*sqrt(2*pi)))*exp(-(1/2)*pow((log(t[i])-mu[TRT[i]])/sigma[TRT[i]],2))/(1-phi((log(t[i])-mu[TRT[i]])/sigma[TRT[i]]))
    H[i] <- -log(1-phi((log(t[i])-mu[TRT[i]])/sigma[TRT[i]]))
  }

  mu[1]~dnorm(0,.001)
  sigma[1]~dunif(0,10)
  mu[2]~dnorm(0,.001)
  sigma[2]~dunif(0,10)

  P = prod(p)


}'

##Run JAGS to fit Weibull
input = list(t = t,event=event,N=N ,ones=rep(1,length(t)),pi = pi,TRT=TRT)
mod<- jags.model(textConnection(model),data=input,n.chains = 2)
update(mod,10000)
out.mod<- coda.samples(mod,variable.names=c('mu','sigma','P'),n.iter=100000)

##Check convergence diagnostics
# print(gelman.diag(out.mod,multivariate=T))

out = as.data.frame(as.matrix(out.mod))
names(out)

##Plot Parametric Fit vs. KM
s = seq(0,100,.1)
mE = mR = rep(NA,length(s))
for(i in 1:length(s)){
  mE[i] = median(exp(log(1-pnorm((log(s[i])-out[,2])/out[,4]))))
  mR[i] = median(exp(log(1-pnorm((log(s[i])-out[,3])/out[,5]))))
}

##Please rename and save this plot for each script
pdf(pdf_path)
KM = survfit(Surv(t,event)~TRT)
plot(KM,conf.int=F,col=c('blue','red'),xlim=c(0,100),lwd=2,xlab='Time (months)',ylab='Overall Survival Probability')
lines(s,mE,col='blue',lwd=2)
lines(s,mR,col='red',lwd=2)
abline(h=.5,lty=2)
legend("bottomleft", legend=c("Placebo","LY"),col=c("blue","red"),horiz=FALSE, lty=1)
dev.off()

##Posterior summaries for true median PFS
median(exp(out[,2]))
median(exp(out[,3]))


############################################################################################

##Inverse CDF for conditional Weibull failure times (given survival to time C)
InvCondCDF = function(u, t0, mu, sigma, beta){
  exp(
    sigma * qnorm(1 - (1 - u)^(1 / beta) * pnorm(-1 * ( log(t0) - mu ) / sigma )) + mu) - t0
}

##Who is censored
R = which(event==0)

##Get the posterior predictive conditional failure times
K = nrow(out)
times = matrix(rep(NA,K*length(R)),length(R),K)
for(i in 1:nrow(out)){
  ##Probability integral transform (PIT) trick
  a = runif(length(R),0,1)
  times[TRT[R]==1,i] = InvCondCDF(a[TRT[R]==1],t[R][TRT[R]==1],out[i,2],out[i,4],1)
  times[TRT[R]==2,i] = InvCondCDF(a[TRT[R]==2],t[R][TRT[R]==2],out[i,3],out[i,5],1)
}

##Take a look at expected additional surival vs. current follow-up
# plot(t[R],apply(times,1,median))

##Add the simulated conditional failure times to the patient's current censoring times
t.append = event.append = matrix(rep(NA,length(t)*nrow(out)),length(t),nrow(out))
AT = rep(NA,nrow(out))
for(i in 1:nrow(out)){
  ##These patients already had events so that part of dataset stays fixed
  t.append[-R,i] = t[-R]
  ##Add the simulated additional PFS to these patients' current PFS censoring time
  t.append[R,i] = t[R]+times[,i]
  E = which(rank(times[,i])<=175)
  Max = max(times[E,i])
  t.append[R,i] =  ifelse(times[,i]>Max,t[R]+Max,t.append[R,i])
  event.append[,i] = event
  event.append[R[E],i] = 1
  AT[i] = Max
}

hist(AT)
median(AT)/12
mean(AT)/12

quantile(AT,c(.025,.975))/12

##Get the KM estimates (this is the actual posterior predictive distribution for that estimator)
HR = mR = mP = Upper = rep(NA,nrow(out))
for(i in 1:(nrow(out))){
  HR[i] = exp(coxph(Surv(t.append[,i],event.append[,i])~TRT)$coef)
  KM = survfit(Surv(t.append[,i],event.append[,i])~TRT)
  mP[i] = quantile(KM,.5)$quantile[1]
  mR[i] = quantile(KM,.5)$quantile[2]
  Upper[i] = summary(coxph(Surv(t.append[,i],event.append[,i])~TRT))$conf.int[4]
  # print(i)
}

# Save below ---------
saveRDS(HR, file = mcmc_path)
file_out = file(text_path, open = "wt")
sink(file_out, type = "output")
sink(file_out, type = "message")

mean(HR)
median(HR)
quantile(HR,c(.025,.975))

mean(HR<.8)
mean(HR<.9)
mean(HR<1)
mean(HR<1.1)
mean(HR<1.2)
mean(HR<1.3)

median(mP)
quantile(mP,c(.025,.975))

median(mR[is.na(mR)==F])
quantile(mR[is.na(mR)==F],c(.025,.975))

mean(Upper<.8)
mean(Upper<.9)
mean(Upper<1)
mean(Upper<1.1)
mean(Upper<1.2)
mean(Upper<1.3)

mean(log(out$P))

closeAllConnections()
