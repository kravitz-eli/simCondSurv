plot_surv = function(
  beta, 
  lambda, 
  nu, 
  data,
  upper_bound
){
  
  ##Plot Parametric Fit vs. KM
  s = seq(0, upper_bound, 0.01)
  
  # control group posterior preditive surivival time
  m_ctrl = rep(NA,length(s))
  # treatment group posterior preditive surivival time
  m_trt = rep(NA, length(s))
  
  for(i in seq_along(s)) {
    ##Curve for control arm
    m_ctrl[i] = median(exp(-lambda*s[i]^nu))
    
    ##Curve for ram arm
    m_trt[i] = median(exp(-lambda*s[i]^nu * exp(beta)))
  }
  
  KM = survfit(Surv(time,event) ~ trt, data = data)
  plot(
    KM,
    conf.int=F,
    col=c('blue','red'),
    xlim= c(0, upper_bound),
    lwd=2,
    xlab='Time (months)',
    ylab='Overall Survival Probability'
  )
  lines(s,m_ctrl,col='blue',lwd=2, lty = 3)
  lines(s,m_trt,col='red',lwd=2, lty = 3)
  abline(h=.5,lty=2)
  legend(
    'bottomleft',
    c('Trt Observed', 'Trt Estimated', 'Ctrl Observed', 'Ctrl Estimated'),
    lwd = 2, 
    col = c('red','red', 'blue', 'blue'),
    lty = c(1, 2, 1, 2))
  
  
  
}