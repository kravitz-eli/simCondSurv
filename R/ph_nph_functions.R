cdf = function(x, ...) useMethod("cond_sample", x)

cond_quantile = function(x, ...) useMethod("cond_quantile", x)

cond_sample.weibull <- function(
  u,
  t0,
  trt = 0,
  log_HR = 0,
  params
){

  with(
    params,
    (-log(1 - u) / (lambda * exp(log_HR * trt)) + t0^nu) ^ (1/nu) - t0
  )
}

u <- runif(nrow(data), 0 , 1)
class(u) <- c(class(u), "weibull", "nph")

cond_quantile.ph = unction(u, t0, trt, params){

  cond_sample(u,
              t0 = time,
              trt = trt,
              params = params)
}

##!!!! Make params carry the class assignment
cond_quantile.nph = function(u, t0, trt, params){


  ctrl_patients = which(trt == 0)
  trt_patients = which(trt == 0)

  ctrl_params = subset(params, select = endsWith(colnames(params), "[1]"))
  names(ctrl_params) = remove_bracket(names(ctrl_param))

  trt_params = subset(params, select = endsWith(colnames(params), "[2]"))
  names(trt_params) = remove_bracket(names(trt_param))

  # THis doesn't work because subsetting u gets rid of the class
  ctrl_times = cond_sample(
    u[ctrl_patients],
    t0 = t0[ctrl_patients],
    trt = trt[ctrl_patients],
    params = ctrl_params
    )

  trt_times = cond_sample(
    u[trt_patients],
    t0 = t0[ctrl_patients],
    trt = trt[ctrl_patients],
    params = ctrl_params
  )

  combined = data.frame()
  combined[ctrl_patients, ] = ctrl_times
  combined[trt_patients, ] = trt_times


}


cdf.ph = function(data, params){

  ctrl_patients = subset(data, trt == 0)
  trt_patients = subset(data, trt == 1)




}


cdf.nph = function(data, params){

  ctrl_patients = subset(data, trt == 0)
  trt_patients = subset(data, trt == 1)




}
