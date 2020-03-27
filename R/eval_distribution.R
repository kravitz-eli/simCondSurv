# Survival functions have different parameterizations than R's built in pdistribution, cdistribution, qdistribution
eval_pdf = function(x, ...) UseMethod("eval_pdf")

eval_cdf = function(x, ...) UseMethod("eval_cdf")
