fit_jags_model = function(x, ...) UseMethod("fit_jags_model", x)

fit_jags_model.default = function(x, ...){
  warning(paste("No methods for distribution:",
                class(x)))

}
