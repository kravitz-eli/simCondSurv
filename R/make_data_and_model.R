# Creates a data_and_model data structure that holds distribution, PH vs NPH,
# posterior samples, paramater names, index of NPH parameters if using
# NPH model,

make_data_and_model <- function(data,
                                model,
                                distribution,
                                prop_haz) {
  data_and_model <- list(
    "time" = data$time,
    "event" = data$event,
    "trt" = data$trt,
    "jags_model" = model,
    "post_params" = tibble::as_tibble(model[[1]])
  )


  # This object has two classes. One indicates PH vs NPH and one indicates distribution
  attr(data_and_model, "class") <- c(ifelse(prop_haz, "PH", "NPH"), distribution)

  data_and_model$param_names <- get_paramater_names(data_and_model)

  data_and_model$ctrl_param_index <- which(endsWith(names(data_and_model$post_params), "[1]"))
  data_and_model$trt_param_index <- which(endsWith(names(data_and_model$post_params), "[2]"))

  return(data_and_model)

}
