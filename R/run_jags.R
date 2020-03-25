#'  Runs a user-supplied JAGS model and returns the output
#'
#' @param distribution Character. name of the distribution of the parametric model
#' @param prop_haz logical. Should a proportional hazards model be fit (TRUE) or a
#' nonproportional hazards model
#' @param time Positive numeric. survival times
#' @param event 0/1 variable. Event happened (\code{event = 1}) or observation is
#' censored (\code{event = 0})
#' @param trt 0/1 variable. Indicates control arm (\code{trt = 0}) or
#' experimental arm (\code{trt = 1})
#' @param chains integer number of MCMC chains to run
#' @param n_iter positive integer number of MCMC iterations, default 100,000
#' @param n_burn positive integer length of burn in, default 1,000
#' @param n_adapt number of iterations for adaption in JAGS
#' @param track_variable_names 	a character vector giving the names of variables to be monitored
#' @param progress.bar should JAGS show a progress bar, "none" (default_), "text" , or "GUI"
#'
#' @return a mcmc.list of posterior draws of parameters in \code{track_variable_names}, contains matrix of size (n_iter x track_num_variables)
#' @export
#'
#' @examples
run_jags <- function(
  distribution,
  prop_haz = TRUE,
  time,
  event,
  trt,
  chains = 2,
  n_iter = 1e3L,
  n_burn = 1e3L,
  n_adapt = 1e3L,
  track_variable_names,
  progress.bar = "none"
) {

  model_file = paste(distribution, ifelse(prop_haz, "PH", "NPH"), sep = "_")

  inits <- replicate(chains, jags_init(), simplify = FALSE)
  jags_model <- rjags::jags.model(
    get_model_file(model_file),
    data = make_jags_data(time, event, trt, prop_haz),
    n.chains = chains,
    inits = jags_init,
    n.adapt = n_adapt,
    quiet = FALSE
  )
  stats::update(jags_model, n.iter = n_burn, progress.bar = "none")
  rjags::coda.samples(
    jags_model,
    variable.names = track_variable_names,
    n.iter = n_iter,
    progress.bar = progress.bar
  )
}

jags_init <- function() {
  list(
    .RNG.name = "base::Mersenne-Twister",
    .RNG.seed = as.integer(runif(n = 1, min = 1L, max = 1e9))
  )
}

get_model_file <- function(model_file) {
  . <- NULL
  fs::path_ext_set(model_file, "jags") %>%
    file.path("inst/jags", .)
  # %>% system.file(package = "simCondSurv", mustWork = TRUE)
}

make_jags_data = function(time, event, trt, prop_haz = TRUE) {

  jags_data = list("t" = time,
                   "event" = event,
                   "TRT" = trt,
                   "ones" = rep(1, length(time)))

  if (!prop_haz) {
    jags_data$TRT = jags_data$TRT + 1
  }

  #!! Add in my new data formatting function once I replace the JAGS files.!!

  return(jags_data)

}
