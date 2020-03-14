#'  Runs a user-supplied JAGS model and returns the output
#'
#' @param model_file string containing name of JAGS model to run. Must be stored in jags/ for now
#' @param data list of parameters and data to pass to JAGS
#' @param chains integer number of MCMC chains to run
#' @param iter integer number of MCMC iterations, default 100,000
#' @param burn integer length of burnin, default 1,000
#' @param track_variable_names 	a character vector giving the names of variables to be monitored
#' @param progress.bar should JAGS show a progress bar, "none" (default_), "text" , or "GUI"
#'
#' @return a mcmc.list of posterior draws of parameters in \code{track_variable_names}, contains matrix of size (iter x track_num_variables)
#' @export
#'
#' @examples
run_jags <- function(
  model_file,
  data,
  chains,
  iter = 1e4L,
  burn = 1e3L,
  track_variable_names,
  progress.bar = "none"
) {
  
  inits <- replicate(chains, jags_init(), simplify = FALSE)
  jags_model <- rjags::jags.model(
    get_model_file(model_file),
    data = data,
    n.chains = chains,
    inits = jags_init,
    n.adapt = burn,
    quiet = TRUE
  )
  stats::update(jags_model, n.iter = burn, progress.bar = "none")
  rjags::coda.samples(
    jags_model,
    variable.names = track_variable_names,
    n.iter = iter,
    progress.bar = progress.bar
  )
}

jags_init <- function() {
  list(
    .RNG.name = "base::Mersenne-Twister",
    .RNG.seed = as.integer(runif(n = 1, min = 1L, max = 1e9))
  )
}

get_model_file <- function(model_name, jags_path = "JAGS") {
  paste0(model_name, ".jags") %>%
  here::here(jags_path, .)
}
