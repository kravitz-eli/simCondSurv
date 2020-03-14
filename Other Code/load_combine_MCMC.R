library(here)
library(magrittr)
library(purrr)
library(dplyr)
library(stringr)

# Get the names of the .Rds file that store the hzard ratios
mcmc_files = Sys.glob(here("MCMC", "*.rds")) %>%
  discard(str_detect, "combined_posterior") # Don't load the previously saved file with all HR's in a data.frame

# Get the names of the models from the .Rds file.
# Ex: home/Documents/exponential_NPH.rds ---> exponential_NPH
model_names = mcmc_files %>%
  basename() %>%
  gsub(pattern=".rds$", replacement="")

# Load the posterior hazard ratios into memory
posterior_HR = map(mcmc_files, readRDS) %>%
  setNames(model_names) %>%
  do.call(cbind, .) %>%
  as.data.frame()

saveRDS(posterior_HR, here("MCMC", "combined_posterior_HR.rds"))
write.csv(posterior_HR, file = here("MCMC", "posterior_HR.csv"), row.names = FALSE)
