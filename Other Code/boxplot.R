### Make a single bxoplot plot that shows the posterior distribution of all the models
##

# Load libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(here)

# Load matrix of posterior samples from each model
df <- as_tibble(readRDS(here("MCMC", "combined_posterior_HR.rds")))

# Put data is a format ggplot can understand, e.g. "tidy" data
df_tidy = df %>%
  # Change the 20,000 x 10 matrix to 200,000 x 2 matrix with indicators for model type
  # pivot_longer(everything()) %>%
  gather(key = "Model", value = "Density") %>%
  # Remove the _ from the model name. Ex: exponential_NPH ---> exponential NPH
  mutate(Model = str_replace(Model, "_", " "))


boxplot = ggplot(df_tidy, aes(x = Model, y = Density)) +
  geom_boxplot(outlier.shape=NA) +
  geom_hline(
    aes(yintercept = 1.3),
    linetype = "dashed",
    size = 1.2
  ) +
  ylim(c(0, 2)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(face = "bold")
  )

ggsave(plot = boxplot, filename = here("Plots", "boxplot.pdf"), dpi = 300)
