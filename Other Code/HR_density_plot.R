# Load libraries --------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(here)

# Load and format data --------------------------------------------------------

# Load matrix of posterior samples from each model
df <- as_tibble(readRDS(here("MCMC", "combined_posterior_HR.rds")))

# Get the Weibull NPH posterior draws, discard everything else
weibull_NPH_HR <- df %>%
  select("Hazard Ratio" = weibull_NPH)

# Create density plot ------------------------------------------------------

# Make a density plot. This won't have a shaded region from 0 to 1 yet
density_plot <- ggplot(
  data = weibull_NPH_HR,
  mapping = aes(x = `Hazard Ratio`)
) +
  # Plot the kernel density estimate
  geom_density(
    color = "blue",
    size = 1.2,
    adjust = 2 # default bandwith too small, make bandwith twice as big
  ) +
  # ggplot adds a line along x-axis to density plot, remove this by adding a
  # white line on top
  geom_hline(
    yintercept = 0,
    color = "white",
    size = 1.2
  ) +
  # Draw line at x = 1
  geom_vline(
    aes(xintercept = 1),
    linetype = "dashed"
  ) +
  # Draw line at x = 1.3
  geom_vline(
    aes(xintercept = 1.3),
    linetype = "dashed"
  ) +
  theme_classic()

# Add a shaded region from 0 to 1 --------
# Get all the ggplot settings
plot_settings <- ggplot_build(density_plot)

# Find the index of data frame where the HR is 0
x1 <- min(which(plot_settings$data[[1]]$x >= 0))
# Find the index where the data frame where HR is 1
x2 <- max(which(plot_settings$data[[1]]$x <= 1))

# Add a shaded region under the curve from 0 to 1
plot_with_shading = density_plot +
  geom_area(
    data = data.frame(
      x = plot_settings$data[[1]]$x[x1:x2],
      y = plot_settings$data[[1]]$y[x1:x2]
    ),
    aes(x = x, y = y),
    fill = "grey",
    alpha = 0.8
  )

ggsave(plot = plot_with_shading,
       file = here("Plots", "weibull_posterior_HR.pdf"))
