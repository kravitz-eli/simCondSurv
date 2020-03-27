make_plots = function(
  params,
  data,
  post_pred_HR,
  criticial_value,
  axis_min_time = 0,
  axis_max_time = 100
){

  survival_plot = plot_over_KM(
    data = data,
    params = params,
    distribution = distribution,
    axis_min_time = axis_min_time,
    axis_max_time = axis_max_time
  )

  HR_plot = plot_hazard_ratio(HR = post_pred_HR)


  return(list("survival_plot" = survival_plot,
              "HR_plot" = HR_plot))

}

plot_over_KM = function(
  data,
  params,
  distribution,
  axis_min_time,
  axis_max_time
){

  #!!! Add check to see if ggplot and survminer are installed!!
  # requireNamespace("survminer", quietly = TRUE) & requireNamespace("ggplot2", quietly = TRUE)`

  # Get a grid of timepoints and calculate the survival function at every point
  time_seq = seq(axis_min_time, axis_max_time, by = 0.5)
  class(params) <- c(distribution, class(params)) # Have to manually set the class for now

  post_median_ctrl = vector('numeric', length(time_seq))
  post_median_trt = vector('numeric', length(time_seq))
  for(i in seq_along(time_seq)) {
    post_median_ctrl[[i]] = median(1 - eval_cdf(params = params,
                                                t = time_seq[[i]],
                                                trt = 0))

    post_median_trt[[i]] = median(1 - eval_cdf(params = params,
                                               t = time_seq[[i]],
                                               trt = 1))
  }


  plotting_grid = data.frame(time_seq, post_median_ctrl, post_median_trt)


  # Plot settings -----------------------------
  color_ctrl = "blue"
  color_trt = "red"

  # Make plot with survminer -----------------------
  KM_fit = survfit(Surv(time, event) ~ trt, data = data)
  KM_plot <- survminer::ggsurvplot(
    KM_fit,
    data = data,
    palette = c(color_ctrl, color_trt),
    ggtheme = theme_bw(),
    xlim = c(axis_min_time, axis_max_time)
  )

  # Edit the survminer plot object

  KM_plot$plot <- KM_plot$plot +
    geom_line(
      data = plotting_grid,
      mapping = aes(x = time_seq, y = post_median_ctrl),
      size = 1.2,
      linetype = "dashed",
      color = color_ctrl
    ) +
    geom_line(
      data = plotting_grid,
      mapping = aes(x = time_seq, y = post_median_trt),
      size = 1.2,
      linetype = "dashed",
      color = color_trt
    ) +
    scale_x_continuous()


}

plot_hazard_ratio = function(HR, critical_value = NULL){

  ggplot(
    data = data.frame("Hazard Ratio" = HR, check.names = FALSE),
    mapping = aes(x = `Hazard Ratio`)
  ) +
    geom_histogram(
      mapping = aes(y=..density..),
      colour = "white",
      alpha = 0.55
      ) +
    stat_density(
      geom = "line",
      size = 1.5,
      adjust = 1.5,
    ) +
    theme_bw()

  # # Shade the area between 0 and the critical value and add a dashed line
  # # at the critical value
  # if (is.null(!critical_value)) {
  #
  #   plot_settings <- ggplot_build(HR_plot)
  #   # Find the index of data frame where the HR is 0
  #   x1 <- min(which(plot_settings$data[[1]]$x >= 0))
  #   # Find the index where the data frame where HR < critical value
  #   x2 <- max(which(plot_settings$data[[1]]$x <= critical_value))
  #
  #   HR_plot +
  #     geom_area(
  #       data = data.frame(
  #         x = plot_settings$data[[1]]$x[x1:x2],
  #         y = plot_settings$data[[1]]$y[x1:x2]
  #       ),
  #       aes(x = x, y = y),
  #       fill = "#8EE5EE",
  #       alpha = 0.8
  #     ) +
  #     geom_linerange(
  #       x = critical_value,
  #       ymin = 0,
  #       ymax = plot_settings$data[[1]]$y[x2],
  #       linetype = "dashed",
  #       size = 1
  #     )
  # }

}
