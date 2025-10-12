library(ggplot2)
library(glue)
library(brms)

source("./utils/helpers.R")

predictive_plot <- function(fit, df, patient_id, time, adj_time_mod = NULL, ylim = NULL) {
  mod_df <- modify_patient(df, patient_id, time, adj_time_mod)
  patient_data <- mod_df[mod_df$PatientID == patient_id, ]

  y_prob <- posterior_epred(fit, newdata = patient_data, draws = 1000)

  y_median <- apply(y_prob, 2, median)
  y_lower <- apply(y_prob, 2, quantile, probs = 0.05)
  y_upper <- apply(y_prob, 2, quantile, probs = 0.95)

  n <- length(y_median)

  plot_df <- data.frame(
    time = patient_data$Time - 1,
    y_median = y_median,
    y_lower = y_lower,
    y_upper = y_upper,
    y_obs = patient_data$Event
  )

  plot_df <- plot_df[plot_df$time + 1 <= time, ]
  plot_df <- plot_df[order(plot_df$time), ]

  plot <- ggplot(plot_df, aes(x = time)) +
    geom_step(aes(y = y_median), direction = "hv", color = "darkblue") +
    geom_step(aes(y = y_lower), direction = "hv", color = "darkblue", alpha = 0.6, linetype = "dashed") +
    geom_step(aes(y = y_upper), direction = "hv", color = "darkblue", alpha = 0.6, linetype = "dashed") +
    scale_x_continuous(breaks = 0:n) +
    labs(y = "Probability of Event", x = "Time (years)") +
    theme_minimal()

  if (!is.null(ylim)) {
    plot <- plot + coord_cartesian(ylim = ylim)
  }

  return(plot)
}
