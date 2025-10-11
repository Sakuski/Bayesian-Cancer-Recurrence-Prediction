library(ggplot2)

time <- seq(0, 40, length.out = 400)

h0_loglogistic <- function(t, p, lambda) {
  num <- lambda * p * (lambda * t)^(p - 1)
  den <- 1 + (lambda * t)^p
  num / den
}

p_shape <- 3.0
lambda_scale <- 0.1
baseline_hazard <- h0_loglogistic(time, p = p_shape, lambda = lambda_scale)

af_decelerated <- 1.8
af_baseline <- 1.0
af_accelerated <- 0.5

hazard_decelerated <- h0_loglogistic(time / af_decelerated, p = p_shape, lambda = lambda_scale) * (1 / af_decelerated)
hazard_baseline <- baseline_hazard
hazard_accelerated <- h0_loglogistic(time / af_accelerated, p = p_shape, lambda = lambda_scale) * (1 / af_accelerated)

plot_data_aft <- data.frame(
  time = rep(time, 3),
  hazard = c(hazard_decelerated, hazard_baseline, hazard_accelerated),
  profile = factor(
    rep(
      c("Low", "Medium", "High"),
      each = length(time)
    ),
    levels = c("High", "Medium", "Low")
  )
)

ggplot(plot_data_aft, aes(x = time, y = hazard, color = profile)) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Time",
    y = "Hazard Rate",
    color = "Risk"
  ) +
  scale_color_manual(values = c(
    "High" = "red",
    "Medium" = "blue",
    "Low" = "green"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25)
  )
