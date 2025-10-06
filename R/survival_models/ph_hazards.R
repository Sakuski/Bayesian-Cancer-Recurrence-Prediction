library(ggplot2)

time <- seq(0, 25, length.out = 300)

shape_param <- 2
scale_param <- 15
baseline_hazard <- dweibull(time, shape = shape_param, scale = scale_param)

hr_profile_1 <- 0.5
hr_profile_2 <- 1.0
hr_profile_3 <- 2.2

plot_data <- data.frame(
  time = rep(time, 3),
  hazard = c(
    baseline_hazard * hr_profile_1,
    baseline_hazard * hr_profile_2,
    baseline_hazard * hr_profile_3
  ),
  profile = factor(
    rep(
      c("Low", "Medium", "High"),
      each = length(time)
    ),
    levels = c("High", "Medium", "Low")
  )
)

ggplot(plot_data, aes(x = time, y = hazard, color = profile)) +
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
    plot.title = element_text(face = "bold")
  )
