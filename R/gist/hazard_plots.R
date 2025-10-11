library(brms)
library(cmdstanr)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

source("./utils/plots.R")

options(brms.backend = "cmdstanr")
options(mc.cores = 8)

gist_sim_long_df <- read_delim("./../../data/gist_sim_long.csv", delim = ",")
gist_sim_short_df <- read_delim("./../../data/gist_sim_short.csv", delim = ",")

gist_sim_short_df <- gist_sim_short_df |>
  mutate(
    SizeScaled = (Size - mean(Size)) / (2 * sd(Size)),
    AgeAtSurgScaled = (AgeAtSurg - mean(AgeAtSurg)) / (2 * sd(AgeAtSurg)),
    MitHPFScaled = (MitHPF - mean(MitHPF)) / (2 * sd(MitHPF))
  )

gist_sim_long_df <- gist_sim_long_df |>
  mutate(
    SizeScaled = (Size - mean(Size)) / (2 * sd(Size)),
    AgeAtSurgScaled = (AgeAtSurg - mean(AgeAtSurg)) / (2 * sd(AgeAtSurg)),
    MitHPFScaled = (MitHPF - mean(MitHPF)) / (2 * sd(MitHPF))
  )

time_grid <- seq(0.001, 11, length.out = 100)

patient_60 <- gist_sim_short_df[gist_sim_short_df$PatientID == 60, ]

patient_no_treat <- tibble(
  AdjTreatm = 0,
  SizeScaled = patient_60$SizeScaled,
  AgeAtSurgScaled = patient_60$AgeAtSurgScaled,
  MitHPFScaled = patient_60$MitHPFScaled,
  GenderMale = patient_60$GenderMale,
  Rupture = patient_60$Rupture,
  Gastric = patient_60$Gastric,
  Censored = 0
)

patient_treat <- tibble(
  AdjTreatm = 1,
  SizeScaled = patient_60$SizeScaled,
  AgeAtSurgScaled = patient_60$AgeAtSurgScaled,
  MitHPFScaled = patient_60$MitHPFScaled,
  GenderMale = patient_60$GenderMale,
  Rupture = patient_60$Rupture,
  Gastric = patient_60$Gastric,
  Censored = 0
)

# Weibull AFT

fit_weibull <- brm(formula = EventTime | cens(Censored) ~ AdjTreatm + s(SizeScaled) + s(AgeAtSurgScaled) + s(MitHPFScaled) + GenderMale + Rupture + Gastric,
                   family = weibull(),
                   prior = c(prior(normal(0, 2), class = b)),
                   chains = 4,
                   cores = 4,
                   threads = 2,
                   control = list(adapt_delta = 0.95, max_treedepth = 15),
                   seed = 6,
                   data = gist_sim_short_df)

plot_weibull_hazard <- function(fit_weibull, new_data) {
  mean_draws_weib <- posterior_epred(fit_weibull, newdata = new_data)
  shape_draws_weib <- posterior_epred(fit_weibull, dpar = "shape", newdata = new_data)
  rate_draws_weib <- gamma(1 + 1 / shape_draws_weib) / mean_draws_weib

  weibull_hazard <- function(t, rate, shape) {
    rate * shape * ((rate * t)^(shape - 1))
  }

  hazard_matrix <- sapply(time_grid, function(t) {
    weibull_hazard(t, rate = rate_draws_weib, shape = shape_draws_weib)
  })

  hazard_long <- hazard_matrix |>
    as_tibble() |>
    setNames(time_grid) |>
    mutate(draw = seq_len(n())) |>
    pivot_longer(
      cols = -draw,
      names_to = "time_point",
      values_to = "hazard",
      names_transform = as.numeric
    )

  hazard_summary <- hazard_long |>
    group_by(time_point) |>
    summarise(
      median_hazard = median(hazard),
      lower_ci = quantile(hazard, 0.025),
      upper_ci = quantile(hazard, 0.975),
      .groups = "drop"
    )

  ggplot(hazard_summary, aes(x = time_point, y = median_hazard)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "skyblue", alpha = 0.5) +
    geom_line(color = "darkblue", linewidth = 1) +
    labs(
      x = "Time (years)",
      y = "Hazard Rate"
    ) +
    coord_cartesian(ylim = c(0, 0.08), xlim = c(0, 10)) +
    theme_minimal()
}

plot_weibull_hazard(fit_weibull, patient_treat) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25)
  )

plot_weibull_hazard(fit_weibull, patient_no_treat) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25)
  )

# Exponential

fit_exponential <- brm(formula = EventTime | cens(Censored) ~ AdjTreatm + s(SizeScaled) + s(AgeAtSurgScaled) + s(MitHPFScaled) + GenderMale + Rupture + Gastric,
                       family = exponential(),
                       prior = c(prior(normal(0, 2), class = b)),
                       chains = 4,
                       cores = 4,
                       threads = 2,
                       control = list(adapt_delta = 0.95, max_treedepth = 15),
                       seed = 5,
                       data = gist_sim_short_df)

plot_exponential_hazard <- function(fit_exponential, new_data) {
  mean_draws_exp <- posterior_epred(fit_exponential, newdata = new_data)
  rate_draws_exp <- 1 / mean_draws_exp

  exponential_hazard <- function(t, rate) {
    rate
  }

  hazard_matrix <- sapply(time_grid, function(t) {
    exponential_hazard(t, rate = rate_draws_exp)
  })

  hazard_long <- hazard_matrix |>
    as_tibble() |>
    setNames(time_grid) |>
    mutate(draw = seq_len(n())) |>
    pivot_longer(
      cols = -draw,
      names_to = "time_point",
      values_to = "hazard",
      names_transform = as.numeric
    )

  hazard_summary <- hazard_long |>
    group_by(time_point) |>
    summarise(
      median_hazard = median(hazard),
      lower_ci = quantile(hazard, 0.025),
      upper_ci = quantile(hazard, 0.975),
      .groups = "drop"
    )

  ggplot(hazard_summary, aes(x = time_point, y = median_hazard)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "skyblue", alpha = 0.5) +
    geom_line(color = "darkblue", linewidth = 1) +
    labs(
      x = "Time (years)",
      y = "Hazard Rate"
    ) +
    coord_cartesian(ylim = c(0, 0.08), xlim = c(0, 10)) +
    theme_minimal()
}

plot_exponential_hazard(fit_exponential, patient_treat) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25)
  )

plot_exponential_hazard(fit_exponential, patient_no_treat) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25)
  )

# Bernoulli model

fit_bernoulli <- brm(formula = Event ~ AdjOn + s(TimeSinceAdjStopped) + s(Time) + s(SizeScaled) + s(AgeAtSurgScaled) + s(MitHPFScaled) + GenderMale + Rupture + Gastric,
                     family = bernoulli(),
                     prior = c(prior(normal(0, 2), class = b)),
                     chains = 4,
                     cores = 4,
                     threads = 2,
                     control = list(adapt_delta = 0.95, max_treedepth = 15),
                     seed = 5,
                     data = gist_sim_long_df)

predictive_plot(fit_bernoulli, gist_sim_long_df, 60, 11, 3, ylim = c(0, 0.08)) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25)
  )

predictive_plot(fit_bernoulli, gist_sim_long_df, 60, 11, 0, ylim = c(0, 0.08)) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25)
  )
