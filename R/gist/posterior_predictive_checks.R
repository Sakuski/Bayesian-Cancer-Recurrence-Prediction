library(brms)
library(cmdstanr)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(bayesplot)

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

yrep_exponential <- posterior_predict(fit_exponential)
ppc_km_overlay(y = gist_sim_short_df$EventTime, yrep = yrep_exponential, status_y = 1 - gist_sim_short_df$Censored)

# Weibull AFT

fit_weibull <- brm(formula = EventTime | cens(Censored) ~ AdjTreatm + s(SizeScaled) + s(AgeAtSurgScaled) + s(MitHPFScaled) + GenderMale + Rupture + Gastric,
                   family = weibull(),
                   prior = c(prior(normal(0, 2), class = b)),
                   chains = 4,
                   cores = 4,
                   threads = 2,
                   control = list(adapt_delta = 0.95, max_treedepth = 15),
                   seed = 5,
                   data = gist_sim_short_df)

yrep_weibull <- posterior_predict(fit_weibull)
ppc_km_overlay(y = gist_sim_short_df$EventTime, yrep = yrep_weibull, status_y = 1 - gist_sim_short_df$Censored)

# Bernoulli

ppc_calibration_pava <- function(y,
                                 p = NULL,
                                 yrep = NULL,
                                 quantiles = 100,
                                 region.method = "resampling",
                                 n.boot = 200,
                                 dot_scale = .25,
                                 fill_alpha = .8,
                                 cep_line_color = "red",
                                 xlim = c(0, 1.01),
                                 ylim = c(0, 1.01)) {
  require("reliabilitydiag")
  require("ggdist")
  require("dplyr")

  if (is.matrix(yrep) && (nrow(yrep) > 1)) {
    bayesplot:::validate_predictions(yrep, length(y))
    if (nrow(yrep) < 100) {
      cli::cli_inform(paste(
        "Computing consistency intervals using",
        nrow(yrep),
        "predictive samples. We recommend using more predictions for more accurate bounds."
      ))
    }
    ep <- colMeans(yrep)
    ep_o <- order(ep)
    ep_s <- sort(ep)
    ceps <- Iso::pava(y[ep_o])
    consistency_intervals <- seq_len(nrow(yrep)) |>
      lapply(\(k) data.frame(
        y = Iso::pava(yrep[k, ep_o]),
        x = ep_s,
        id_ = ep_o
      )) |>
      bind_rows() |>
      group_by(id_) |>
      summarise(
        upper = quantile(y, .95),
        lower = quantile(y, .05),
        x = mean(x)
      ) |>
      arrange(x)
  } else {
    ceps <- Iso::pava(y[order(p)])
    consistency_intervals <- reliabilitydiag::reliabilitydiag(
      y = y,
      x = p,
      region.method = region.method,
      region.position = "diagonal",
      n.boot = n.boot
    )$x$regions |>
      reframe(
        lower = rep(lower, n),
        upper = rep(upper, n),
        x = rep(x, n)
      ) |>
      select(lower, upper, x)
  }
  ggplot(consistency_intervals) +
    aes(
      x = x,
      y = ceps,
      ymin = lower,
      ymax = upper
    ) +
    stat_dots(
      aes(x = x),
      quantiles = 100,
      height = dot_scale,
      scale = 1,
      shape = 19,
      colour = color_scheme_get()$mid,
      inherit.aes = FALSE
    ) +
    geom_ribbon(
      alpha = fill_alpha,
      fill = color_scheme_get()$mid
    ) +
    geom_abline(slope = 1, intercept = 0, col = "black", lty = 2, alpha = .3) +
    geom_line(colour = cep_line_color, linewidth = 1) +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    xlab("Predicted probability") +
    ylab("CEP") +
    theme(
      panel.grid = element_line(colour = "gray", linewidth = .2),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
}

yrep_bernoulli <- posterior_predict(fit_bernoulli)
ppc_calibration_pava(y = gist_sim_long_df$Event, yrep = yrep_bernoulli)
ppc_calibration_pava(y = gist_sim_long_df$Event, yrep = yrep_bernoulli, xlim = c(0, 0.26), ylim = c(0, 0.26))
