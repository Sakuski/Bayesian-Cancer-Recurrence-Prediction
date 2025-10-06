library(brms)
library(ggplot2)
library(bayesplot)

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

# Balanced binary data

set.seed(5)
n_obs <- 250
predictor <- rnorm(n_obs, mean = 0, sd = 1.5)
true_prob <- plogis(-0.4 + 0.8 * predictor)
observed_y <- rbinom(n_obs, size = 1, prob = true_prob)
sim_data <- data.frame(observed_y, predictor)

fit_brms <- brm(
  formula = observed_y ~ predictor,
  data = sim_data,
  family = bernoulli(),
  seed = 5,
  iter = 2000,
  chains = 4
)

yrep_brms <- posterior_predict(fit_brms)

bars_plot <- ppc_bars(y = sim_data$observed_y, yrep = yrep_brms)
bars_plot <- bars_plot +
  scale_x_continuous(breaks = c(0, 1), labels = c("0.0", "1.0"))
bars_plot

ppc_calibration_pava(y = sim_data$observed_y, yrep = yrep_brms)

# Binary data where most values are near 0

set.seed(5)
n_obs <- 500
predictor <- rnorm(n_obs, mean = 0, sd = 1)
true_intercept <- -2.5
true_slope <- 0.7
true_prob <- plogis(true_intercept + true_slope * predictor)
observed_y <- rbinom(n_obs, size = 1, prob = true_prob)
low_prob_data <- data.frame(observed_y, predictor)

fit_low_prob <- brm(
  formula = observed_y ~ predictor,
  data = low_prob_data,
  family = bernoulli(),
  seed = 5,
  iter = 2000,
  chains = 4,
  silent = 2
)

yrep_low_prob <- posterior_predict(fit_low_prob)

ppc_calibration_pava(y = low_prob_data$observed_y, yrep = yrep_low_prob)
ppc_calibration_pava(y = low_prob_data$observed_y, yrep = yrep_low_prob, xlim = c(0, 0.31), ylim = c(0, 0.31))
