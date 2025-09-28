library(brms)
library(cmdstanr)
library(readr)
library(dplyr)
library(loo)

options(brms.backend = "cmdstanr")
options(mc.cores = 8)

gist_sim_short_df <- read_delim("./../../data/gist_sim_short.csv", delim = ",")

gist_sim_short_df <- gist_sim_short_df |>
  mutate(
    SizeScaled = (Size - mean(Size)) / (2 * sd(Size)),
    AgeAtSurgScaled = (AgeAtSurg - mean(AgeAtSurg)) / (2 * sd(AgeAtSurg)),
    MitHPFScaled = (MitHPF - mean(MitHPF)) / (2 * sd(MitHPF))
  )

gist_sim_short_df <- gist_sim_short_df |>
  mutate(EventBefore5 = if_else(EventTime <= 5, 1, 0))

true_outcomes <- gist_sim_short_df$EventBefore5

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

mean_draws_weib <- posterior_epred(fit_weibull)
shape_draws_weib <- posterior_epred(fit_weibull, dpar = "shape")
rate_draws_weib <- gamma(1 + 1 / shape_draws_weib) / mean_draws_weib

prob_event_before_5_weib <- 1 - exp(-((5 * rate_draws_weib)^shape_draws_weib))

log_score_matrix_weib <- matrix(NA, nrow = nrow(prob_event_before_5_weib), ncol = ncol(prob_event_before_5_weib))

for (i in seq_len(ncol(prob_event_before_5_weib))) {
  p_i <- prob_event_before_5_weib[, i]
  if (true_outcomes[i] == 1) {
    log_score_matrix_weib[, i] <- log(p_i)
  } else {
    log_score_matrix_weib[, i] <- log(1 - p_i)
  }
}

chain_ids_weib <- rep(1:fit_weibull$fit@sim$chains, each = fit_weibull$fit@sim$iter - fit_weibull$fit@sim$warmup)
r_eff_weib <- relative_eff(exp(log_score_matrix_weib), chain_id = chain_ids_weib)

loo_weib <- loo(log_score_matrix_weib, r_eff = r_eff_weib)

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

mean_draws_exp <- posterior_epred(fit_exponential)
rate_draws_exp <- 1 / mean_draws_exp

prob_event_before_5_exp <- 1 - exp(-5 * rate_draws_exp)

log_score_matrix_exp <- matrix(NA, nrow = nrow(prob_event_before_5_exp), ncol = ncol(prob_event_before_5_exp))

for (i in seq_len(ncol(prob_event_before_5_exp))) {
  p_i <- prob_event_before_5_exp[, i]
  if (true_outcomes[i] == 1) {
    log_score_matrix_exp[, i] <- log(p_i)
  } else {
    log_score_matrix_exp[, i] <- log(1 - p_i)
  }
}

chain_ids_exp <- rep(1:fit_exponential$fit@sim$chains, each = fit_exponential$fit@sim$iter - fit_exponential$fit@sim$warmup)
r_eff_exp <- relative_eff(exp(log_score_matrix_exp), chain_id = chain_ids_exp)

loo_exp <- loo(log_score_matrix_exp, r_eff = r_eff_exp)

# Comparison

loo_compare(
  list(weibull = loo_weib, exponential = loo_exp)
)
