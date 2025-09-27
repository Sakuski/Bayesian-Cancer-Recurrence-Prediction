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

# Exponential

fit_exponential <- brm(formula = EventTime | cens(Censored) ~ AdjTreatm + s(SizeScaled) + s(AgeAtSurgScaled) + s(MitHPFScaled) + GenderMale + Rupture + Gastric,
                       family = exponential(),
                       prior = c(prior(normal(0, 2), class = b)),
                       chains = 4,
                       cores = 4,
                       threads = 2,
                       iter = 4000,
                       control = list(adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE),
                       seed = 5,
                       data = gist_sim_short_df)

loo_fit_exponential <- loo(fit_exponential)

# Cox PH

fit_cox <- brm(formula = EventTime | cens(Censored) ~ AdjTreatm + s(SizeScaled) + s(AgeAtSurgScaled) + s(MitHPFScaled) + GenderMale + Rupture + Gastric,
               family = cox(),
               prior = c(prior(normal(0, 2), class = b)),
               chains = 4,
               cores = 4,
               threads = 2,
               iter = 4000,
               control = list(adapt_delta = 0.99),
               save_pars = save_pars(all = TRUE),
               seed = 5,
               data = gist_sim_short_df)

loo_fit_cox <- loo(fit_cox)

# Weibull AFT

fit_weibull <- brm(formula = EventTime | cens(Censored) ~ AdjTreatm + s(SizeScaled) + s(AgeAtSurgScaled) + s(MitHPFScaled) + GenderMale + Rupture + Gastric,
                   family = weibull(),
                   prior = c(prior(normal(0, 2), class = b)),
                   chains = 4,
                   cores = 4,
                   threads = 2,
                   iter = 4000,
                   control = list(adapt_delta = 0.95, max_treedepth = 15),
                   save_pars = save_pars(all = TRUE),
                   seed = 5,
                   data = gist_sim_short_df)

loo_fit_weibull <- loo(fit_weibull)

# Comparison

loo_compare(
  list(exponential = loo_fit_exponential, cox = loo_fit_cox, weibull = loo_fit_weibull)
)
