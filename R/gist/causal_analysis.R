library(brms)
library(cmdstanr)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)

source("./utils/helpers.R")
source("./utils/plots.R")
source("./utils/ros_helpers.R")

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

fit_propensity_score <- brm(formula = AdjTreatm ~ s(SizeScaled) + s(AgeAtSurgScaled) + s(MitHPFScaled) + GenderMale + Rupture + Gastric,
                            family = bernoulli(),
                            prior = c(prior(normal(0, 2), class = b)),
                            chains = 4,
                            cores = 4,
                            threads = 2,
                            control = list(adapt_delta = 0.99, max_treedepth = 15),
                            data = gist_sim_short_df)

summary(fit_propensity_score, robust = TRUE)
prior_summary(fit_propensity_score)

linpred <- posterior_linpred(fit_propensity_score, transform = FALSE)
pscores <- apply(plogis(linpred), 2, mean)
matches_wr <- matching(z = gist_sim_short_df$AdjTreatm, score = pscores, replace = TRUE)
wts_wr <- matches_wr$cnts

covs <- c("SizeScaled", "AgeAtSurgScaled", "MitHPFScaled", "GenderMale", "Rupture", "Gastric")
cov_names <- c("SizeScaled", "AgeAtSurgScaled", "MitHPFScaled", "GenderMale", "Rupture", "Gastric")
num_names <- c("SizeScaled", "AgeAtSurgScaled", "MitHPFScaled")
bin_names <- c("GenderMale", "Rupture", "Gastric")

bal_wr <- balance(rawdat = as.data.frame(gist_sim_short_df[, covs]), gist_sim_short_df$AdjTreatm,
                  matched <- matches_wr$cnts, estimand = "ATT")
plot.balance(bal_wr, longcovnames = num_names, which.cov = "cont")
plot.balance(bal_wr, longcovnames = bin_names, which.cov = "binary")

sum(pscores[gist_sim_short_df$AdjTreatm == 1] > max(pscores[gist_sim_short_df$AdjTreatm == 0]))
sum(pscores[gist_sim_short_df$AdjTreatm == 0] < min(pscores[gist_sim_short_df$AdjTreatm == 1]))

logit <- function(p) log(p / (1 - p))
logit_pscores <- logit(pscores)

ps_df <- data.frame(
  logit_pscore = logit_pscores,
  treat = gist_sim_short_df$AdjTreatm,
  weight = matches_wr$cnts
)

ggplot(ps_df, aes(x = logit_pscore, color = factor(treat))) +
  geom_histogram(
    aes(y = ..density..), bins = 30,
    fill = NA, position = "identity", size = 1
  ) +
  scale_color_manual(
    values = c("gray70", "black"),
    labels = c("Control", "Treated")
  ) +
  theme_minimal() +
  labs(title = "Before Matching",
       x = "logit propensity scores",
       y = "Density",
       color = "") +
  theme(legend.position = "top")

ggplot(ps_df, aes(x = logit_pscore, color = factor(treat), weight = weight)) +
  geom_histogram(
    aes(y = ..density..), bins = 30,
    fill = NA, position = "identity", size = 1
  ) +
  scale_color_manual(
    values = c("gray70", "black"),
    labels = c("Control", "Treated")
  ) +
  theme_minimal() +
  labs(title = "After Matching",
       x = "logit propensity scores",
       y = "Density",
       color = "") +
  theme(legend.position = "top")

data_expanded_df <- gist_sim_short_df[rep(seq_len(nrow(gist_sim_short_df)), matches_wr$cnts), ]
data_expanded_df$PatientID <- seq_len(nrow(data_expanded_df))

data_expanded_long_df <- data_expanded_df |>
  rowwise() |>
  mutate(data = list(tibble(Time = 1:EventTime))) |>
  unnest(data) |>
  ungroup() |>

  mutate(
    AdjOn = if_else(AdjTreatm == 1 & Time <= 3, 1, 0)
  ) |>

  group_by(PatientID) |>
  mutate(
    TimeSinceAdjStopped = {
      out <- integer(n())
      counter <- 0
      for (i in seq_len(n())) {
        if (AdjOn[i] == 1) {
          counter <- 0
        } else {
          if (i == 1 || AdjOn[i - 1] == 1) {
            counter <- 1
          } else {
            counter <- counter + 1
          }
        }
        out[i] <- counter
      }
      out
    },

    Event = if_else(Time == max(Time), 1 - Censored, 0)
  ) |>
  ungroup() |>

  select(PatientID, Size, AgeAtSurg, MitHPF, GenderMale, Rupture, Gastric, Time, AdjOn, TimeSinceAdjStopped, Event)

data_expanded_long_df <- data_expanded_long_df |>
  mutate(
    SizeScaled = (Size - mean(Size)) / (2 * sd(Size)),
    AgeAtSurgScaled = (AgeAtSurg - mean(AgeAtSurg)) / (2 * sd(AgeAtSurg)),
    MitHPFScaled = (MitHPF - mean(MitHPF)) / (2 * sd(MitHPF))
  )

fit_bernoulli_treat <- brm(formula = Event ~ AdjOn + s(TimeSinceAdjStopped) + s(Time) + s(SizeScaled) + s(AgeAtSurgScaled) + s(MitHPFScaled) + GenderMale + Rupture + Gastric,
                           family = bernoulli(),
                           prior = c(prior(normal(0, 2), class = b)),
                           chains = 4,
                           cores = 4,
                           threads = 2,
                           control = list(adapt_delta = 0.99, max_treedepth = 15),
                           data = data_expanded_long_df)

summary(fit_bernoulli_treat, robust = TRUE)
prior_summary(fit_bernoulli_treat)

treated_df <- gist_sim_long_df |>
  group_by(PatientID) |>
  filter(any(AdjOn == 1)) |>
  ungroup()

avg_surv_prob_treat <- average_cumulative_survival_probability(fit_bernoulli_treat, treated_df, 10, 3)
avg_surv_prob_treat

avg_surv_prob_no_treat <- average_cumulative_survival_probability(fit_bernoulli_treat, treated_df, 10, 0)
avg_surv_prob_no_treat

avg_treat_effect <- average_treatment_effect(fit_bernoulli_treat, treated_df, 10)
avg_treat_effect

survival_comparison_plot <- plot_survival_comparison(fit_bernoulli_treat, treated_df, 5)
survival_comparison_plot
