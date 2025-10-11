library(survival)
library(brms)
library(cmdstanr)
library(ggplot2)

ignore <- with(rotterdam, recur == 0 & death == 1 & rtime < dtime)
rfs <- with(rotterdam, ifelse(recur == 1 | ignore, recur, death))
rfstime <- with(rotterdam, ifelse(recur == 1 | ignore, rtime, dtime))

rotterdam_df <- rotterdam
rotterdam_df$rfs <- rfs
rotterdam_df$rfstime <- rfstime
rotterdam_df$rfsmonths <- rotterdam_df$rfstime / 30.44

# Event time in days

fit_weibull_days <- brm(formula = rfstime | cens(1 - rfs) ~ 1,
                        family = weibull(),
                        chains = 4,
                        cores = 4,
                        threads = 2,
                        seed = 5,
                        data = rotterdam_df)

loo_fit_weibull_days <- loo(fit_weibull_days)

pointwise_elpd_days <- loo_fit_weibull_days$pointwise[, "elpd_loo"]
rotterdam_df$elpd_days <- pointwise_elpd_days

df_elpd_days <- data.frame(
  Observation = seq_along(pointwise_elpd_days),
  elpd = pointwise_elpd_days
)

ggplot(df_elpd_days, aes(x = elpd)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_vline(aes(xintercept = mean(elpd)), color = "red", linetype = "dashed") +
  labs(
    x = "Pointwise ELPD",
    y = "Count"
  ) +
  coord_cartesian(xlim = c(-10, 0), ylim = c(0, 1000)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25)
  )

# Event time in months

fit_weibull_months <- brm(formula = rfsmonths | cens(1 - rfs) ~ 1,
                          family = weibull(),
                          chains = 4,
                          cores = 4,
                          threads = 2,
                          seed = 5,
                          data = rotterdam_df)

loo_fit_weibull_months <- loo(fit_weibull_months)

pointwise_elpd_months <- loo_fit_weibull_months$pointwise[, "elpd_loo"]
rotterdam_df$elpd_months <- pointwise_elpd_months

df_elpd_months <- data.frame(
  Observation = seq_along(pointwise_elpd_months),
  elpd = pointwise_elpd_months
)

ggplot(df_elpd_months, aes(x = elpd)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_vline(aes(xintercept = mean(elpd)), color = "red", linetype = "dashed") +
  labs(
    x = "Pointwise ELPD",
    y = "Count"
  ) +
  coord_cartesian(xlim = c(-10, 0), ylim = c(0, 1000)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25)
  )
