library(survival)
library(brms)
library(cmdstanr)
library(dplyr)

ignore <- with(rotterdam, recur == 0 & death == 1 & rtime < dtime)
rfs <- with(rotterdam, ifelse(recur == 1 | ignore, recur, death))
rfstime <- with(rotterdam, ifelse(recur == 1 | ignore, rtime, dtime))

rotterdam_df <- rotterdam
rotterdam_df$rfs <- rfs
rotterdam_df$rfstime <- rfstime

rotterdam_df <- rotterdam_df[rotterdam_df$rfs == 1, ]

fit_weibull <- brm(formula = rfstime ~ 1,
                   family = weibull(),
                   chains = 4,
                   cores = 4,
                   threads = 2,
                   data = rotterdam_df)

# Intervals plot
pp_check(fit_weibull, type = "intervals", ndraws = 50)

# PIT-ECDF plot
pp_check(fit_weibull, type = "pit_ecdf")

# Kaplan-Meier overlay plot
pp_check(fit_weibull, status_y = rotterdam_df$rfs, type = "km_overlay")
