library(survival)
library(brms)
library(cmdstanr)
library(bayesplot)

rotterdam_df <- rotterdam
rotterdam_df$rfs <- pmax(rotterdam$recur, rotterdam$death)
rotterdam_df$rfstime <- with(rotterdam, ifelse(recur == 1, rtime, dtime))

fit_weibull <- brm(formula = rfstime | cens(1 - rfs) ~ 1,
                   family = weibull(),
                   chains = 4,
                   cores = 4,
                   threads = 2,
                   seed = 5,
                   data = rotterdam_df)

set.seed(5)
yrep <- posterior_predict(fit_weibull, ndraws = 1000)

# Unreasonable amount of extrapolation
ppc_km_overlay(y = rotterdam_df$rfstime, yrep = yrep, status_y = rotterdam_df$rfs, extrapolation_factor = Inf) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 25)
  )

# Controlled extrapolation
ppc_km_overlay(y = rotterdam_df$rfstime, yrep = yrep, status_y = rotterdam_df$rfs) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 25)
  )

yimp <- sapply(seq_len(nrow(rotterdam)), \(i) {
                                               ifelse(rotterdam_df$rfs[i] == 0, sample(yrep[yrep[, i] > rotterdam_df$rfstime[i], i], size = 1), rotterdam_df$rfstime[i])})

# Controlled extrapolation with imputation (requires modifying ppc_km_overlay to plomt survival curve for yimp)
ppc_km_overlay(y = rotterdam_df$rfstime, yrep = yrep, yimp = yimp, status_y = rotterdam_df$rfs) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 25)
  )
