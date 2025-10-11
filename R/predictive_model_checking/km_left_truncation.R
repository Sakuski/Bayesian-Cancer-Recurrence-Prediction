library(brms)
library(cmdstanr)
library(survival)
library(ggplot2)
library(bayesplot)

data <- local({
  set.seed(1901)
  n <- 1000
  y <- rweibull(n, shape = 1.2, scale = 1) 
  v <- rexp(n)
  data.frame(time = y, v = v)
})

data_trunc <- data[data$time >= data$v, ]

fit_brms <- brm(formula = time | trunc(lb = v, ub = Inf) ~ 1,
                family = weibull(),
                data = data_trunc,
                chains = 1,
                iter = 2000,
                seed = 5,
                backend = "cmdstanr",
                refresh = 0)

set.seed(5)
pp <- posterior_predict(fit_brms, draw_ids = 1:5)

# Without left-truncation
ppc_km_overlay(data_trunc$time, pp, status_y = rep(1, nrow(data_trunc))) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 25)
  )

# With left-truncation
ppc_km_overlay(data_trunc$time, pp, status_y = rep(1, nrow(data_trunc)), left_truncation_y = data_trunc$v) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 25)
  )
