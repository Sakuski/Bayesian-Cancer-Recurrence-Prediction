library(ggplot2)
library(glue)
library(brms)

source("./utils/helpers.R")

predictive_plot <- function(fit, df, patient_id, time, adj_time_mod = NULL, ylim = NULL) {
  mod_df <- modify_patient(df, patient_id, time, adj_time_mod)
  patient_data <- mod_df[mod_df$PatientID == patient_id, ]

  y_prob <- posterior_epred(fit, newdata = patient_data, draws = 1000)

  y_median <- apply(y_prob, 2, median)
  y_lower <- apply(y_prob, 2, quantile, probs = 0.05)
  y_upper <- apply(y_prob, 2, quantile, probs = 0.95)

  n <- length(y_median)

  plot_df <- data.frame(
    time = patient_data$Time - 1,
    y_median = y_median,
    y_lower = y_lower,
    y_upper = y_upper,
    y_obs = patient_data$Event
  )

  plot_df <- plot_df[plot_df$time + 1 <= time, ]
  plot_df <- plot_df[order(plot_df$time), ]

  plot <- ggplot(plot_df, aes(x = time)) +
    geom_step(aes(y = y_median), direction = "hv", color = "darkblue") +
    geom_step(aes(y = y_lower), direction = "hv", color = "darkblue", alpha = 0.6, linetype = "dashed") +
    geom_step(aes(y = y_upper), direction = "hv", color = "darkblue", alpha = 0.6, linetype = "dashed") +
    scale_x_continuous(breaks = 0:n) +
    labs(y = "Probability of Event", x = "Time (years)") +
    theme_minimal()

  if (!is.null(ylim)) {
    plot <- plot + coord_cartesian(ylim = ylim)
  }

  return(plot)
}

plot_survival_comparison <- function(fit, df, time) {
  patient_ids <- unique(df$PatientID)
  n_draws <- 100

  results_df <- data.frame(
    PatientID = patient_ids,
    WithTreatment = NA_real_,
    WithoutTreatment = NA_real_
  )

  for (i in seq_along(patient_ids)) {
    pid <- patient_ids[i]

    # Scenario 1: Adjuvant treatment for 3 years
    treated_df <- modify_patient(df, patient_id = pid, obs_time = time, adj_time = 3)
    treated_probs <- posterior_epred(fit, newdata = treated_df, ndraws = n_draws)
    treated_surv <- t(apply(1 - treated_probs, 1, cumprod))
    treated_median <- median(treated_surv[, time])

    # Scenario 2: No adjuvant treatment
    untreated_df <- modify_patient(df, patient_id = pid, obs_time = time, adj_time = 0)
    untreated_probs <- posterior_epred(fit, newdata = untreated_df, ndraws = n_draws)
    untreated_surv <- t(apply(1 - untreated_probs, 1, cumprod))
    untreated_median <- median(untreated_surv[, time])

    results_df$WithTreatment[i] <- round(treated_median, 4)
    results_df$WithoutTreatment[i] <- round(untreated_median, 4)
  }

  ggplot(results_df, aes(x = results_df$WithoutTreatment, y = results_df$WithTreatment)) +
    geom_point(alpha = 0.7, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_x_continuous(limits = c(0.0, 1.0)) +
    scale_y_continuous(limits = c(0.0, 1.0)) +
    labs(
      x = "Survival Probability Without Treatment",
      y = "Survival Probability With Treatment",
      title = "Heterogeneity in Treatment Effect on Survival",
      subtitle = glue("Each point is a patient {time} years after surgery")
    ) +
    coord_fixed() +
    theme_minimal()
}
