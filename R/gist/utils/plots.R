library(ggplot2)
library(glue)
library(brms)

source("./utils/helpers.R")

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
