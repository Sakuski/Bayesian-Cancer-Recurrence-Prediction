library(dplyr)

modify_patient <- function(df, patient_id, obs_time, adj_time = NULL) {
  mod_df <- df

  patient_rows <- mod_df[mod_df$PatientID == patient_id, ]
  n_rows <- nrow(patient_rows)

  # Add more rows if wanted observation time requires that
  if (n_rows < obs_time) {
    n_add <- obs_time - n_rows
    row_template <- patient_rows[1, ]
    max_time <- max(patient_rows$Time, na.rm = TRUE)
    new_rows <- do.call("rbind", replicate(n_add, row_template, simplify = FALSE))
    new_rows$Time <- seq(max_time + 1, length.out = n_add)
    mod_df <- rbind(mod_df, new_rows)
  }

  if (!is.null(adj_time)) {
    mod_df$AdjOn[mod_df$PatientID == patient_id] <- ifelse(
      mod_df$Time[mod_df$PatientID == patient_id] <= adj_time, 1, 0
    )
  }

  mod_df <- mod_df |>
    group_by(PatientID) |>
    arrange(Time, .by_group = TRUE) |>
    mutate(
      AdjStopFlag = if_else(AdjOn == 0 & cumsum(AdjOn) == max(cumsum(AdjOn)), 1, 0),
      StopTime = if_else(AdjStopFlag == 1, Time, NA),
      StopTime = if (all(is.na(StopTime))) NA else min(StopTime, na.rm = TRUE),
      TimeSinceAdjStopped = if_else(Time >= StopTime, Time - StopTime + 1, 0),
      TimeBeforeAdjStopped = if_else(Time <= StopTime, StopTime - Time, 0)
    ) |>
    ungroup()

  return(mod_df[mod_df$PatientID == patient_id, ])
}

average_cumulative_survival_probability <- function(fit, df, time, adj_time_mod = NULL) {
  patient_ids <- unique(df$PatientID)
  n_patients <- length(patient_ids)
  n_draws <- 10

  survival_list <- vector("list", n_patients)

  for (i in seq_along(patient_ids)) {
    pid <- patient_ids[i]
    patient_df <- modify_patient(df, patient_id = pid, obs_time = time, adj_time = adj_time_mod)
    y_prob <- posterior_epred(fit, newdata = patient_df, ndraws = n_draws)
    cum_surv <- t(apply(1 - y_prob, 1, cumprod))
    survival_list[[i]] <- cum_surv
  }

  survival_array <- array(unlist(survival_list), dim = c(n_draws, time, n_patients))
  avg_surv <- apply(survival_array, c(1, 2), mean)

  median_surv <- apply(avg_surv, 2, median)
  lower_surv  <- apply(avg_surv, 2, quantile, probs = 0.05)
  upper_surv  <- apply(avg_surv, 2, quantile, probs = 0.95)

  result_df <- data.frame(
    Time = 1:time,
    Median = round(median_surv, 4),
    Lower  = round(lower_surv, 4),
    Upper  = round(upper_surv, 4)
  )

  result_df
}

average_treatment_effect <- function(fit, df, time) {
  patient_ids <- unique(df$PatientID)
  n_patients <- length(patient_ids)
  n_draws <- 10

  effect_list <- vector("list", n_patients)

  for (i in seq_along(patient_ids)) {
    pid <- patient_ids[i]

    # Scenario 1: Adjuvant treatment ON for 3 years
    treated_df <- modify_patient(df, patient_id = pid, obs_time = time, adj_time = 3)
    treated_probs <- posterior_epred(fit, newdata = treated_df, ndraws = n_draws)
    treated_surv <- t(apply(1 - treated_probs, 1, cumprod))

    # Scenario 2: Adjuvant treatment OFF
    untreated_df <- modify_patient(df, patient_id = pid, obs_time = time, adj_time = 0)
    untreated_probs <- posterior_epred(fit, newdata = untreated_df, ndraws = n_draws)
    untreated_surv <- t(apply(1 - untreated_probs, 1, cumprod))

    # Difference: treated - untreated
    effect_surv <- treated_surv - untreated_surv

    effect_list[[i]] <- effect_surv
  }

  effect_array <- array(unlist(effect_list), dim = c(n_draws, time, n_patients))
  avg_effect <- apply(effect_array, c(1, 2), mean)

  median_effect <- apply(avg_effect, 2, median)
  lower_effect  <- apply(avg_effect, 2, quantile, probs = 0.05)
  upper_effect  <- apply(avg_effect, 2, quantile, probs = 0.95)

  result_df <- data.frame(
    Time   = 1:time,
    Median = round(median_effect, 4),
    Lower  = round(lower_effect, 4),
    Upper  = round(upper_effect, 4)
  )

  result_df
}
