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
