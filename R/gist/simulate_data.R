library(brms)
library(readr)
library(dplyr)
library(synthpop)

set.seed(5)

# Path to file containing original data.
# Original data is not publicly available.
original_data_file_path <- ""

# Path to .RData file containing the model fit to original data.
# This file is not publicly available.
original_model_r_data_path <- ""

load(original_model_r_data_path)

gist_original_df <- read_delim(original_data_file_path, delim = ",")
gist_original_df <- gist_original_df |>
  select(AdjTreatm, Size, AgeAtSurg, MitHPF, GenderMale, Rupture, Gastric)

syn_data <- syn(gist_original_df)
gist_sim_df <- syn_data$syn
compare(syn_data, gist_original_df)

gist_sim_df$PatientID <- seq_len(nrow(gist_sim_df))
gist_sim_df <- gist_sim_df |>
  mutate(
    SizeScaled = (Size - mean(Size)) / (2 * sd(Size)),
    AgeAtSurgScaled = (AgeAtSurg - mean(AgeAtSurg)) / (2 * sd(AgeAtSurg)),
    MitHPFScaled = (MitHPF - mean(MitHPF)) / (2 * sd(MitHPF))
  )

gist_sim_long_df <- data.frame(
  PatientID = integer(),
  Size = numeric(),
  SizeScaled = numeric(),
  AgeAtSurg = numeric(),
  AgeAtSurgScaled = numeric(),
  MitHPF = numeric(),
  MitHPFScaled = numeric(),
  GenderMale = numeric(),
  Rupture = numeric(),
  Gastric = numeric(),
  Event = numeric(),
  Time = numeric(),
  AdjOn = numeric(),
  TimeSinceAdjStopped = numeric()
)

for (i in seq_len(nrow(gist_sim_df))) {
  row <- gist_sim_df[i, ]

  for (j in 1:10) {
    new_row <- data.frame(
      PatientID = row$PatientID,
      Size = row$Size,
      SizeScaled = row$SizeScaled,
      AgeAtSurg = row$AgeAtSurg,
      AgeAtSurgScaled = row$AgeAtSurgScaled,
      MitHPF = row$MitHPF,
      MitHPFScaled = row$MitHPFScaled,
      GenderMale = row$GenderMale,
      Rupture = row$Rupture,
      Gastric = row$Gastric,
      Event = 0,
      Time = j,
      AdjOn = if_else(row$AdjTreatm == 1 && j <= 3, 1, 0),
      TimeSinceAdjStopped = if_else(
                                    row$AdjTreatm == 1 && j > 3,
                                    j - 3,
                                    if_else(row$AdjTreatm == 1, 0, j))
    )

    event <- posterior_predict(fit_bernoulli, newdata = new_row, ndraws = 1)
    new_row$Event <- event
    gist_sim_long_df <- rbind(gist_sim_long_df, new_row)

    if (new_row$Event == 1) {
      break
    }
  }
}

gist_sim_short_df <- gist_sim_long_df |>
  group_by(PatientID) |>
  summarise(
    Size = first(Size),
    SizeScaled = first(SizeScaled),
    AgeAtSurg = first(AgeAtSurg),
    AgeAtSurgScaled = first(AgeAtSurgScaled),
    MitHPF = first(MitHPF),
    MitHPFScaled = first(MitHPFScaled),
    GenderMale = first(GenderMale),
    Rupture = first(Rupture),
    Gastric = first(Gastric),
    EventTime = max(Time),
    Censored = 1 - Event[which.max(Time)],
    AdjTreatm = AdjOn[Time == 1],
    .groups = "drop"
  )

gist_sim_long_file_df <- gist_sim_long_df |>
  select(-SizeScaled, -AgeAtSurgScaled, -MitHPFScaled)
write.csv(gist_sim_long_file_df, "./../../data/gist_sim_long.csv",
          row.names = FALSE)

gist_sim_short_file_df <- gist_sim_short_df |>
  select(-SizeScaled, -AgeAtSurgScaled, -MitHPFScaled)
write.csv(gist_sim_short_file_df, "./../../data/gist_sim_short.csv",
          row.names = FALSE)
