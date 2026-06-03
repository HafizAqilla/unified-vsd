#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

workbook_path <- file.path(
  "outputs",
  "pre_surgery_patient_workbook",
  "MasterSheet_Hasil.xlsx"
)

out_dir <- file.path("outputs", "pre_surgery_patient_workbook", "r_stats")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rmse_raw <- read_excel(workbook_path, sheet = "RMSE Summary", skip = 3)

rmse <- rmse_raw %>%
  filter(!is.na(`Patient ID`)) %>%
  transmute(
    patient_id = `Patient ID`,
    calibration_status = `Calibration Status`,
    target_tier_profile = `Target Tier Profile`,
    inclusion_flag_summary = `Inclusion Flag Summary`,
    baseline_primary_rmse = as.numeric(`Baseline Primary RMSE`),
    best_primary_rmse = as.numeric(`Best/Calibrated Primary RMSE`),
    rmse_improvement_abs = as.numeric(`RMSE Improvement Abs`),
    rmse_improvement_pct = as.numeric(`RMSE Improvement %`),
    baseline_full_rmse = as.numeric(`Baseline Full RMSE`),
    best_full_rmse = as.numeric(`Best/Calibrated Full RMSE`),
    worst_baseline_metric = `Worst Baseline Error Metric`,
    worst_baseline_error_pct = as.numeric(`Worst Baseline Error %`),
    worst_best_metric = `Worst Best-Model Error Metric`,
    worst_best_error_pct = as.numeric(`Worst Best-Model Error %`),
    pass_fail_review_flag = `Pass/Fail Review Flag`,
    clinical_notes = `Clinical Consistency / Notes`
  )

summarise_numeric <- function(x) {
  c(
    n = sum(!is.na(x)),
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

numeric_summary_row <- function(metric_name, x) {
  stats <- as.list(summarise_numeric(x))
  tibble(metric = metric_name, !!!stats)
}

cohort_summary <- bind_rows(
  numeric_summary_row("baseline_primary_rmse", rmse$baseline_primary_rmse),
  numeric_summary_row("best_primary_rmse", rmse$best_primary_rmse),
  numeric_summary_row("rmse_improvement_pct", rmse$rmse_improvement_pct),
  numeric_summary_row("baseline_full_rmse", rmse$baseline_full_rmse),
  numeric_summary_row("best_full_rmse", rmse$best_full_rmse)
)

status_counts <- rmse %>%
  count(calibration_status, name = "n") %>%
  arrange(desc(n))

flag_counts <- rmse %>%
  count(pass_fail_review_flag, name = "n") %>%
  arrange(desc(n))

patient_long <- rmse %>%
  select(patient_id, baseline_primary_rmse, best_primary_rmse) %>%
  pivot_longer(
    cols = c(baseline_primary_rmse, best_primary_rmse),
    names_to = "stage",
    values_to = "primary_rmse"
  ) %>%
  mutate(
    stage = recode(
      stage,
      baseline_primary_rmse = "Baseline",
      best_primary_rmse = "Best/Calibrated"
    )
  )

write_csv(rmse, file.path(out_dir, "rmse_summary_table.csv"))
write_csv(cohort_summary, file.path(out_dir, "cohort_numeric_summary.csv"))
write_csv(status_counts, file.path(out_dir, "calibration_status_counts.csv"))
write_csv(flag_counts, file.path(out_dir, "review_flag_counts.csv"))

plot_primary <- ggplot(patient_long, aes(x = reorder(patient_id, primary_rmse), y = primary_rmse, fill = stage)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Baseline" = "#9aa5b1", "Best/Calibrated" = "#1f77b4")) +
  labs(
    title = "Pre-Surgery Primary RMSE by Patient",
    x = NULL,
    y = "Primary governed RMSE",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "primary_rmse_by_patient.png"),
  plot = plot_primary,
  width = 9,
  height = 5.5,
  dpi = 300
)

cat("Workbook:", workbook_path, "\n\n")
cat("Patient count:", nrow(rmse), "\n\n")
cat("Calibration status counts:\n")
print(status_counts)
cat("\nReview flag counts:\n")
print(flag_counts)
cat("\nNumeric summary:\n")
print(cohort_summary)
cat("\nFiles written to:", out_dir, "\n")
