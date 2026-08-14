#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
out_dir <- file.path(root, "03_outputs", "16_exploratory_sensitivity_summary")
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

implementations <- data.table(
  sensitivity = c("primary", "unweighted", "interval_1_to_3_years", "exclude_explicit_proxy"),
  directory = c(
    "11_exploratory_upgrade_point",
    "15_exploratory_sensitivity_unweighted",
    "15_exploratory_sensitivity_interval_1_to_3_years",
    "15_exploratory_sensitivity_exclude_explicit_proxy"
  )
)
expected_implementations <- nrow(implementations)

required_files <- unlist(lapply(implementations$directory, function(d) {
  file.path(root, "03_outputs", d, c(
    "low_education_wealth_compensation_point_screen.csv",
    "recovery_phase_model_contrasts.csv",
    "model_convergence_qc.csv"
  ))
}))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing sensitivity outputs: ", paste(missing_files, collapse = "; "))

lowedu_parts <- list()
phase_parts <- list()
qc_parts <- list()
for (i in seq_len(nrow(implementations))) {
  sensitivity_name <- implementations$sensitivity[[i]]
  path <- file.path(root, "03_outputs", implementations$directory[[i]])

  lowedu <- fread(file.path(path, "low_education_wealth_compensation_point_screen.csv"))
  lowedu[, sensitivity := sensitivity_name]
  lowedu_parts[[length(lowedu_parts) + 1L]] <- lowedu

  phase <- fread(file.path(path, "recovery_phase_model_contrasts.csv"))
  phase[, sensitivity := sensitivity_name]
  phase_parts[[length(phase_parts) + 1L]] <- phase

  qc <- fread(file.path(path, "model_convergence_qc.csv"))
  qc[, sensitivity := sensitivity_name]
  qc_parts[[length(qc_parts) + 1L]] <- qc
}

lowedu <- rbindlist(lowedu_parts, fill = TRUE)
phase <- rbindlist(phase_parts, fill = TRUE)
qc <- rbindlist(qc_parts, fill = TRUE)

primary <- lowedu[sensitivity == "primary", .(
  cohort,
  primary_dfle_gap = dfle_gap,
  primary_recovery_relapse = recovery_relapse_contribution_years
)]
lowedu <- merge(lowedu, primary, by = "cohort", all.x = TRUE)
lowedu[, `:=`(
  dfle_gap_direction_adverse = dfle_gap < 0,
  recovery_relapse_direction_adverse = recovery_relapse_contribution_years < 0,
  recovery_relapse_direction_matches_primary =
    sign(recovery_relapse_contribution_years) == sign(primary_recovery_relapse),
  point_magnitude_adverse_ge_0_50y = recovery_relapse_contribution_years <= -0.50
)]
lowedu[, .sensitivity_order := match(sensitivity, implementations$sensitivity)]
setorder(lowedu, cohort, .sensitivity_order)
lowedu[, .sensitivity_order := NULL]
fwrite(lowedu, file.path(out_dir, "low_education_wealth_sensitivity_point_estimates.csv"))

robustness <- lowedu[, .(
  implementations = .N,
  all_dfle_gaps_adverse = all(dfle_gap_direction_adverse),
  all_recovery_relapse_directions_match_primary = all(recovery_relapse_direction_matches_primary),
  all_recovery_relapse_directions_adverse = all(recovery_relapse_direction_adverse),
  minimum_absolute_recovery_relapse_years = min(abs(recovery_relapse_contribution_years)),
  maximum_absolute_recovery_relapse_years = max(abs(recovery_relapse_contribution_years)),
  all_implementations_meet_adverse_0_50y = all(point_magnitude_adverse_ge_0_50y)
), by = cohort]
robustness[, cohort_robustness_gate :=
  implementations == expected_implementations & all_dfle_gaps_adverse &
  all_recovery_relapse_directions_adverse]
setorder(robustness, cohort)
fwrite(robustness, file.path(out_dir, "low_education_wealth_sensitivity_robustness.csv"))

phase[, contrast_direction := fifelse(
  contrast %in% c("low_vs_high_early_recovery", "low_vs_high_sustained_recovery"),
  fifelse(hazard_ratio > 1, "low_wealth_more_relapse", "low_wealth_less_relapse"),
  fifelse(hazard_ratio > 1, "stronger_low_wealth_gradient_after_sustained_recovery", "weaker_low_wealth_gradient_after_sustained_recovery")
)]
phase[, .sensitivity_order := match(sensitivity, implementations$sensitivity)]
setorder(phase, cohort, contrast, .sensitivity_order)
phase[, .sensitivity_order := NULL]
fwrite(phase, file.path(out_dir, "recovery_phase_sensitivity_point_estimates.csv"))

phase_robustness <- phase[
  contrast %in% c("low_vs_high_early_recovery", "low_vs_high_sustained_recovery"),
  .(
    implementations = .N,
    all_hazard_ratios_above_one = all(hazard_ratio > 1),
    minimum_hazard_ratio = min(hazard_ratio),
    maximum_hazard_ratio = max(hazard_ratio)
  ),
  by = .(cohort, contrast)
]
fwrite(phase_robustness, file.path(out_dir, "recovery_phase_sensitivity_robustness.csv"))

qc_summary <- qc[, .(
  models = .N,
  all_converged = all(converged),
  all_coefficients_finite = all(is.finite(coefficient_abs_max)),
  maximum_abs_coefficient = max(coefficient_abs_max),
  warning_models = sum(warning_n > 0)
), by = .(sensitivity, cohort)]
qc_summary[, .sensitivity_order := match(sensitivity, implementations$sensitivity)]
setorder(qc_summary, .sensitivity_order, cohort)
qc_summary[, .sensitivity_order := NULL]
fwrite(qc_summary, file.path(out_dir, "exploratory_sensitivity_model_qc.csv"))

decision <- data.table(
  complete_four_implementations = all(robustness$implementations == expected_implementations),
  cohorts_with_directional_robustness = robustness[cohort_robustness_gate == TRUE, .N],
  cohorts_with_full_magnitude_robustness = robustness[all_implementations_meet_adverse_0_50y == TRUE, .N],
  low_education_wealth_directionally_robust_in_at_least_two_cohorts =
    robustness[cohort_robustness_gate == TRUE, .N] >= 2L,
  all_sensitivity_models_converged = all(qc_summary$all_converged & qc_summary$all_coefficients_finite)
)
fwrite(decision, file.path(out_dir, "exploratory_sensitivity_decision.csv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, "sessionInfo_exploratory_sensitivity_summary.txt")
)
print(robustness)
print(decision)
