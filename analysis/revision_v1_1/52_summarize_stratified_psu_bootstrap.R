#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
boot_dir <- file.path(revision_root, "03_outputs", "03_lowedu_stratified_psu_bootstrap_final")
household_dir <- file.path(revision_root, "03_outputs", "04_lowedu_inference")
out_dir <- file.path(revision_root, "03_outputs", "14_stratified_psu_inference")
log_dir <- file.path(revision_root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- c("ELSA", "HRS")
target <- as.integer(Sys.getenv("D4_REV_TARGET_BOOTSTRAP", unset = "500"))
qfun <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE, type = 6))

qc <- rbindlist(lapply(cohorts, function(cohort_name) {
  fread(file.path(boot_dir, paste0("bootstrap_", tolower(cohort_name), "_qc.csv")))
}), fill = TRUE)
metrics <- rbindlist(lapply(cohorts, function(cohort_name) {
  fread(file.path(boot_dir, paste0("bootstrap_", tolower(cohort_name), "_metrics.csv")))
}), fill = TRUE)
qc[, status := tolower(status)]
valid <- qc[status == "valid", .(cohort, replicate)]
setorder(valid, cohort, replicate)
valid[, inference_order := seq_len(.N), by = cohort]
valid_use <- valid[inference_order <= target, .(cohort, replicate)]
valid_qc <- qc[, .(
  attempts = .N,
  valid = sum(status == "valid"),
  failed = sum(status != "valid"),
  design_missing_people_excluded = max(design_missing_people_excluded, na.rm = TRUE),
  sampled_clusters_median = median(sampled_clusters_unique, na.rm = TRUE)
), by = cohort]
if (any(valid_qc$valid < target)) {
  stop("Fewer than target valid PSU replicates: ", paste(valid_qc[, paste0(cohort, "=", valid)], collapse = "; "))
}

metrics <- metrics[valid_use, on = .(cohort, replicate), nomatch = 0]
summary_ci <- metrics[, .(
  psu_bootstrap_mean = mean(estimate),
  psu_bootstrap_sd = sd(estimate),
  psu_ci_low = qfun(estimate, 0.025),
  psu_ci_high = qfun(estimate, 0.975)
), by = .(cohort, metric)]

household <- fread(file.path(household_dir, "bootstrap_percentile_ci_summary.csv"))[
  cohort %chin% cohorts,
  .(
    cohort, metric, point_estimate,
    household_ci_low = ci_low,
    household_ci_high = ci_high,
    household_bootstrap_sd = bootstrap_sd
  )
]
comparison <- household[summary_ci, on = .(cohort, metric)]
comparison[, `:=`(
  household_ci_width = household_ci_high - household_ci_low,
  psu_ci_width = psu_ci_high - psu_ci_low,
  psu_to_household_width_ratio = (psu_ci_high - psu_ci_low) / (household_ci_high - household_ci_low)
)]

key_metrics <- c(
  "population_fdfle_gap", "population_recovery_relapse",
  "population_contribution_initial_state", "population_contribution_onset",
  "population_contribution_recovery", "population_contribution_relapse",
  "population_contribution_post_difficulty_mortality",
  "population_contribution_pre_difficulty_mortality"
)
stability <- metrics[metric %chin% key_metrics, {
  reps <- sort(unique(replicate))
  half <- reps[seq_len(floor(length(reps) / 2L))]
  full_ci <- c(qfun(estimate, 0.025), qfun(estimate, 0.975))
  half_ci <- c(qfun(estimate[replicate %in% half], 0.025), qfun(estimate[replicate %in% half], 0.975))
  .(
    half_replicates = length(half),
    full_replicates = length(reps),
    half_ci_low = half_ci[[1L]], half_ci_high = half_ci[[2L]],
    full_ci_low = full_ci[[1L]], full_ci_high = full_ci[[2L]],
    max_endpoint_drift_years = max(abs(half_ci - full_ci))
  )
}, by = .(cohort, metric)]

fwrite(valid_qc, file.path(out_dir, "stratified_psu_bootstrap_qc.csv"))
fwrite(summary_ci, file.path(out_dir, "stratified_psu_percentile_ci.csv"))
fwrite(comparison, file.path(out_dir, "household_vs_stratified_psu_ci_comparison.csv"))
fwrite(stability, file.path(out_dir, "stratified_psu_monte_carlo_stability.csv"))

cat("\nHousehold versus stratified-PSU inference for key endpoints:\n")
print(comparison[metric %chin% c("population_fdfle_gap", "population_recovery_relapse")])
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_52_stratified_psu_summary.txt"))

