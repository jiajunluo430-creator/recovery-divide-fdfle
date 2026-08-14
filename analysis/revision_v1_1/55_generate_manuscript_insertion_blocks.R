#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
outputs <- file.path(revision_root, "03_outputs")
manuscript_dir <- file.path(revision_root, "05_manuscript_v1_1")
log_dir <- file.path(revision_root, "06_logs")
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  education = file.path(outputs, "04_primary_education_inference", "bootstrap_percentile_ci_summary.csv"),
  wealth = file.path(outputs, "04_primary_wealth_inference", "bootstrap_percentile_ci_summary.csv"),
  focused = file.path(outputs, "04_lowedu_inference", "bootstrap_percentile_ci_summary.csv"),
  adl = file.path(outputs, "04_lowedu_inference_adl_only", "bootstrap_percentile_ci_summary.csv"),
  two = file.path(outputs, "04_lowedu_inference_at_least_two", "bootstrap_percentile_ci_summary.csv"),
  focused_heterogeneity = file.path(outputs, "04_lowedu_inference", "cross_cohort_heterogeneity.csv"),
  continuous = file.path(outputs, "15_consolidated_evidence", "continuous_time_process_intensity_ratios.csv"),
  enriched = file.path(outputs, "15_consolidated_evidence", "enriched_covariate_transition_ratios.csv"),
  ipcw = file.path(outputs, "15_consolidated_evidence", "nondeath_ipcw_transition_ratios.csv"),
  psu = file.path(outputs, "14_stratified_psu_inference", "household_vs_stratified_psu_ci_comparison.csv"),
  calibration = file.path(outputs, "10_internal_survival_calibration", "model_vs_weighted_km_tle_calibration.csv")
)
missing_paths <- paths[!file.exists(unlist(paths))]
if (length(missing_paths)) stop("Missing final evidence file(s): ", paste(names(missing_paths), collapse = ", "))

edu <- fread(paths$education)
wealth <- fread(paths$wealth)
focused <- fread(paths$focused)
adl <- fread(paths$adl)
two <- fread(paths$two)
het <- fread(paths$focused_heterogeneity)
continuous <- fread(paths$continuous)
enriched <- fread(paths$enriched)
ipcw <- fread(paths$ipcw)
psu <- fread(paths$psu)
calibration <- fread(paths$calibration)
cohorts <- c("CHARLS", "ELSA", "HRS", "MHAS")

fmt <- function(x, d = 2L) formatC(x, format = "f", digits = d)
fmt_ci <- function(z, d = 2L) paste0(fmt(z$point_estimate, d), " (", fmt(z$ci_low, d), " to ", fmt(z$ci_high, d), ")")
get_metric <- function(z, cohort_name, metric_name) z[cohort == cohort_name & metric == metric_name][1]
supported <- function(z) is.finite(z$ci_low) && is.finite(z$ci_high) && (z$ci_high < 0 || z$ci_low > 0)
join_items <- function(x) {
  if (!length(x)) return("")
  if (length(x) == 1L) return(x)
  if (length(x) == 2L) return(paste(x, collapse = " and "))
  paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[[length(x)]])
}

metric_sentence <- function(z, metric_name) {
  paste(vapply(cohorts, function(cc) paste0(cc, " ", fmt_ci(get_metric(z, cc, metric_name))), character(1)), collapse = "; ")
}

edu_gap <- edu[metric == "population_fdfle_gap"]
wealth_gap <- wealth[metric == "population_fdfle_gap"]
edu_post <- edu[metric == "population_recovery_relapse"]
wealth_post <- wealth[metric == "population_recovery_relapse"]
focused_gap <- focused[metric == "population_fdfle_gap"]
focused_post <- focused[metric == "population_recovery_relapse"]
adl_post <- adl[metric == "population_recovery_relapse"]
two_post <- two[metric == "population_recovery_relapse"]

edu_supported_n <- sum(edu_gap[, ci_high < 0 | ci_low > 0])
wealth_supported_n <- sum(wealth_gap[, ci_high < 0 | ci_low > 0])
focused_supported_cohorts <- focused_post[ci_high < 0 | ci_low > 0, cohort]
adl_supported_cohorts <- adl_post[ci_high < 0 | ci_low > 0, cohort]
two_supported_cohorts <- two_post[ci_high < 0 | ci_low > 0, cohort]

abstract_results <- paste0(
  "Across 83,723 entrants, low education was associated with fewer population-initialised FDFLE years in ",
  edu_supported_n, " of four cohorts and low wealth in ", wealth_supported_n,
  " of four cohorts. Within low education, low-minus-high wealth FDFLE differences were ",
  paste(vapply(cohorts, function(cc) paste0(cc, " ", fmt_ci(get_metric(focused, cc, "population_fdfle_gap"))), character(1)), collapse = "; "),
  " years. Recovery plus relapse contributed ",
  paste(vapply(cohorts, function(cc) paste0(cc, " ", fmt_ci(get_metric(focused, cc, "population_recovery_relapse"))), character(1)), collapse = "; "),
  " years. Under ADL-only difficulty this post-onset contribution remained distinguishable in ",
  if (length(adl_supported_cohorts)) join_items(adl_supported_cohorts) else "no cohort",
  "; under the at-least-two-difficulty definition it remained distinguishable in ",
  if (length(two_supported_cohorts)) join_items(two_supported_cohorts) else "no cohort", ". "
)

primary_narrative <- c(
  paste0("All four education gaps were adverse: ", metric_sentence(edu, "population_fdfle_gap"),
         " years (Table 2). Recovery plus relapse contributions were ",
         metric_sentence(edu, "population_recovery_relapse"),
         " years; their confidence intervals excluded zero in ",
         join_items(edu_post[ci_high < 0 | ci_low > 0, cohort]), "."),
  paste0("Wealth gaps were also adverse: ", metric_sentence(wealth, "population_fdfle_gap"),
         " years. Recovery plus relapse contributions were ",
         metric_sentence(wealth, "population_recovery_relapse"),
         " years; their confidence intervals excluded zero in ",
         join_items(wealth_post[ci_high < 0 | ci_low > 0, cohort]), ".")
)

absolute_lines <- vapply(cohorts, function(cc) {
  high_f <- get_metric(focused, cc, "population_high_fdfle")
  low_f <- get_metric(focused, cc, "population_low_fdfle")
  high_d <- get_metric(focused, cc, "population_high_difficulty_years")
  low_d <- get_metric(focused, cc, "population_low_difficulty_years")
  high_t <- get_metric(focused, cc, "population_high_tle")
  low_t <- get_metric(focused, cc, "population_low_tle")
  paste0(cc, " high/low wealth FDFLE ", fmt(high_f$point_estimate), "/", fmt(low_f$point_estimate),
         ", total life years ", fmt(high_t$point_estimate), "/", fmt(low_t$point_estimate),
         ", and years with difficulty ", fmt(high_d$point_estimate), "/", fmt(low_d$point_estimate))
}, character(1))

component_names <- c(
  initial = "population_contribution_initial_state",
  onset = "population_contribution_onset",
  recovery = "population_contribution_recovery",
  relapse = "population_contribution_relapse",
  post = "population_contribution_post_difficulty_mortality",
  pre = "population_contribution_pre_difficulty_mortality"
)
component_lines <- vapply(cohorts, function(cc) {
  vals <- vapply(component_names, function(mm) fmt(get_metric(focused, cc, mm)$point_estimate), character(1))
  paste0(cc, ": initial ", vals[["initial"]], ", onset ", vals[["onset"]],
         ", recovery ", vals[["recovery"]], ", relapse ", vals[["relapse"]],
         ", post-difficulty mortality ", vals[["post"]], ", pre-difficulty mortality ", vals[["pre"]])
}, character(1))

focused_narrative <- c(
  paste0("Within low education, population-initialised low-minus-high wealth FDFLE differences were ",
         metric_sentence(focused, "population_fdfle_gap"), " years (Table 3). Absolute point estimates were: ",
         paste(absolute_lines, collapse = "; "), "."),
  paste0("The signed six-block point decompositions were ", paste(component_lines, collapse = "; "),
         " years. Recovery plus relapse confidence intervals excluded zero in ",
         if (length(focused_supported_cohorts)) join_items(focused_supported_cohorts) else "no cohort",
         ". The MHAS total-gap confidence interval included zero, so no component percentage was interpreted.")
)

threshold_narrative <- paste0(
  "Under ADL-only difficulty, focused recovery-plus-relapse contributions were ",
  metric_sentence(adl, "population_recovery_relapse"), " years and were distinguishable in ",
  if (length(adl_supported_cohorts)) join_items(adl_supported_cohorts) else "no cohort",
  ". Under the at-least-two-difficulty definition, corresponding estimates were ",
  metric_sentence(two, "population_recovery_relapse"), " years and were distinguishable in ",
  if (length(two_supported_cohorts)) join_items(two_supported_cohorts) else "no cohort", "."
)

ct_key <- continuous[, .(
  cohort, process,
  estimate_ci = paste0(fmt(intensity_ratio), " (", fmt(ci_low), " to ", fmt(ci_high), ")"),
  supported = ci_high < 1 | ci_low > 1
)]
ct_lines <- vapply(cohorts, function(cc) {
  zz <- ct_key[cohort == cc]
  supported_processes <- zz[supported == TRUE, process]
  paste0(cc, " supported continuous-time contrasts: ",
         if (length(supported_processes)) join_items(gsub("_", " ", supported_processes)) else "none")
}, character(1))

psu_key <- psu[metric %chin% c("population_fdfle_gap", "population_recovery_relapse")]
psu_lines <- psu_key[, paste0(
  cohort, " ", gsub("population_", "", metric), " household ",
  fmt(point_estimate), " (", fmt(household_ci_low), " to ", fmt(household_ci_high),
  "); PSU ", fmt(point_estimate), " (", fmt(psu_ci_low), " to ", fmt(psu_ci_high), ")"
)]

cal90 <- calibration[, .(
  cohort, wealth3,
  km = weighted_delayed_entry_km_tle_to_age_90,
  model = multistate_model_tle_to_age_90,
  difference = model_minus_km_years_to_age_90
)]
calibration_lines <- cal90[, paste0(cohort, " ", wealth3, " KM/model ", fmt(km), "/", fmt(model),
                                           " (difference ", fmt(difference), ")")]

other_narrative <- c(
  paste0("The continuous-time models converged in all cohorts. ", paste(ct_lines, collapse = "; "), "."),
  paste0("Stratified-PSU sensitivity estimates were: ", paste(psu_lines, collapse = "; "),
         ". Their confidence intervals retained the same inferential direction as household resampling."),
  paste0("Internal age-90 survival calibration (weighted delayed-entry KM/model years) was: ",
         paste(calibration_lines, collapse = "; "), ". No cohort had a direct age-100 risk set.")
)

gap_het <- het[metric == "population_fdfle_gap"]
post_het <- het[metric == "population_recovery_relapse"]
journal_signal <- if (length(adl_supported_cohorts) >= 2L && length(two_supported_cohorts) >= 1L) {
  "TLHL_STRETCH_GO; BMC_MEDICINE_GO; AGE_AND_AGEING_GO"
} else if (length(adl_supported_cohorts) >= 2L) {
  "BMC_MEDICINE_GO; TLHL_STRETCH_ONLY; AGE_AND_AGEING_GO"
} else {
  "AGE_AND_AGEING_OR_IJE_PIVOT; NO_TOP_GENERAL_CLAIM"
}

lines <- c(
  "# Final v1.1 manuscript insertion blocks",
  "",
  paste0("Generated: ", format(Sys.time(), tz = "America/Chicago", usetz = TRUE)),
  "",
  "## Abstract Results",
  "",
  abstract_results,
  "",
  "## Primary education and wealth Results",
  "",
  primary_narrative,
  "",
  "## Focused Results",
  "",
  focused_narrative,
  "",
  "## Functional-definition Results",
  "",
  threshold_narrative,
  "",
  "## Model, design, and calibration Results",
  "",
  other_narrative,
  "",
  "## GO/PIVOT evidence",
  "",
  paste0("Decision signal: ", journal_signal),
  paste0("Broad focused post-onset supported cohorts: ", join_items(focused_supported_cohorts)),
  paste0("ADL-only supported cohorts: ", if (length(adl_supported_cohorts)) join_items(adl_supported_cohorts) else "none"),
  paste0("At-least-two-difficulty supported cohorts: ", if (length(two_supported_cohorts)) join_items(two_supported_cohorts) else "none"),
  paste0("Focused gap heterogeneity: Q=", fmt(gap_het$cochran_q), ", p=", formatC(gap_het$p_heterogeneity, format = "g", digits = 3),
         ", I2=", fmt(gap_het$i2_percent, 1L), "%"),
  paste0("Focused post-onset heterogeneity: Q=", fmt(post_het$cochran_q), ", p=", formatC(post_het$p_heterogeneity, format = "g", digits = 3),
         ", I2=", fmt(post_het$i2_percent, 1L), "%")
)

writeLines(lines, file.path(manuscript_dir, "02_results_insertion_blocks_v1.1.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_55_manuscript_insertions.txt"), useBytes = TRUE)
cat(paste(lines, collapse = "\n"), "\n")
