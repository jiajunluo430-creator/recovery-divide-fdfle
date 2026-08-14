#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
outputs <- file.path(revision_root, "03_outputs")
out_dir <- file.path(outputs, "15_consolidated_evidence")
log_dir <- file.path(revision_root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("53_consolidate_revision_evidence_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

cohorts <- c("CHARLS", "ELSA", "HRS", "MHAS")
point_dirs <- c(
  primary = "02_revision_point",
  adl_only = "02_revision_point_adl_only",
  at_least_two = "02_revision_point_at_least_two",
  incident_inception = "02_revision_point_incident_inception",
  entry_age_60_64 = "02_revision_point_entry_age_60_64",
  pre_covid = "02_revision_point_pre_covid",
  common_interval_1_5_3_25 = "02_revision_point_common_interval",
  two_year_interval_1_75_2_25 = "02_revision_point_two_year_interval",
  common_calendar = "02_revision_point_common_calendar",
  mortality_common_waves = "02_revision_point_mortality_common_waves",
  piecewise_exponential = "02_revision_point_piecewise_exponential",
  age110_hazards_frozen_95 = "02_revision_point_age110_tail",
  age110_hazards_frozen_100 = "02_revision_point_age100_hazard",
  wealth_not_after_first_difficulty = "02_revision_point_wealth_pre_difficulty",
  unweighted = "02_revision_point_unweighted",
  self_report_only_where_observed = "02_revision_point_self_report_only"
)

read_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  fread(path)
}

point_parts <- list()
for (sensitivity_name in names(point_dirs)) {
  d <- file.path(outputs, point_dirs[[sensitivity_name]])
  gap <- read_if_exists(file.path(d, "low_high_fdfle_gaps.csv"))
  shapley <- read_if_exists(file.path(d, "shapley_decomposition_point.csv"))
  life <- read_if_exists(file.path(d, "life_expectancy_point_estimates.csv"))
  if (is.null(gap) || is.null(shapley)) next
  g <- gap[
    module == "wealth_within_low_education" & estimand == "population_initialised",
    .(cohort, gap_years = low_minus_high_gap)
  ]
  rr <- shapley[
    module == "wealth_within_low_education" & estimand == "population_initialised" &
      block %chin% c("recovery", "relapse"),
    .(recovery_relapse_years = sum(contribution_years)),
    by = cohort
  ]
  comp <- dcast(
    shapley[
      module == "wealth_within_low_education" & estimand == "population_initialised",
      .(cohort, block, contribution_years)
    ],
    cohort ~ block,
    value.var = "contribution_years"
  )
  z <- Reduce(function(x, y) merge(x, y, by = "cohort", all = TRUE), list(g, rr, comp))
  if (!is.null(life)) {
    levels <- dcast(
      life[
        module == "wealth_within_low_education" & estimand == "population_initialised" &
          ses %chin% c("high", "low") &
          metric %chin% c("total_life_expectancy", "fdfle", "difficulty_years", "residual_alive_at_end"),
        .(cohort, key = paste(ses, metric, sep = "_"), estimate)
      ],
      cohort ~ key,
      value.var = "estimate"
    )
    z <- merge(z, levels, by = "cohort", all.x = TRUE)
  }
  z[, sensitivity := sensitivity_name]
  point_parts[[sensitivity_name]] <- z
}
focused_sensitivity <- rbindlist(point_parts, use.names = TRUE, fill = TRUE)
setcolorder(focused_sensitivity, c("sensitivity", "cohort", setdiff(names(focused_sensitivity), c("sensitivity", "cohort"))))
focused_sensitivity[, `:=`(
  sensitivity_order = match(sensitivity, names(point_dirs)),
  cohort_order = match(cohort, cohorts)
)]
setorder(focused_sensitivity, sensitivity_order, cohort_order)
focused_sensitivity[, c("sensitivity_order", "cohort_order") := NULL]

primary_dir <- file.path(outputs, "02_revision_point")
primary_gap <- fread(file.path(primary_dir, "low_high_fdfle_gaps.csv"))[
  estimand == "population_initialised"
]
primary_shapley <- fread(file.path(primary_dir, "shapley_decomposition_point.csv"))[
  estimand == "population_initialised"
]
primary_life <- fread(file.path(primary_dir, "life_expectancy_point_estimates.csv"))[
  estimand == "population_initialised" & ses %chin% c("high", "low")
]
primary_rr <- primary_shapley[block %chin% c("recovery", "relapse"), .(
  recovery_relapse_years = sum(contribution_years)
), by = .(cohort, module)]
primary_point <- primary_gap[, .(
  cohort, module, high_fdfle, low_fdfle, gap_years = low_minus_high_gap
)][primary_rr, on = .(cohort, module)]

inference_dirs <- c(
  primary_education = "04_primary_education_inference",
  primary_wealth = "04_primary_wealth_inference",
  wealth_within_low_education = "04_lowedu_inference",
  focused_adl_only = "04_lowedu_inference_adl_only",
  focused_at_least_two = "04_lowedu_inference_at_least_two"
)
inference_parts <- list()
for (module_name in names(inference_dirs)) {
  path <- file.path(outputs, inference_dirs[[module_name]], "bootstrap_percentile_ci_summary.csv")
  z <- read_if_exists(path)
  if (is.null(z)) next
  z[, analysis := module_name]
  inference_parts[[module_name]] <- z
}
bootstrap_inference <- rbindlist(inference_parts, use.names = TRUE, fill = TRUE)
if (nrow(bootstrap_inference)) {
  setcolorder(bootstrap_inference, c("analysis", "cohort", "metric", setdiff(names(bootstrap_inference), c("analysis", "cohort", "metric"))))
}

confirmed <- fread(file.path(outputs, "05_confirmed_state_sensitivity", "confirmed_state_cluster_robust_estimates.csv"))
msm <- rbindlist(lapply(tolower(cohorts), function(x) {
  fread(file.path(outputs, "07_continuous_time_panel_msm", paste0("process_intensity_ratios_", x, ".csv")))
}), use.names = TRUE, fill = TRUE)
enriched <- rbindlist(lapply(tolower(cohorts), function(x) {
  fread(file.path(outputs, "08_enriched_covariate_sensitivity", paste0("enriched_transition_results_", x, ".csv")))
}), use.names = TRUE, fill = TRUE)
ipcw <- rbindlist(lapply(tolower(cohorts), function(x) {
  fread(file.path(outputs, "09_nondeath_observation_ipcw", paste0("ipcw_transition_results_", x, ".csv")))
}), use.names = TRUE, fill = TRUE)

flow_dir <- file.path(outputs, "06_flow_missingness_support")
focused_dir <- file.path(outputs, "11_focused_baseline_support")
calibration_dir <- file.path(outputs, "10_internal_survival_calibration")
riskset_dir <- file.path(outputs, "01_revision_risksets")

key_file_ledger <- data.table(
  manuscript_domain = c(
    "Primary point estimates", "Binding focused uncertainty", "Primary education uncertainty",
    "Primary wealth uncertainty", "ADL-only uncertainty", "At-least-two uncertainty",
    "Functional construct mapping", "Sample flow and missingness", "Focused baseline and support",
    "Continuous-time panel model", "Confirmed states", "Enriched covariates", "Non-death IPCW",
    "Internal TLE calibration", "Wealth timing and interval audit"
  ),
  source = c(
    file.path(primary_dir, "low_high_fdfle_gaps.csv"),
    file.path(outputs, "04_lowedu_inference", "bootstrap_percentile_ci_summary.csv"),
    file.path(outputs, "04_primary_education_inference", "bootstrap_percentile_ci_summary.csv"),
    file.path(outputs, "04_primary_wealth_inference", "bootstrap_percentile_ci_summary.csv"),
    file.path(outputs, "04_lowedu_inference_adl_only", "bootstrap_percentile_ci_summary.csv"),
    file.path(outputs, "04_lowedu_inference_at_least_two", "bootstrap_percentile_ci_summary.csv"),
    file.path(outputs, "13_data_dictionary", "functional_item_mapping_by_cohort_wave.csv"),
    file.path(flow_dir, "sample_flow_by_module_v1_1.csv"),
    file.path(focused_dir, "focused_loweducation_wealth_transition_support.csv"),
    file.path(outputs, "07_continuous_time_panel_msm"),
    file.path(outputs, "05_confirmed_state_sensitivity", "confirmed_state_cluster_robust_estimates.csv"),
    file.path(outputs, "08_enriched_covariate_sensitivity"),
    file.path(outputs, "09_nondeath_observation_ipcw"),
    file.path(calibration_dir, "model_vs_weighted_km_tle_calibration.csv"),
    file.path(riskset_dir, "wealth_measurement_timing_audit.csv")
  )
)
key_file_ledger[, exists := file.exists(source) | dir.exists(source)]

fwrite(focused_sensitivity, file.path(out_dir, "focused_point_sensitivity_matrix.csv"))
fwrite(primary_point, file.path(out_dir, "primary_education_wealth_point_summary.csv"))
fwrite(primary_life, file.path(out_dir, "primary_absolute_levels_long.csv"))
fwrite(primary_shapley, file.path(out_dir, "primary_six_block_decomposition_long.csv"))
fwrite(bootstrap_inference, file.path(out_dir, "available_household_bootstrap_inference_long.csv"))
fwrite(confirmed, file.path(out_dir, "confirmed_state_transition_ratios.csv"))
fwrite(msm, file.path(out_dir, "continuous_time_process_intensity_ratios.csv"))
fwrite(enriched, file.path(out_dir, "enriched_covariate_transition_ratios.csv"))
fwrite(ipcw, file.path(out_dir, "nondeath_ipcw_transition_ratios.csv"))
fwrite(key_file_ledger, file.path(out_dir, "manuscript_evidence_file_ledger.csv"))

focused_main <- bootstrap_inference[
  analysis == "wealth_within_low_education" &
    metric %chin% c("population_fdfle_gap", "population_recovery_relapse"),
  .(cohort, metric, point_estimate, ci_low, ci_high, valid_replicates)
]
primary_main <- bootstrap_inference[
  analysis %chin% c("primary_education", "primary_wealth") &
    metric %chin% c("population_fdfle_gap", "population_recovery_relapse"),
  .(analysis, cohort, metric, point_estimate, ci_low, ci_high, valid_replicates)
]
report <- c(
  "# Consolidated revision evidence snapshot",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Focused secondary binding household-bootstrap endpoints",
  "",
  paste(capture.output(print(focused_main)), collapse = "\n"),
  "",
  "## Available primary education/wealth household-bootstrap endpoints",
  "",
  if (nrow(primary_main)) paste(capture.output(print(primary_main)), collapse = "\n") else "Pending.",
  "",
  "## Continuous-time panel intensity ratios",
  "",
  paste(capture.output(print(msm)), collapse = "\n"),
  "",
  "## Missing evidence-file checks",
  "",
  paste(capture.output(print(key_file_ledger[exists == FALSE])), collapse = "\n")
)
writeLines(report, file.path(out_dir, "consolidated_revision_evidence_snapshot.md"), useBytes = TRUE)

cat(paste(report, collapse = "\n"), "\n")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_53_consolidated_evidence.txt"))
