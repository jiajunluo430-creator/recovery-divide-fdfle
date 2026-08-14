#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
outputs <- file.path(revision_root, "03_outputs")
out_dir <- file.path(revision_root, "04_tables_v1_1")
supp_dir <- file.path(out_dir, "supplementary")
source_dir <- file.path(out_dir, "source_data")
log_dir <- file.path(revision_root, "06_logs")
for (d in c(out_dir, supp_dir, source_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  education = file.path(outputs, "04_primary_education_inference", "bootstrap_percentile_ci_summary.csv"),
  wealth = file.path(outputs, "04_primary_wealth_inference", "bootstrap_percentile_ci_summary.csv"),
  focused = file.path(outputs, "04_lowedu_inference", "bootstrap_percentile_ci_summary.csv"),
  adl = file.path(outputs, "04_lowedu_inference_adl_only", "bootstrap_percentile_ci_summary.csv"),
  two = file.path(outputs, "04_lowedu_inference_at_least_two", "bootstrap_percentile_ci_summary.csv")
)
missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required)) {
  stop("Required final inference files are missing: ", paste(names(missing_required), missing_required, collapse = " | "))
}

fmt_num <- function(x, digits = 2L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}
fmt_ci <- function(est, low, high, digits = 2L) {
  paste0(fmt_num(est, digits), " (", fmt_num(low, digits), " to ", fmt_num(high, digits), ")")
}

baseline <- fread(file.path(outputs, "06_flow_missingness_support", "baseline_characteristics_corrected_v1_1.csv"))
gate <- fread(file.path(root, "03_outputs", "01_gate0", "gate0_event_counts_by_cohort.csv"))
flow <- fread(file.path(outputs, "06_flow_missingness_support", "sample_flow_by_module_v1_1.csv"))
death_loss <- fread(file.path(outputs, "06_flow_missingness_support", "death_vs_non_death_loss_summary_v1_1.csv"))

flow_final <- flow[
  , .SD[stage_order == max(stage_order)],
  by = .(cohort, module)
]
flow_wide <- dcast(
  flow_final[
    module %chin% c(
      "primary_education",
      "primary_wealth",
      "focused_low_education_wealth"
    ),
    .(cohort, key = fifelse(
      module == "primary_education", "primary_education_people",
      fifelse(module == "primary_wealth", "primary_wealth_people", "focused_lowedu_wealth_people")
    ), people)
  ],
  cohort ~ key,
  value.var = "people"
)
gate_for_merge <- gate[, setdiff(names(gate), "country"), with = FALSE]
table1_numeric <- Reduce(function(x, y) merge(x, y, by = "cohort", all = TRUE), list(
  baseline, gate_for_merge, flow_wide,
  death_loss[, .(cohort, unknown_vital_or_observation_intervals, primary_supported_deaths)]
))
table1_display <- table1_numeric[, .(
  Cohort = cohort,
  Country = country,
  `Eligible entrants` = eligible_entrants,
  `Mean age, years (SD)` = paste0(fmt_num(age_mean, 1), " (", fmt_num(age_sd, 1), ")"),
  `Women among known sex, n/N (%)` = paste0(women_n, "/", sex_known_n, " (", fmt_num(women_percent_among_known_sex, 1), ")"),
  `Entry functional difficulty, n (%)` = paste0(baseline_difficulty_n, " (", fmt_num(baseline_difficulty_percent, 1), ")"),
  `Primary education analysis, people` = primary_education_people,
  `Primary wealth analysis, people` = primary_wealth_people,
  `Focused low-education wealth, people` = focused_lowedu_wealth_people,
  `Onset / recovery / relapse events` = paste(incident_events, recovery_events, relapse_events, sep = " / "),
  `Verified deaths / post-difficulty deaths` = paste(total_death_events, disabled_death_events, sep = " / "),
  `Unknown vital or observation intervals` = unknown_vital_or_observation_intervals
)]

read_inf <- function(path, module_label) {
  z <- fread(path)
  z[, module := module_label]
  z
}
primary_inf <- rbindlist(list(
  read_inf(required_files[["education"]], "Education"),
  read_inf(required_files[["wealth"]], "Wealth")
), use.names = TRUE, fill = TRUE)
table2_numeric <- dcast(
  primary_inf[metric %chin% c("population_fdfle_gap", "population_recovery_relapse")],
  cohort + module ~ metric,
  value.var = c("point_estimate", "ci_low", "ci_high")
)
table2_numeric[, post_onset_supported :=
  ci_high_population_recovery_relapse < 0 | ci_low_population_recovery_relapse > 0]
table2_display <- table2_numeric[, .(
  Cohort = cohort,
  Exposure = module,
  `Low-minus-high FDFLE difference, years (95% CI)` = fmt_ci(
    point_estimate_population_fdfle_gap, ci_low_population_fdfle_gap, ci_high_population_fdfle_gap
  ),
  `Recovery-plus-relapse contribution, years (95% CI)` = fmt_ci(
    point_estimate_population_recovery_relapse,
    ci_low_population_recovery_relapse,
    ci_high_population_recovery_relapse
  ),
  `Post-onset CI excludes 0` = fifelse(post_onset_supported, "Yes", "No")
)]

focused <- fread(required_files[["focused"]])
level_metrics <- c(
  "population_high_tle", "population_low_tle",
  "population_high_fdfle", "population_low_fdfle",
  "population_high_difficulty_years", "population_low_difficulty_years",
  "population_high_fdfle_percent", "population_low_fdfle_percent"
)
focused_levels_numeric <- focused[metric %chin% level_metrics]
focused_levels_display <- focused_levels_numeric[, .(
  Cohort = cohort,
  Group = fifelse(grepl("_high_", metric), "High wealth", "Low wealth"),
  Outcome = fcase(
    grepl("_tle$", metric), "Total life years through age 100",
    grepl("_fdfle$", metric), "Functional-difficulty-free years through age 100",
    grepl("difficulty_years$", metric), "Years with functional difficulty through age 100",
    grepl("fdfle_percent$", metric), "FDFLE as percentage of truncated life years"
  ),
  `Estimate (95% CI)` = fmt_ci(point_estimate, ci_low, ci_high)
)]
focused_levels_display[, group_order := match(Group, c("High wealth", "Low wealth"))]
focused_levels_display[, outcome_order := match(Outcome, c(
  "Total life years through age 100",
  "Functional-difficulty-free years through age 100",
  "Years with functional difficulty through age 100",
  "FDFLE as percentage of truncated life years"
))]
focused_levels_display[, cohort_order := match(Cohort, c("CHARLS", "ELSA", "HRS", "MHAS"))]
setorder(focused_levels_display, cohort_order, group_order, outcome_order)
focused_levels_display[, c("cohort_order", "group_order", "outcome_order") := NULL]

component_names <- c(
  population_fdfle_gap = "Total low-minus-high FDFLE difference",
  population_contribution_initial_state = "Initial-state composition",
  population_contribution_onset = "Onset",
  population_contribution_recovery = "Recovery",
  population_contribution_relapse = "Relapse",
  population_contribution_post_difficulty_mortality = "Post-difficulty mortality",
  population_contribution_pre_difficulty_mortality = "Pre-difficulty mortality",
  population_recovery_relapse = "Recovery plus relapse"
)
focused_components_numeric <- focused[metric %chin% names(component_names)]
focused_components_display <- focused_components_numeric[, .(
  Cohort = cohort,
  Component = unname(component_names[metric]),
  `Contribution, years (95% CI)` = fmt_ci(point_estimate, ci_low, ci_high)
)]
focused_components_display[, component_order := match(Component, unname(component_names))]
focused_components_display[, cohort_order := match(Cohort, c("CHARLS", "ELSA", "HRS", "MHAS"))]
setorder(focused_components_display, cohort_order, component_order)
focused_components_display[, c("cohort_order", "component_order") := NULL]

threshold_inf <- rbindlist(list(
  read_inf(required_files[["focused"]], "Primary: any of nine difficulties"),
  read_inf(required_files[["adl"]], "ADL-only"),
  read_inf(required_files[["two"]], "At least two of nine difficulties")
), use.names = TRUE, fill = TRUE)
threshold_table <- threshold_inf[
  metric %chin% c("population_fdfle_gap", "population_recovery_relapse"),
  .(
    cohort, definition = module, metric, point_estimate, ci_low, ci_high,
    estimate_ci = fmt_ci(point_estimate, ci_low, ci_high),
    ci_excludes_zero = ci_high < 0 | ci_low > 0,
    valid_replicates
  )
]

copy_files <- list(
  TableS02_confirmed_state = file.path(outputs, "05_confirmed_state_sensitivity", "confirmed_state_cluster_robust_estimates.csv"),
  TableS03_continuous_time = file.path(outputs, "15_consolidated_evidence", "continuous_time_process_intensity_ratios.csv"),
  TableS04_enriched_covariates = file.path(outputs, "15_consolidated_evidence", "enriched_covariate_transition_ratios.csv"),
  TableS05_nondeath_ipcw = file.path(outputs, "15_consolidated_evidence", "nondeath_ipcw_transition_ratios.csv"),
  TableS06_point_sensitivities = file.path(outputs, "15_consolidated_evidence", "focused_point_sensitivity_matrix.csv"),
  TableS07_flow = file.path(outputs, "06_flow_missingness_support", "sample_flow_by_module_v1_1.csv"),
  TableS08_missingness = file.path(outputs, "06_flow_missingness_support", "functional_item_proxy_weight_missingness_v1_1.csv"),
  TableS09_wealth_timing = file.path(outputs, "01_revision_risksets", "wealth_measurement_timing_audit.csv"),
  TableS10_interval_distribution = file.path(outputs, "01_revision_risksets", "interval_distribution_by_ses.csv"),
  TableS11_focused_baseline = file.path(outputs, "11_focused_baseline_support", "focused_loweducation_wealth_baseline_characteristics.csv"),
  TableS12_focused_support = file.path(outputs, "11_focused_baseline_support", "focused_loweducation_wealth_transition_support.csv"),
  TableS13_internal_calibration = file.path(outputs, "10_internal_survival_calibration", "model_vs_weighted_km_tle_calibration.csv"),
  TableS14_item_mapping = file.path(outputs, "13_data_dictionary", "functional_item_mapping_by_cohort_wave.csv"),
  TableS15_state_dictionary = file.path(outputs, "13_data_dictionary", "history_state_dictionary.csv")
)
missing_copy <- copy_files[!file.exists(unlist(copy_files))]
if (length(missing_copy)) stop("Missing supplementary source files: ", paste(names(missing_copy), collapse = ", "))

fwrite(table1_numeric, file.path(source_dir, "Table1_numeric_source.csv"))
fwrite(table1_display, file.path(out_dir, "Table1_cohort_flow_and_support.csv"))
fwrite(table2_numeric, file.path(source_dir, "Table2_numeric_source.csv"))
fwrite(table2_display, file.path(out_dir, "Table2_primary_education_wealth_results.csv"))
fwrite(focused_levels_numeric, file.path(source_dir, "Table3A_numeric_source.csv"))
fwrite(focused_levels_display, file.path(out_dir, "Table3A_focused_absolute_levels.csv"))
fwrite(focused_components_numeric, file.path(source_dir, "Table3B_numeric_source.csv"))
fwrite(focused_components_display, file.path(out_dir, "Table3B_focused_six_block_decomposition.csv"))
fwrite(threshold_table, file.path(supp_dir, "TableS01_threshold_sensitivity_household_bootstrap.csv"))

for (nm in names(copy_files)) {
  z <- fread(copy_files[[nm]])
  fwrite(z, file.path(supp_dir, paste0(nm, ".csv")))
}

inference_modules <- list(
  primary_education = dirname(required_files[["education"]]),
  primary_wealth = dirname(required_files[["wealth"]]),
  focused_broad = dirname(required_files[["focused"]]),
  focused_adl_only = dirname(required_files[["adl"]]),
  focused_at_least_two = dirname(required_files[["two"]])
)
read_module_file <- function(module_name, filename) {
  z <- fread(file.path(inference_modules[[module_name]], filename))
  z[, module := module_name]
  z
}
heterogeneity_all <- rbindlist(lapply(
  names(inference_modules), read_module_file,
  filename = "cross_cohort_heterogeneity.csv"
), use.names = TRUE, fill = TRUE)
pairwise_all <- rbindlist(lapply(
  names(inference_modules), read_module_file,
  filename = "cross_cohort_pairwise_contrasts.csv"
), use.names = TRUE, fill = TRUE)
monte_carlo_all <- rbindlist(lapply(
  names(inference_modules), read_module_file,
  filename = "bootstrap_monte_carlo_stability.csv"
), use.names = TRUE, fill = TRUE)
replicate_qc_all <- rbindlist(lapply(
  names(inference_modules), read_module_file,
  filename = "bootstrap_replicate_qc.csv"
), use.names = TRUE, fill = TRUE)

primary_component_metrics <- c(
  "population_fdfle_gap",
  "population_contribution_initial_state",
  "population_contribution_onset",
  "population_contribution_recovery",
  "population_contribution_relapse",
  "population_contribution_post_difficulty_mortality",
  "population_contribution_pre_difficulty_mortality",
  "population_recovery_relapse"
)
primary_decomposition <- primary_inf[
  metric %chin% primary_component_metrics,
  .(
    cohort, exposure = module, metric,
    point_estimate, ci_low, ci_high, valid_replicates,
    estimate_ci = fmt_ci(point_estimate, ci_low, ci_high)
  )
]

psu_comparison_file <- file.path(outputs, "14_stratified_psu_inference", "household_vs_stratified_psu_ci_comparison.csv")
if (!file.exists(psu_comparison_file)) {
  stop("Missing stratified-PSU comparison file: ", psu_comparison_file)
}
psu_comparison <- fread(psu_comparison_file)

fwrite(heterogeneity_all, file.path(supp_dir, "TableS16_cross_cohort_heterogeneity.csv"))
fwrite(pairwise_all, file.path(supp_dir, "TableS17_all_pairwise_independent_cohort_contrasts.csv"))
fwrite(monte_carlo_all, file.path(supp_dir, "TableS18_bootstrap_monte_carlo_stability.csv"))
fwrite(replicate_qc_all, file.path(supp_dir, "TableS19_bootstrap_replicate_qc.csv"))
fwrite(primary_decomposition, file.path(supp_dir, "TableS20_primary_six_block_decomposition.csv"))
fwrite(psu_comparison, file.path(supp_dir, "TableS21_household_vs_stratified_psu_inference.csv"))

notes <- c(
  "# Revision v1.1 table titles and notes",
  "",
  "## Table 1. Cohort flow, baseline characteristics, and transition support",
  "",
  "Entry functional difficulty uses the primary nine-item broad definition. Women are expressed among participants with known sex; ELSA had 706 entrants with missing sex, explaining the denominator. Unknown vital or observation intervals were never coded as death.",
  "",
  "## Table 2. Primary education and wealth differences in population-initialised functional-difficulty-free life expectancy",
  "",
  "All differences and contributions are low minus high. Negative values denote fewer functional-difficulty-free years in the lower socioeconomic group. Estimates are truncated to ages 60-100. Confidence intervals use household-cluster percentile bootstrap with replicate-specific wealth cut-points.",
  "",
  "## Table 3A. Absolute truncated life years within low education by wealth",
  "",
  "The focused wealth analysis is secondary. Total life years, FDFLE, and years with functional difficulty are model-based remaining years accumulated from age 60 through age 100. They are not complete lifetime expectancies beyond age 100.",
  "",
  "## Table 3B. Six-block decomposition of the focused low-education wealth difference",
  "",
  "The six Shapley components sum to the low-minus-high FDFLE difference apart from rounding. Replacements are descriptive and are not causal mediation or intervention effects. Percentages are not shown when the total-gap confidence interval includes zero."
)
writeLines(notes, file.path(out_dir, "TABLE_TITLES_AND_NOTES.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_54_revision_tables.txt"))
cat("Revision v1.1 tables built in ", out_dir, "\n", sep = "")
