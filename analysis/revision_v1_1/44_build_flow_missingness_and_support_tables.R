options(stringsAsFactors = FALSE, width = 220, warn = 1)

if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
derived_dir <- file.path(revision_root, "02_derived")
point_dir <- file.path(revision_root, "03_outputs", "02_revision_point")
out_dir <- file.path(revision_root, "03_outputs", "06_flow_missingness_support")
log_dir <- file.path(revision_root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

panel <- as.data.table(readRDS(file.path(derived_dir, "individual_wave_panel_revision_v1_1.rds")))
function_risk <- as.data.table(readRDS(file.path(
  derived_dir, "formal_function_transition_riskset_revision_v1_1.rds"
)))
mortality_risk <- as.data.table(readRDS(file.path(
  derived_dir, "formal_mortality_riskset_revision_v1_1.rds"
)))
person_ses <- as.data.table(readRDS(file.path(derived_dir, "person_fixed_ses_revision_v1_1.rds")))
cohorts <- c("CHARLS", "ELSA", "HRS", "MHAS")

baseline <- panel[
  observed_state == TRUE & is.finite(age) & age >= 60 & age <= 95
][order(cohort, person_id, interview_time), .SD[1L], by = .(cohort, person_id)]
baseline <- person_ses[baseline, on = .(cohort, person_id)]

baseline_table <- baseline[, .(
  country = first(country),
  eligible_entrants = .N,
  age_mean = mean(age),
  age_sd = sd(age),
  women_n = sum(female == 1, na.rm = TRUE),
  sex_known_n = sum(!is.na(female)),
  women_percent_among_known_sex = 100 * mean(female == 1, na.rm = TRUE),
  sex_missing_n = sum(is.na(female)),
  baseline_difficulty_n = sum(difficulty == 1, na.rm = TRUE),
  baseline_difficulty_percent = 100 * mean(difficulty == 1, na.rm = TRUE),
  proxy_status_available_n = sum(!is.na(proxy)),
  proxy_interview_n = sum(proxy == 1, na.rm = TRUE),
  proxy_percent_among_available = if (sum(!is.na(proxy))) 100 * mean(proxy == 1, na.rm = TRUE) else NA_real_,
  positive_entry_weight_n = sum(is.finite(respondent_weight) & respondent_weight > 0),
  positive_entry_weight_percent = 100 * mean(is.finite(respondent_weight) & respondent_weight > 0),
  education_observed_n = sum(!is.na(education3_fixed)),
  education_missing_n = sum(is.na(education3_fixed)),
  wealth_observed_n = sum(!is.na(wealth_entry_value)),
  wealth_missing_n = sum(is.na(wealth_entry_value)),
  household_id_available_n = sum(!is.na(wealth_entry_household_id) & nzchar(wealth_entry_household_id))
), by = cohort]
baseline_table[, women_display := sprintf(
  "%d/%d (%.1f%%)", women_n, sex_known_n, women_percent_among_known_sex
)]

flow_parts <- list()
add_flow <- function(cohort_name, module_name, stages) {
  previous <- NA_integer_
  for (i in seq_along(stages)) {
    ids <- unique(stages[[i]]$ids)
    n <- length(ids)
    flow_parts[[length(flow_parts) + 1L]] <<- data.table(
      cohort = cohort_name,
      module = module_name,
      stage_order = i,
      stage = stages[[i]]$label,
      people = n,
      excluded_from_previous_stage = if (is.na(previous)) NA_integer_ else previous - n
    )
    previous <- n
  }
}

for (cohort_name in cohorts) {
  base_ids <- baseline[cohort == cohort_name, person_id]
  fr <- function_risk[cohort == cohort_name]
  formal_ids <- intersect(base_ids, unique(fr$person_id))

  add_flow(cohort_name, "primary_education", list(
    list(label = "Observed functional state at age 60-95", ids = base_ids),
    list(label = "At least one formal adjacent functional interval", ids = formal_ids),
    list(label = "Education observed and eligible", ids = intersect(
      formal_ids, unique(fr[education_eligible == TRUE, person_id])
    ))
  ))

  wealth_observed_ids <- intersect(formal_ids, person_ses[
    cohort == cohort_name & !is.na(wealth_entry_value), person_id
  ])
  wealth_time_ids <- intersect(formal_ids, unique(fr[wealth_time_eligible == TRUE, person_id]))
  wealth_extreme_ids <- intersect(
    wealth_time_ids,
    person_ses[cohort == cohort_name & wealth3 %chin% c("high", "low"), person_id]
  )
  add_flow(cohort_name, "primary_wealth", list(
    list(label = "Observed functional state at age 60-95", ids = base_ids),
    list(label = "At least one formal adjacent functional interval", ids = formal_ids),
    list(label = "Wealth observed", ids = wealth_observed_ids),
    list(label = "At least one interval beginning on/after first wealth measure", ids = wealth_time_ids),
    list(label = "High or low within-wave wealth tertile", ids = wealth_extreme_ids)
  ))

  low_education_ids <- intersect(formal_ids, person_ses[
    cohort == cohort_name & education3_fixed == "low", person_id
  ])
  focused_wealth_observed <- intersect(low_education_ids, wealth_observed_ids)
  focused_time <- intersect(low_education_ids, wealth_time_ids)
  focused_extreme <- intersect(low_education_ids, wealth_extreme_ids)
  focused_household <- intersect(
    focused_extreme,
    person_ses[
      cohort == cohort_name & !is.na(wealth_entry_household_id) & nzchar(wealth_entry_household_id),
      person_id
    ]
  )
  add_flow(cohort_name, "focused_low_education_wealth", list(
    list(label = "At least one formal adjacent functional interval", ids = formal_ids),
    list(label = "Low education", ids = low_education_ids),
    list(label = "Wealth observed", ids = focused_wealth_observed),
    list(label = "At least one interval beginning on/after first wealth measure", ids = focused_time),
    list(label = "High or low within-wave wealth tertile", ids = focused_extreme),
    list(label = "Household cluster identifier available", ids = focused_household)
  ))
}
sample_flow <- rbindlist(flow_parts)

item_vars <- c("dressa", "batha", "eata", "beda", "toilta", "shopa", "mealsa", "medsa", "moneya")
age_eligible_inwave <- panel[in_wave == 1 & is.finite(age) & age >= 60 & age <= 95]
missing_parts <- list()
for (field in c(item_vars, "difficulty", "proxy", "respondent_weight")) {
  missing_parts[[length(missing_parts) + 1L]] <- age_eligible_inwave[, .(
    field = field,
    denominator_interviews = .N,
    missing_or_unavailable_n = sum(is.na(get(field))),
    missing_or_unavailable_percent = 100 * mean(is.na(get(field)))
  ), by = cohort]
}
interview_missingness <- rbindlist(missing_parts)
interview_missingness[, interpretation := fifelse(
  field == "proxy",
  "Unavailable may reflect cohort/wave non-collection; not imputed as self interview",
  fifelse(field == "respondent_weight", "Nonpositive weights audited separately", "Item nonresponse")
)]

ses_missingness <- baseline[, .(
  entrants = .N,
  education_missing_n = sum(is.na(education3_fixed)),
  education_missing_percent = 100 * mean(is.na(education3_fixed)),
  wealth_missing_n = sum(is.na(wealth_entry_value)),
  wealth_missing_percent = 100 * mean(is.na(wealth_entry_value)),
  household_missing_among_wealth_observed_n = sum(
    !is.na(wealth_entry_value) & (is.na(wealth_entry_household_id) | !nzchar(wealth_entry_household_id))
  ),
  design_stratum_available_among_wealth_observed_n = sum(
    !is.na(wealth_entry_value) & !is.na(wealth_entry_design_stratum)
  ),
  design_psu_available_among_wealth_observed_n = sum(
    !is.na(wealth_entry_value) & !is.na(wealth_entry_design_psu)
  )
), by = cohort]

vital_status_intervals <- mortality_risk[, .(
  intervals = .N,
  people = uniqueN(person_id),
  deaths = sum(death_event == TRUE, na.rm = TRUE),
  person_years = sum(interval_years, na.rm = TRUE)
), by = .(
  cohort,
  mortality_exclusion,
  mortality_outcome_known
)]
vital_status_summary <- mortality_risk[, .(
  intervals_total = .N,
  people_total = uniqueN(person_id),
  verified_death_intervals = sum(mortality_outcome_known == TRUE & death_event == TRUE),
  known_alive_or_observed_intervals = sum(mortality_outcome_known == TRUE & death_event == FALSE),
  unknown_vital_or_observation_intervals = sum(mortality_outcome_known == FALSE),
  primary_supported_intervals = sum(primary_mortality_interval == TRUE),
  primary_supported_deaths = sum(primary_mortality_interval == TRUE & death_event == TRUE)
), by = cohort]

transition_support <- fread(file.path(point_dir, "transition_model_qc.csv"))[, .(
  cohort,
  module,
  process,
  intervals,
  people,
  events,
  person_years,
  converged,
  warnings
)]

fwrite(baseline_table, file.path(out_dir, "baseline_characteristics_corrected_v1_1.csv"))
fwrite(sample_flow, file.path(out_dir, "sample_flow_by_module_v1_1.csv"))
fwrite(interview_missingness, file.path(out_dir, "functional_item_proxy_weight_missingness_v1_1.csv"))
fwrite(ses_missingness, file.path(out_dir, "ses_and_design_missingness_v1_1.csv"))
fwrite(vital_status_intervals, file.path(out_dir, "vital_status_interval_audit_v1_1.csv"))
fwrite(vital_status_summary, file.path(out_dir, "death_vs_non_death_loss_summary_v1_1.csv"))
fwrite(transition_support, file.path(out_dir, "transition_event_support_all_primary_modules_v1_1.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_44_flow_missingness.txt"))

cat("Corrected baseline table:\n")
print(baseline_table)
cat("\nFocused sample flow:\n")
print(sample_flow[module == "focused_low_education_wealth"])
cat("\nDeath/non-death loss summary:\n")
print(vital_status_summary)
