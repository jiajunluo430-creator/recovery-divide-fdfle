#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 41_build_threshold_revision_risksets.R adl_only|at_least_two")
definition <- tolower(args[[1L]])
stopifnot(definition %in% c("adl_only", "at_least_two"))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
base_derived <- file.path(revision_root, "02_derived")
derived_dir <- file.path(base_derived, "sensitivity", definition)
out_dir <- file.path(revision_root, "03_outputs", paste0("04_threshold_support_", definition))
log_dir <- file.path(revision_root, "06_logs")
for (d in c(derived_dir, out_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

panel <- as.data.table(readRDS(file.path(base_derived, "individual_wave_panel_revision_v1_1.rds")))
person_ses <- as.data.table(readRDS(file.path(base_derived, "person_fixed_ses_revision_v1_1.rds")))
setorder(panel, cohort, person_id, wave)

adl_items <- c("dressa", "batha", "eata", "beda", "toilta")
if (definition == "adl_only") {
  item_matrix <- as.matrix(panel[, ..adl_items])
  complete_alt <- rowSums(!is.na(item_matrix)) == length(adl_items)
  count_alt <- rowSums(item_matrix, na.rm = FALSE)
  difficulty_alt <- fifelse(complete_alt, as.integer(count_alt > 0), NA_integer_)
} else {
  complete_alt <- panel$function_complete9
  count_alt <- panel$difficulty_count
  difficulty_alt <- fifelse(complete_alt, as.integer(count_alt >= 2), NA_integer_)
}

panel[, `:=`(
  difficulty_primary = difficulty,
  difficulty_count_alt = count_alt,
  function_complete_alt = complete_alt,
  difficulty = difficulty_alt
)]
panel[, observed_state := !is.na(difficulty) & in_wave == 1]
panel[, ever_difficulty_prior := shift(cummax(fifelse(observed_state & difficulty == 1, 1L, 0L)), fill = 0L), by = .(cohort, person_id)]
panel[, recovery_observation := as.integer(observed_state & difficulty == 0 & ever_difficulty_prior == 1L)]
panel[, ever_recovery_prior := shift(cummax(recovery_observation), fill = 0L), by = .(cohort, person_id)]
panel[, history_state := fifelse(!observed_state, NA_character_,
  fifelse(difficulty == 0,
    fifelse(ever_difficulty_prior == 1L, "R1", "I0"),
    fifelse(ever_recovery_prior == 1L, "D2", "D1")
  ))]

obs <- panel[observed_state == TRUE]
setorder(obs, cohort, person_id, interview_time, wave)
obs[, `:=`(
  dest_wave_living = shift(wave, type = "lead"),
  dest_time_living = shift(interview_time, type = "lead"),
  dest_state_living = shift(history_state, type = "lead")
), by = .(cohort, person_id)]
study_end <- panel[, .(study_end_time = max(interview_time, na.rm = TRUE)), by = cohort]
intervals <- study_end[obs, on = "cohort"]
intervals[, death_before_next := !is.na(death_time) & death_time > interview_time &
  death_time <= study_end_time & (is.na(dest_time_living) | death_time < dest_time_living)]
intervals[, destination := fifelse(death_before_next, "DEAD", dest_state_living)]
intervals[, destination_time := fifelse(death_before_next, death_time, dest_time_living)]
intervals[, destination_wave := fifelse(death_before_next, NA_integer_, dest_wave_living)]
intervals[, `:=`(
  interval_years = destination_time - interview_time,
  origin_state = history_state,
  origin_wave = wave,
  origin_age = age,
  origin_year = nominal_year
)]
intervals[, long_mhas_gap := !is.na(destination_wave) & cohort == "MHAS" & origin_wave == 2L & destination_wave == 3L]
intervals[, primary_interval := !is.na(destination) & !is.na(interval_years) &
  interval_years >= 1 & interval_years <= 4 & !long_mhas_gap &
  !is.na(origin_age) & origin_age >= 60 & origin_age <= 95 & !living_after_death]
allowed_pairs <- c(
  "I0->I0", "I0->D1", "I0->DEAD",
  "D1->D1", "D1->R1", "D1->DEAD",
  "R1->R1", "R1->D2", "R1->DEAD",
  "D2->D2", "D2->R1", "D2->DEAD"
)
intervals[, transition := paste0(origin_state, "->", destination)]
intervals[, impossible_transition := !is.na(destination) & !transition %in% allowed_pairs]

function_risk <- intervals[primary_interval == TRUE & impossible_transition == FALSE]
function_risk <- merge(
  function_risk,
  person_ses[, .(
    cohort, person_id, education3_fixed, wealth3, wealth_entry_wave,
    wealth_entry_time, wealth_entry_age, wealth_entry_value, wealth_entry_weight,
    household_id = wealth_entry_household_id,
    design_stratum = wealth_entry_design_stratum,
    design_psu = wealth_entry_design_psu
  )],
  by = c("cohort", "person_id"), all.x = TRUE, sort = FALSE
)
function_risk[, `:=`(
  wealth_time_eligible = !is.na(wealth3) & !is.na(wealth_entry_time) & interview_time >= wealth_entry_time - 1e-8,
  education_eligible = !is.na(education3_fixed),
  event_onset = as.integer(origin_state == "I0" & destination == "D1"),
  event_recovery_d1 = as.integer(origin_state == "D1" & destination == "R1"),
  event_recovery_d2 = as.integer(origin_state == "D2" & destination == "R1"),
  event_recovery = as.integer(origin_state %in% c("D1", "D2") & destination == "R1"),
  event_relapse = as.integer(origin_state == "R1" & destination == "D2")
)]

# Rebuild mortality risk sets with the alternative history states.
origin <- panel[observed_state == TRUE, .(
  cohort, country, person_id, origin_wave = wave, next_wave = wave + 1L,
  origin_time = interview_time, origin_age = age, origin_state = history_state,
  origin_year = nominal_year, female, proxy, origin_weight = respondent_weight,
  death_time, death_source, living_after_death
)]
next_wave <- panel[, .(
  cohort, person_id, next_wave = wave, next_time = interview_time,
  next_observed_state = observed_state, next_in_wave = in_wave, next_iwstat = iwstat
)]
mortality <- merge(origin, next_wave, by = c("cohort", "person_id", "next_wave"), all = FALSE, sort = FALSE)
mortality[, death_event := !is.na(death_time) & death_time > origin_time & death_time <= next_time]
mortality[, known_alive_boundary := next_observed_state == TRUE |
  (!is.na(next_in_wave) & next_in_wave == 1) | (!is.na(next_iwstat) & next_iwstat %in% c(1, 4))]
mortality[, mortality_outcome_known := death_event | known_alive_boundary]
mortality[, endpoint_time := fifelse(death_event, death_time, next_time)]
mortality[, interval_years := endpoint_time - origin_time]
mortality[, long_mhas_gap := cohort == "MHAS" & origin_wave == 2L & next_wave == 3L]
mortality[, mortality_window_supported :=
  (cohort == "CHARLS" & origin_wave %in% 1:4) |
  (cohort == "HRS" & origin_wave %in% 2:15) |
  (cohort == "ELSA" & origin_wave %in% 1:5) |
  (cohort == "MHAS" & origin_wave %in% 1:5)]
mortality[, primary_mortality_interval := mortality_window_supported & mortality_outcome_known &
  !is.na(interval_years) & interval_years > 0 & interval_years <= 4 & !long_mhas_gap &
  !is.na(origin_age) & origin_age >= 60 & origin_age <= 95 & !living_after_death]
mortality <- merge(
  mortality,
  person_ses[, .(
    cohort, person_id, education3_fixed, wealth3, wealth_entry_wave,
    wealth_entry_time, wealth_entry_age, wealth_entry_value, wealth_entry_weight,
    household_id = wealth_entry_household_id,
    design_stratum = wealth_entry_design_stratum,
    design_psu = wealth_entry_design_psu
  )],
  by = c("cohort", "person_id"), all.x = TRUE, sort = FALSE
)
mortality[, `:=`(
  wealth_time_eligible = !is.na(wealth3) & !is.na(wealth_entry_time) & origin_time >= wealth_entry_time - 1e-8,
  education_eligible = !is.na(education3_fixed),
  event_death = as.integer(death_event)
)]

panel_ses <- merge(
  panel,
  person_ses[, .(cohort, person_id, education3_fixed, wealth3, wealth_entry_time, wealth_entry_household_id)],
  by = c("cohort", "person_id"), all.x = TRUE, sort = FALSE
)
entry_education <- panel_ses[
  observed_state == TRUE & age >= 60 & age <= 95 & !is.na(education3_fixed)
][order(cohort, person_id, interview_time, wave), .SD[1L], by = .(cohort, person_id)]
entry_wealth <- panel_ses[
  observed_state == TRUE & age >= 60 & age <= 95 & !is.na(wealth3) &
    !is.na(wealth_entry_time) & interview_time >= wealth_entry_time - 1e-8
][order(cohort, person_id, interview_time, wave), .SD[1L], by = .(cohort, person_id)]
initial_state_data <- rbindlist(list(
  entry_education[, .(
    cohort, country, person_id, exposure = "education", level = education3_fixed,
    entry_wave = wave, entry_time = interview_time, entry_age = age,
    entry_difficulty = difficulty, female, proxy, entry_weight = respondent_weight,
    household_id, design_stratum, design_psu
  )],
  entry_wealth[, .(
    cohort, country, person_id, exposure = "wealth", level = wealth3,
    entry_wave = wave, entry_time = interview_time, entry_age = age,
    entry_difficulty = difficulty, female, proxy, entry_weight = respondent_weight,
    household_id = wealth_entry_household_id, design_stratum, design_psu
  )]
), use.names = TRUE, fill = TRUE)

support <- function_risk[wealth_time_eligible == TRUE & education3_fixed == "low", .(
  intervals = .N, people = uniqueN(person_id),
  onset_events = sum(event_onset), recovery_events = sum(event_recovery),
  relapse_events = sum(event_relapse)
), by = cohort]

saveRDS(function_risk, file.path(derived_dir, "formal_function_transition_riskset_revision_v1_1.rds"), compress = "xz")
saveRDS(mortality, file.path(derived_dir, "formal_mortality_riskset_revision_v1_1.rds"), compress = "xz")
saveRDS(initial_state_data, file.path(derived_dir, "initial_state_data_revision_v1_1.rds"), compress = "xz")
fwrite(support, file.path(out_dir, "low_education_wealth_event_support.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, paste0("sessionInfo_41_", definition, ".txt")))
cat("Built threshold risk sets: ", definition, "\n", sep = "")
print(support)

