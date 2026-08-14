#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
derived_root <- file.path(root, "02_derived")
threshold_derived_root <- file.path(derived_root, "threshold_sensitivity")
out_root <- file.path(root, "03_outputs", "09_threshold_sensitivity")
log_dir <- file.path(root, "06_logs")
for (d in c(threshold_derived_root, out_root, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("08b_build_threshold_sensitivity_risksets_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

panel_base <- as.data.table(readRDS(file.path(derived_root, "individual_wave_panel_gate0.rds")))
person_ses <- as.data.table(readRDS(file.path(derived_root, "person_fixed_ses_internal.rds")))
setorder(panel_base, cohort, person_id, interview_time, wave)

all_items <- c("dressa", "batha", "eata", "beda", "toilta", "shopa", "mealsa", "medsa", "moneya")
adl_items <- all_items[1:5]
allowed_pairs <- c(
  "I0->I0", "I0->D1", "I0->DEAD",
  "D1->D1", "D1->R1", "D1->DEAD",
  "R1->R1", "R1->D2", "R1->DEAD",
  "D2->D2", "D2->R1", "D2->DEAD"
)
thresholds <- c("adl_only", "at_least_two", "permissive_partial")

define_difficulty <- function(panel, threshold_name) {
  all_matrix <- as.matrix(panel[, ..all_items])
  if (threshold_name == "adl_only") {
    adl_matrix <- as.matrix(panel[, ..adl_items])
    complete <- rowSums(!is.na(adl_matrix)) == length(adl_items)
    count <- rowSums(adl_matrix, na.rm = FALSE)
    return(fifelse(complete, as.integer(count > 0), NA_integer_))
  }
  if (threshold_name == "at_least_two") {
    complete <- rowSums(!is.na(all_matrix)) == length(all_items)
    count <- rowSums(all_matrix, na.rm = FALSE)
    return(fifelse(complete, as.integer(count >= 2), NA_integer_))
  }
  observed_n <- rowSums(!is.na(all_matrix))
  any_difficulty <- rowSums(all_matrix == 1, na.rm = TRUE) > 0
  all_observed_clear <- observed_n == length(all_items) & rowSums(all_matrix == 1, na.rm = TRUE) == 0
  fifelse(any_difficulty, 1L, fifelse(all_observed_clear, 0L, NA_integer_))
}

build_threshold <- function(threshold_name) {
  cat("\nBuilding threshold: ", threshold_name, "\n", sep = "")
  panel <- copy(panel_base)
  panel[, difficulty_threshold := define_difficulty(panel, threshold_name)]
  panel[, observed_threshold := in_wave == 1 & !is.na(difficulty_threshold)]
  panel[, ever_difficulty_prior_threshold := shift(
    cummax(fifelse(observed_threshold & difficulty_threshold == 1, 1L, 0L)), fill = 0L
  ), by = .(cohort, person_id)]
  panel[, recovery_observation_threshold := as.integer(
    observed_threshold & difficulty_threshold == 0 & ever_difficulty_prior_threshold == 1L
  )]
  panel[, ever_recovery_prior_threshold := shift(
    cummax(recovery_observation_threshold), fill = 0L
  ), by = .(cohort, person_id)]
  panel[, history_state_threshold := fifelse(!observed_threshold, NA_character_,
    fifelse(difficulty_threshold == 0,
      fifelse(ever_difficulty_prior_threshold == 1L, "R1", "I0"),
      fifelse(ever_recovery_prior_threshold == 1L, "D2", "D1")
    ))]

  obs <- panel[observed_threshold == TRUE]
  setorder(obs, cohort, person_id, interview_time, wave)
  obs[, `:=`(
    dest_wave_living = shift(wave, type = "lead"),
    dest_time_living = shift(interview_time, type = "lead"),
    dest_age_living = shift(age, type = "lead"),
    dest_state_living = shift(history_state_threshold, type = "lead")
  ), by = .(cohort, person_id)]
  study_end <- panel[, .(study_end_time = max(interview_time, na.rm = TRUE)), by = cohort]
  intervals <- study_end[obs, on = "cohort"]
  intervals[, death_before_next := !is.na(death_time) & death_time > interview_time &
    death_time <= study_end_time & (is.na(dest_time_living) | death_time < dest_time_living)]
  intervals[, destination := fifelse(death_before_next, "DEAD", dest_state_living)]
  intervals[, destination_time := fifelse(death_before_next, death_time, dest_time_living)]
  intervals[, destination_wave := fifelse(death_before_next, NA_integer_, dest_wave_living)]
  intervals[, interval_years := destination_time - interview_time]
  intervals[, `:=`(
    origin_state = history_state_threshold,
    origin_wave = wave,
    origin_age = age,
    origin_year = nominal_year
  )]
  intervals[, long_mhas_gap := !is.na(destination_wave) & cohort == "MHAS" & origin_wave == 2L & destination_wave == 3L]
  intervals[, primary_interval := !is.na(destination) & !is.na(interval_years) &
    interval_years >= 1 & interval_years <= 4 & !long_mhas_gap &
    !is.na(origin_age) & origin_age >= 60 & origin_age <= 95 & !living_after_death]
  intervals[, transition := paste0(origin_state, "->", destination)]
  intervals[, impossible_transition := !is.na(destination) & !transition %in% allowed_pairs]

  function_risk <- intervals[primary_interval == TRUE & impossible_transition == FALSE]
  function_risk <- merge(
    function_risk,
    person_ses[, .(cohort, person_id, education3_fixed, wealth3, wealth_entry_wave)],
    by = c("cohort", "person_id"), all.x = TRUE, sort = FALSE
  )
  function_risk[, `:=`(
    event_onset = as.integer(origin_state == "I0" & destination == "D1"),
    event_recovery_d1 = as.integer(origin_state == "D1" & destination == "R1"),
    event_recovery_d2 = as.integer(origin_state == "D2" & destination == "R1"),
    event_relapse = as.integer(origin_state == "R1" & destination == "D2")
  )]

  origin <- panel[observed_threshold == TRUE, .(
    cohort, country, person_id,
    origin_wave = wave, next_wave = wave + 1L,
    origin_time = interview_time, origin_age = age,
    origin_state = history_state_threshold, origin_year = nominal_year,
    female, proxy, origin_weight = respondent_weight,
    death_time, death_source, living_after_death
  )]
  next_wave <- panel[, .(
    cohort, person_id, next_wave = wave,
    next_time = interview_time,
    next_observed_state = observed_threshold,
    next_in_wave = in_wave, next_iwstat = iwstat
  )]
  mortality_risk <- merge(
    origin, next_wave, by = c("cohort", "person_id", "next_wave"),
    all.x = FALSE, all.y = FALSE, sort = FALSE
  )
  mortality_risk[, death_event := !is.na(death_time) & death_time > origin_time & death_time <= next_time]
  mortality_risk[, known_alive_boundary := next_observed_state == TRUE |
    (!is.na(next_in_wave) & next_in_wave == 1) |
    (!is.na(next_iwstat) & next_iwstat %in% c(1, 4))]
  mortality_risk[, mortality_outcome_known := death_event | known_alive_boundary]
  mortality_risk[, endpoint_time := fifelse(death_event, death_time, next_time)]
  mortality_risk[, interval_years := endpoint_time - origin_time]
  mortality_risk[, long_mhas_gap := cohort == "MHAS" & origin_wave == 2L & next_wave == 3L]
  mortality_risk[, mortality_window_supported :=
    (cohort == "CHARLS" & origin_wave %in% 1:4) |
    (cohort == "HRS" & origin_wave %in% 2:15) |
    (cohort == "ELSA" & origin_wave %in% 1:5) |
    (cohort == "MHAS" & origin_wave %in% 1:5)]
  mortality_risk[, primary_mortality_interval :=
    mortality_window_supported & mortality_outcome_known &
    !is.na(interval_years) & interval_years > 0 & interval_years <= 4 &
    !long_mhas_gap & !is.na(origin_age) & origin_age >= 60 & origin_age <= 95 &
    !living_after_death]
  mortality_risk <- merge(
    mortality_risk,
    person_ses[, .(cohort, person_id, education3_fixed, wealth3, wealth_entry_wave)],
    by = c("cohort", "person_id"), all.x = TRUE, sort = FALSE
  )
  mortality_risk[, event_death := as.integer(death_event)]

  gate <- function_risk[, .(
    incident_events = sum(origin_state == "I0" & destination == "D1"),
    incident_persons = uniqueN(person_id[origin_state == "I0" & destination == "D1"]),
    recovery_events = sum(origin_state %in% c("D1", "D2") & destination == "R1"),
    recovery_persons = uniqueN(person_id[origin_state %in% c("D1", "D2") & destination == "R1"]),
    relapse_events = sum(origin_state == "R1" & destination == "D2"),
    relapse_persons = uniqueN(person_id[origin_state == "R1" & destination == "D2"]),
    disabled_deaths = sum(origin_state %in% c("D1", "D2") & destination == "DEAD"),
    total_deaths = sum(destination == "DEAD"),
    impossible_transitions = sum(impossible_transition, na.rm = TRUE)
  ), by = .(cohort, country)]
  gate[, threshold := threshold_name]

  mortality_gate <- mortality_risk[primary_mortality_interval == TRUE, .(
    supported_intervals = .N, supported_persons = uniqueN(person_id),
    deaths = sum(event_death),
    disabled_deaths = sum(event_death & origin_state %in% c("D1", "D2")),
    pre_disability_deaths = sum(event_death & origin_state == "I0"),
    recovered_state_deaths = sum(event_death & origin_state == "R1")
  ), by = .(cohort, country)]
  mortality_gate[, threshold := threshold_name]

  missingness <- panel[, .(
    wave_rows = .N,
    observed_threshold_states = sum(observed_threshold),
    missing_threshold_states = sum(in_wave == 1 & !observed_threshold),
    difficulty_states = sum(observed_threshold & difficulty_threshold == 1),
    independent_states = sum(observed_threshold & difficulty_threshold == 0)
  ), by = .(cohort, wave)]
  missingness[, threshold := threshold_name]

  threshold_dir <- file.path(threshold_derived_root, threshold_name)
  dir.create(threshold_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(function_risk, file.path(threshold_dir, "formal_function_transition_riskset.rds"), compress = "xz")
  saveRDS(mortality_risk, file.path(threshold_dir, "formal_mortality_riskset.rds"), compress = "xz")
  list(gate = gate, mortality_gate = mortality_gate, missingness = missingness)
}

results <- lapply(thresholds, build_threshold)
gate <- rbindlist(lapply(results, `[[`, "gate"), fill = TRUE)
mortality_gate <- rbindlist(lapply(results, `[[`, "mortality_gate"), fill = TRUE)
missingness <- rbindlist(lapply(results, `[[`, "missingness"), fill = TRUE)
setcolorder(gate, c("threshold", setdiff(names(gate), "threshold")))
setcolorder(mortality_gate, c("threshold", setdiff(names(mortality_gate), "threshold")))
setcolorder(missingness, c("threshold", setdiff(names(missingness), "threshold")))
fwrite(gate, file.path(out_root, "threshold_gate0_event_counts.csv"))
fwrite(mortality_gate, file.path(out_root, "threshold_mortality_gate.csv"))
fwrite(missingness, file.path(out_root, "threshold_wave_state_missingness.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_threshold_risksets.txt"))

cat("\nThreshold Gate-0 events:\n")
print(gate)
cat("\nThreshold mortality gates:\n")
print(mortality_gate)
cat("Completed: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
