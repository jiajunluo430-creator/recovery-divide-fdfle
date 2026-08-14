#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
derived_dir <- file.path(root, "02_derived")
out_dir <- file.path(root, "03_outputs", "10_exploratory_upgrade_support")
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("13_exploratory_upgrade_support_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

panel <- as.data.table(readRDS(file.path(derived_dir, "individual_wave_panel_gate0.rds")))
function_risk <- as.data.table(readRDS(file.path(derived_dir, "formal_function_transition_riskset.rds")))
mortality_risk <- as.data.table(readRDS(file.path(derived_dir, "formal_mortality_riskset.rds")))
person_ses <- as.data.table(readRDS(file.path(derived_dir, "person_fixed_ses_internal.rds")))

cohort_order <- c("CHARLS", "HRS", "ELSA", "MHAS")
exposure_order <- c("education", "wealth")
process_order <- c("onset", "recovery", "relapse", "death_pre", "death_post")

function_risk[, event_recovery := as.integer(origin_state %in% c("D1", "D2") & destination == "R1")]
function_risk[, sex := fifelse(female == 1, "female", fifelse(female == 0, "male", NA_character_))]
mortality_risk[, sex := fifelse(female == 1, "female", fifelse(female == 0, "male", NA_character_))]

person_ses[, joint_ses9 := fifelse(
  is.na(education3_fixed) | is.na(wealth3), NA_character_,
  paste0("edu_", education3_fixed, "__wealth_", wealth3)
)]
person_ses[, joint_ses4 := fcase(
  education3_fixed == "high" & wealth3 == "high", "high_high",
  education3_fixed == "low" & wealth3 == "high", "lowedu_highwealth",
  education3_fixed == "high" & wealth3 == "low", "highedu_lowwealth",
  education3_fixed == "low" & wealth3 == "low", "low_low",
  default = NA_character_
)]

function_risk <- merge(
  function_risk,
  person_ses[, .(cohort, person_id, joint_ses9, joint_ses4)],
  by = c("cohort", "person_id"), all.x = TRUE, sort = FALSE
)
mortality_risk <- merge(
  mortality_risk,
  person_ses[, .(cohort, person_id, joint_ses9, joint_ses4)],
  by = c("cohort", "person_id"), all.x = TRUE, sort = FALSE
)

# Recovery phase is tied to immediately adjacent scheduled observed waves.
setorder(panel, cohort, person_id, wave)
panel[, `:=`(
  previous_wave = shift(wave),
  previous_observed = shift(observed_state),
  previous_history_state = shift(history_state)
), by = .(cohort, person_id)]
panel[, recovery_phase := fcase(
  observed_state == TRUE & history_state == "R1" &
    previous_wave == wave - 1L & previous_observed == TRUE &
    previous_history_state %in% c("D1", "D2"), "early_recovery",
  observed_state == TRUE & history_state == "R1" &
    previous_wave == wave - 1L & previous_observed == TRUE &
    previous_history_state == "R1", "sustained_recovery",
  observed_state == TRUE & history_state == "R1", "unclassified",
  default = NA_character_
)]
phase_map <- unique(panel[
  observed_state == TRUE & history_state == "R1",
  .(cohort, person_id, origin_wave = wave, recovery_phase)
])
function_risk <- merge(
  function_risk, phase_map,
  by = c("cohort", "person_id", "origin_wave"), all.x = TRUE, sort = FALSE
)
mortality_risk <- merge(
  mortality_risk, phase_map,
  by = c("cohort", "person_id", "origin_wave"), all.x = TRUE, sort = FALSE
)

process_specs <- list(
  onset = list(source = "function", origins = "I0", event = "event_onset"),
  recovery = list(source = "function", origins = c("D1", "D2"), event = "event_recovery"),
  relapse = list(source = "function", origins = "R1", event = "event_relapse"),
  death_pre = list(source = "mortality", origins = "I0", event = "event_death"),
  death_post = list(source = "mortality", origins = c("D1", "R1", "D2"), event = "event_death")
)

get_process_data <- function(process_name) {
  spec <- process_specs[[process_name]]
  if (spec$source == "function") {
    z <- copy(function_risk[origin_state %in% spec$origins])
  } else {
    z <- copy(mortality_risk[
      primary_mortality_interval == TRUE & origin_state %in% spec$origins
    ])
  }
  z[, event := as.integer(get(spec$event))]
  z
}

sex_support <- rbindlist(lapply(process_order, function(process_name) {
  z <- get_process_data(process_name)
  rbindlist(lapply(exposure_order, function(exposure_name) {
    level_col <- if (exposure_name == "education") "education3_fixed" else "wealth3"
    z[!is.na(sex) & get(level_col) %in% c("low", "high"), .(
      intervals = .N,
      persons = uniqueN(person_id),
      events = sum(event),
      person_years = sum(interval_years)
    ), by = .(cohort, sex, level = get(level_col))][, `:=`(
      exposure = exposure_name,
      process = process_name
    )]
  }))
}), fill = TRUE)
setcolorder(sex_support, c("cohort", "exposure", "sex", "level", "process", "intervals", "persons", "events", "person_years"))

joint_person_support <- person_ses[!is.na(joint_ses4), .(
  persons = uniqueN(person_id)
), by = .(cohort, joint_ses4)]

joint_process_support <- rbindlist(lapply(process_order, function(process_name) {
  z <- get_process_data(process_name)
  z[!is.na(joint_ses4), .(
    intervals = .N,
    persons = uniqueN(person_id),
    events = sum(event),
    person_years = sum(interval_years)
  ), by = .(cohort, joint_ses4)][, process := process_name]
}), fill = TRUE)
setcolorder(joint_process_support, c("cohort", "joint_ses4", "process", "intervals", "persons", "events", "person_years"))

phase_function_support <- rbindlist(lapply(exposure_order, function(exposure_name) {
  level_col <- if (exposure_name == "education") "education3_fixed" else "wealth3"
  function_risk[
    origin_state == "R1" & recovery_phase %in% c("early_recovery", "sustained_recovery") &
      get(level_col) %in% c("low", "high"),
    .(
      intervals = .N,
      persons = uniqueN(person_id),
      relapses = sum(event_relapse),
      person_years = sum(interval_years)
    ),
    by = .(cohort, recovery_phase, level = get(level_col))
  ][, exposure := exposure_name]
}), fill = TRUE)

phase_mortality_support <- rbindlist(lapply(exposure_order, function(exposure_name) {
  level_col <- if (exposure_name == "education") "education3_fixed" else "wealth3"
  mortality_risk[
    primary_mortality_interval == TRUE & origin_state == "R1" &
      recovery_phase %in% c("early_recovery", "sustained_recovery") &
      get(level_col) %in% c("low", "high"),
    .(
      intervals = .N,
      persons = uniqueN(person_id),
      deaths = sum(event_death),
      person_years = sum(interval_years)
    ),
    by = .(cohort, recovery_phase, level = get(level_col))
  ][, exposure := exposure_name]
}), fill = TRUE)

age_bands <- c("60-69", "70-79", "80-89", "90-95")
age_support <- rbindlist(lapply(c("recovery", "relapse"), function(process_name) {
  z <- get_process_data(process_name)
  z[, age_band := cut(
    origin_age,
    breaks = c(60, 70, 80, 90, 96),
    right = FALSE,
    labels = age_bands
  )]
  rbindlist(lapply(exposure_order, function(exposure_name) {
    level_col <- if (exposure_name == "education") "education3_fixed" else "wealth3"
    z[!is.na(age_band) & get(level_col) %in% c("low", "high"), .(
      intervals = .N,
      persons = uniqueN(person_id),
      events = sum(event),
      person_years = sum(interval_years)
    ), by = .(cohort, age_band, level = get(level_col))][, `:=`(
      exposure = exposure_name,
      process = process_name
    )]
  }))
}), fill = TRUE)

support_gate <- rbindlist(list(
  sex_support[process %in% c("recovery", "relapse"), .(
    decisive_cells = .N,
    min_recovery_events = min(events[process == "recovery"]),
    min_relapse_events = min(events[process == "relapse"]),
    supported = min(events[process == "recovery"]) >= 100 & min(events[process == "relapse"]) >= 50
  ), by = cohort][, module := "sex"],
  joint_process_support[process %in% c("recovery", "relapse"), .(
    decisive_cells = .N,
    min_recovery_events = min(events[process == "recovery"]),
    min_relapse_events = min(events[process == "relapse"]),
    supported = min(events[process == "recovery"]) >= 100 & min(events[process == "relapse"]) >= 50
  ), by = cohort][, module := "joint_ses"],
  phase_function_support[, .(
    decisive_cells = .N,
    min_recovery_events = NA_integer_,
    min_relapse_events = min(relapses),
    supported = min(relapses) >= 50
  ), by = cohort][, module := "recovery_phase"]
), fill = TRUE)

module_decision <- support_gate[, .(
  supported_cohorts = sum(supported),
  total_cohorts = .N,
  decision = ifelse(sum(supported) >= 3, "SUPPORTED_FOR_POINT_MODEL", "DESCRIPTIVE_OR_STOP")
), by = module]

fwrite(sex_support, file.path(out_dir, "sex_process_support.csv"))
fwrite(joint_person_support, file.path(out_dir, "joint_ses_person_support.csv"))
fwrite(joint_process_support, file.path(out_dir, "joint_ses_process_support.csv"))
fwrite(phase_function_support, file.path(out_dir, "recovery_phase_relapse_support.csv"))
fwrite(phase_mortality_support, file.path(out_dir, "recovery_phase_mortality_support.csv"))
fwrite(age_support, file.path(out_dir, "age_band_recovery_relapse_support.csv"))
fwrite(support_gate, file.path(out_dir, "module_support_gate_by_cohort.csv"))
fwrite(module_decision, file.path(out_dir, "module_support_decision.csv"))

cat("Exploratory module support decisions:\n")
print(module_decision)
cat("\nSupport by cohort:\n")
print(support_gate[order(module, cohort)])
cat("\nJoint SES person support:\n")
print(joint_person_support[order(cohort, joint_ses4)])

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_exploratory_upgrade_support.txt"))
cat("Completed: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
