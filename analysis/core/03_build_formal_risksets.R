#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
derived_dir <- file.path(root, "02_derived")
out_dir <- file.path(root, "03_outputs", "03_formal_risksets")
log_dir <- file.path(root, "06_logs")
for (d in c(derived_dir, out_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("03_build_formal_risksets_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

panel_path <- file.path(derived_dir, "individual_wave_panel_gate0.rds")
interval_path <- file.path(derived_dir, "transition_intervals_gate0.rds")
stopifnot(file.exists(panel_path), file.exists(interval_path))
panel <- as.data.table(readRDS(panel_path))
intervals <- as.data.table(readRDS(interval_path))

last_nonmissing <- function(z) {
  idx <- which(!is.na(z))
  if (length(idx)) z[idx[[length(idx)]]] else z[NA_integer_]
}

weighted_cut <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) / sum(w) >= p)[1L]]
}

# Fixed education and first valid wealth at/after age 60.
setorder(panel, cohort, person_id, interview_time, wave)
education_person <- panel[, .(
  education3_fixed = last_nonmissing(education3)
), by = .(cohort, country, person_id)]

wealth_entry <- panel[
  in_wave == 1 & !is.na(age) & age >= 60 & is.finite(wealth)
][order(cohort, person_id, interview_time, wave), .SD[1L], by = .(cohort, country, person_id)]
wealth_entry <- wealth_entry[, .(
  cohort, country, person_id,
  wealth_entry_wave = wave,
  wealth_entry_time = interview_time,
  wealth_entry_age = age,
  wealth_entry_value = wealth,
  wealth_entry_weight = respondent_weight
)]

wealth_cutpoints <- wealth_entry[, {
  weighted_ok <- any(is.finite(wealth_entry_weight) & wealth_entry_weight > 0)
  ww <- if (weighted_ok) wealth_entry_weight else rep(1, .N)
  list(
    n = .N,
    positive_weight_n = sum(is.finite(wealth_entry_weight) & wealth_entry_weight > 0),
    cut_method = if (weighted_ok) "weighted_empirical" else "unweighted_fallback",
    q33 = weighted_cut(wealth_entry_value, ww, 1 / 3),
    q67 = weighted_cut(wealth_entry_value, ww, 2 / 3)
  )
}, by = .(cohort, country, wealth_entry_wave)]

wealth_entry <- wealth_cutpoints[wealth_entry,
  on = .(cohort, country, wealth_entry_wave), nomatch = 0]
wealth_entry[, wealth3 := fifelse(
  !is.na(wealth_entry_value) & wealth_entry_value <= q33, "low",
  fifelse(!is.na(wealth_entry_value) & wealth_entry_value <= q67, "middle",
    fifelse(!is.na(wealth_entry_value), "high", NA_character_))
)]

person_ses <- merge(
  education_person,
  wealth_entry[, .(
    cohort, country, person_id, wealth3, wealth_entry_wave,
    wealth_entry_time, wealth_entry_age, wealth_entry_value,
    wealth_entry_weight, cut_method, q33, q67
  )],
  by = c("cohort", "country", "person_id"), all.x = TRUE, sort = FALSE
)

# Observed living-state risk intervals for onset, recovery, and relapse.
function_risk <- intervals[
  primary_interval == TRUE & impossible_transition == FALSE
]
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

# Adjacent scheduled-wave mortality risk set. Unknown status is not coded alive.
origin <- panel[observed_state == TRUE, .(
  cohort, country, person_id,
  origin_wave = wave, next_wave = wave + 1L,
  origin_time = interview_time, origin_age = age,
  origin_state = history_state, origin_year = nominal_year,
  female, proxy, origin_weight = respondent_weight,
  death_time, death_source, living_after_death
)]
next_wave <- panel[, .(
  cohort, person_id, next_wave = wave,
  next_time = interview_time,
  next_observed_state = observed_state,
  next_in_wave = in_wave,
  next_iwstat = iwstat
)]
mortality_risk <- merge(
  origin, next_wave,
  by = c("cohort", "person_id", "next_wave"), all.x = FALSE, all.y = FALSE,
  sort = FALSE
)
mortality_risk[, death_event := !is.na(death_time) & death_time > origin_time & death_time <= next_time]
mortality_risk[, known_alive_boundary :=
  next_observed_state == TRUE |
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
mortality_risk[, mortality_exclusion := fifelse(!mortality_window_supported, "unsupported_cohort_wave",
  fifelse(!mortality_outcome_known, "unknown_vital_or_observation_status",
    fifelse(is.na(interval_years), "missing_interval",
      fifelse(interval_years <= 0, "nonpositive_interval",
        fifelse(interval_years > 4, "interval_gt_4y",
          fifelse(long_mhas_gap, "prespecified_mhas_2003_2012_gap",
            fifelse(is.na(origin_age), "missing_origin_age",
              fifelse(origin_age < 60, "origin_age_lt_60",
                fifelse(origin_age > 95, "origin_age_gt_95",
                  fifelse(living_after_death, "living_after_recorded_death", "included"))))))))))]
mortality_risk <- merge(
  mortality_risk,
  person_ses[, .(cohort, person_id, education3_fixed, wealth3, wealth_entry_wave)],
  by = c("cohort", "person_id"), all.x = TRUE, sort = FALSE
)
mortality_risk[, event_death := as.integer(death_event)]

# Aggregate data feasibility and exposure support before any model fitting.
ses_summary <- function(level_col, exposure_name) {
  z <- person_ses[, .(persons = uniqueN(person_id)),
    by = .(cohort, country, level = get(level_col))]
  z[, exposure := exposure_name]
  setcolorder(z, c("cohort", "country", "exposure", "level", "persons"))
  z
}

function_summary <- function(level_col, exposure_name, origin_name, event_col, process_name) {
  z <- function_risk[origin_state == origin_name, .(
    intervals = .N,
    persons = uniqueN(person_id),
    events = sum(get(event_col)),
    person_years = sum(interval_years)
  ), by = .(cohort, level = get(level_col))]
  z[, `:=`(exposure = exposure_name, process = process_name)]
  setcolorder(z, c("cohort", "exposure", "level", "process", "intervals", "persons", "events", "person_years"))
  z
}

mortality_summary <- function(level_col, exposure_name) {
  z <- mortality_risk[primary_mortality_interval == TRUE, .(
    intervals = .N,
    persons = uniqueN(person_id),
    deaths = sum(event_death),
    person_years = sum(interval_years)
  ), by = .(cohort, origin_state, level = get(level_col))]
  z[, exposure := exposure_name]
  setcolorder(z, c("cohort", "origin_state", "exposure", "level", "intervals", "persons", "deaths", "person_years"))
  z
}

ses_distribution <- rbindlist(list(
  ses_summary("education3_fixed", "education"),
  ses_summary("wealth3", "wealth")
), use.names = TRUE, fill = TRUE)

function_events_by_ses <- rbindlist(lapply(c("education", "wealth"), function(exposure_name) {
  level_col <- if (exposure_name == "education") "education3_fixed" else "wealth3"
  rbindlist(list(
    function_summary(level_col, exposure_name, "I0", "event_onset", "onset"),
    function_summary(level_col, exposure_name, "D1", "event_recovery_d1", "recovery_d1"),
    function_summary(level_col, exposure_name, "D2", "event_recovery_d2", "recovery_d2"),
    function_summary(level_col, exposure_name, "R1", "event_relapse", "relapse")
  ))
}), use.names = TRUE, fill = TRUE)

mortality_events_by_ses <- rbindlist(list(
  mortality_summary("education3_fixed", "education"),
  mortality_summary("wealth3", "wealth")
), use.names = TRUE, fill = TRUE)

mortality_support_audit <- mortality_risk[, .(
  intervals = .N,
  persons = uniqueN(person_id),
  deaths = sum(death_event, na.rm = TRUE),
  disabled_deaths = sum(death_event & origin_state %in% c("D1", "D2"), na.rm = TRUE)
), by = .(cohort, mortality_exclusion)]

formal_gate <- mortality_risk[primary_mortality_interval == TRUE, .(
  supported_intervals = .N,
  supported_persons = uniqueN(person_id),
  deaths = sum(event_death),
  disabled_deaths = sum(event_death & origin_state %in% c("D1", "D2")),
  pre_disability_deaths = sum(event_death & origin_state == "I0"),
  recovered_state_deaths = sum(event_death & origin_state == "R1")
), by = .(cohort, country)]

saveRDS(person_ses, file.path(derived_dir, "person_fixed_ses_internal.rds"), compress = "xz")
saveRDS(function_risk, file.path(derived_dir, "formal_function_transition_riskset.rds"), compress = "xz")
saveRDS(mortality_risk, file.path(derived_dir, "formal_mortality_riskset.rds"), compress = "xz")
fwrite(wealth_cutpoints, file.path(out_dir, "wealth_entry_wave_cutpoints.csv"))
fwrite(ses_distribution, file.path(out_dir, "person_ses_distribution.csv"))
fwrite(function_events_by_ses, file.path(out_dir, "function_events_by_ses.csv"))
fwrite(mortality_events_by_ses, file.path(out_dir, "supported_mortality_events_by_ses.csv"))
fwrite(mortality_support_audit, file.path(out_dir, "mortality_support_and_exclusion_audit.csv"))
fwrite(formal_gate, file.path(out_dir, "formal_mortality_gate.csv"))

cat("Formal mortality gate:\n")
print(formal_gate)
cat("\nMortality support audit:\n")
print(mortality_support_audit[order(cohort, mortality_exclusion)])
cat("\nSES distributions:\n")
print(ses_distribution[!is.na(level)][order(cohort, exposure, level)])

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_formal_risksets.txt"))
cat("Completed: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
