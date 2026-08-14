#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
  library(tidyselect)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
derived_dir <- file.path(revision_root, "02_derived")
out_dir <- file.path(revision_root, "03_outputs", "01_revision_risksets")
log_dir <- file.path(revision_root, "06_logs")
for (d in c(derived_dir, out_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("38_build_revision_panel_and_risksets_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

base_derived <- file.path(root, "02_derived")
panel <- as.data.table(readRDS(file.path(base_derived, "individual_wave_panel_gate0.rds")))
intervals <- as.data.table(readRDS(file.path(base_derived, "transition_intervals_gate0.rds")))
mortality <- as.data.table(readRDS(file.path(base_derived, "formal_mortality_riskset.rds")))

source_paths <- c(
  CHARLS = Sys.getenv("RECOVERY_DIVIDE_CHARLS_FILE", unset = ""),
  ELSA = Sys.getenv("RECOVERY_DIVIDE_ELSA_FILE", unset = ""),
  HRS = Sys.getenv("RECOVERY_DIVIDE_HRS_FILE", unset = ""),
  MHAS = Sys.getenv("RECOVERY_DIVIDE_MHAS_FILE", unset = "")
)
stopifnot(all(file.exists(source_paths)))

id_chr <- function(x) {
  if (is.character(x)) return(trimws(x))
  out <- rep(NA_character_, length(x))
  ok <- is.finite(as.numeric(x))
  out[ok] <- sprintf("%.0f", as.numeric(x)[ok])
  out
}

read_selected <- function(path, variables) {
  meta <- read_dta(path, n_max = 0)
  available <- intersect(variables, names(meta))
  missing <- setdiff(variables, names(meta))
  if (length(missing)) cat("Missing requested variables in ", path, ": ", paste(missing, collapse = ", "), "\n", sep = "")
  as.data.table(read_dta(path, col_select = all_of(available)))
}

# Build a wave-level household/design map. Household is the binding resampling
# unit. PSU/stratum fields are retained only where the harmonised source exposes
# them; their absence is never replaced by a fabricated survey design.
map_parts <- list()

charls_vars <- c(
  "ID", "hhidc", "householdID_w1", "householdID",
  "hh1atotb", "hh2atotb", "h3atotb", "h4atotb"
)
charls <- read_selected(source_paths[["CHARLS"]], charls_vars)
charls[, person_id := id_chr(ID)]
charls[, household_fixed := fifelse(
  !is.na(hhidc) & nzchar(trimws(hhidc)), trimws(hhidc),
  fifelse(!is.na(householdID_w1) & nzchar(trimws(householdID_w1)), trimws(householdID_w1), substr(person_id, 1L, pmax(nchar(person_id) - 3L, 1L)))
)]
charls_wealth_vars <- c("hh1atotb", "hh2atotb", "h3atotb", "h4atotb")
map_parts[["CHARLS"]] <- rbindlist(lapply(seq_along(charls_wealth_vars), function(w) {
  data.table(
    cohort = "CHARLS", person_id = charls$person_id, wave = as.integer(w),
    wealth_corrected = as.numeric(charls[[charls_wealth_vars[[w]]]]),
    household_id = paste0("CHARLS::", charls$household_fixed),
    design_stratum = NA_character_, design_psu = NA_character_
  )
}))

elsa_hh <- paste0("hh", 1:10, "hhidc")
elsa_strat <- paste0("r", 1:10, "strat")
elsa_clust <- paste0("r", 1:10, "clust")
elsa <- read_selected(source_paths[["ELSA"]], c("idauniq", elsa_hh, elsa_strat, elsa_clust))
elsa[, person_id := id_chr(idauniq)]
map_parts[["ELSA"]] <- rbindlist(lapply(1:10, function(w) {
  hh <- id_chr(elsa[[elsa_hh[[w]]]])
  stratum <- id_chr(elsa[[elsa_strat[[w]]]])
  psu <- id_chr(elsa[[elsa_clust[[w]]]])
  data.table(
    cohort = "ELSA", person_id = elsa$person_id, wave = as.integer(w),
    wealth_corrected = NA_real_,
    household_id = fifelse(!is.na(hh) & nzchar(hh), paste0("ELSA::", hh), NA_character_),
    design_stratum = fifelse(!is.na(stratum) & nzchar(stratum), paste0("ELSA::", stratum), NA_character_),
    design_psu = fifelse(!is.na(psu) & nzchar(psu), paste0("ELSA::", stratum, "::", psu), NA_character_)
  )
}))

hrs <- read_selected(source_paths[["HRS"]], c("hhidpn", "hhid", "raestrat", "raehsamp"))
hrs[, person_id := id_chr(hhidpn)]
hrs_hh <- id_chr(hrs$hhid)
hrs_stratum <- id_chr(hrs$raestrat)
hrs_half <- id_chr(hrs$raehsamp)
hrs_map_person <- data.table(
  cohort = "HRS", person_id = hrs$person_id,
  household_id = fifelse(!is.na(hrs_hh) & nzchar(hrs_hh), paste0("HRS::", hrs_hh), NA_character_),
  design_stratum = fifelse(!is.na(hrs_stratum) & nzchar(hrs_stratum), paste0("HRS::", hrs_stratum), NA_character_),
  design_psu = fifelse(
    !is.na(hrs_stratum) & nzchar(hrs_stratum) & !is.na(hrs_half) & nzchar(hrs_half),
    paste0("HRS::", hrs_stratum, "::", hrs_half), NA_character_
  )
)

mhas <- read_selected(source_paths[["MHAS"]], c("rahhidnp", "unhhid"))
mhas[, person_id := id_chr(rahhidnp)]
mhas_hh <- id_chr(mhas$unhhid)
mhas_map_person <- data.table(
  cohort = "MHAS", person_id = mhas$person_id,
  household_id = fifelse(!is.na(mhas_hh) & nzchar(mhas_hh), paste0("MHAS::", mhas_hh), NA_character_),
  design_stratum = NA_character_, design_psu = NA_character_
)

wave_map <- rbindlist(map_parts, use.names = TRUE, fill = TRUE)
stopifnot(!anyDuplicated(wave_map, by = c("cohort", "person_id", "wave")))

panel[, c("household_id", "design_stratum", "design_psu") := NA_character_]
panel[wave_map, on = .(cohort, person_id, wave), `:=`(
  household_id = i.household_id,
  design_stratum = i.design_stratum,
  design_psu = i.design_psu
)]
panel[hrs_map_person, on = .(cohort, person_id), `:=`(
  household_id = i.household_id,
  design_stratum = i.design_stratum,
  design_psu = i.design_psu
)]
panel[mhas_map_person, on = .(cohort, person_id), `:=`(
  household_id = i.household_id,
  design_stratum = i.design_stratum,
  design_psu = i.design_psu
)]

# Correct the documented CHARLS wave-1/2 mapping defect while retaining the
# already correct wave-3/4 values from the same harmonised release.
panel[wave_map[cohort == "CHARLS"], on = .(cohort, person_id, wave), wealth := i.wealth_corrected]

# Carry household/design identifiers to adapters not present in the harmonised
# core (ELSA wave 11 and CHARLS wave 5) from the nearest prior observed record.
setorder(panel, cohort, person_id, wave)
locf_character <- function(x) {
  last_value <- NA_character_
  for (i in seq_along(x)) {
    if (!is.na(x[[i]]) && nzchar(x[[i]])) {
      last_value <- x[[i]]
    } else if (!is.na(last_value)) {
      x[[i]] <- last_value
    }
  }
  x
}
panel[, household_id := locf_character(household_id), by = .(cohort, person_id)]
panel[, design_stratum := locf_character(design_stratum), by = .(cohort, person_id)]
panel[, design_psu := locf_character(design_psu), by = .(cohort, person_id)]
panel[is.na(household_id) | !nzchar(household_id), household_id := paste0(cohort, "::PERSON::", person_id)]
panel[, household_fallback_to_person := grepl("::PERSON::", household_id, fixed = TRUE)]

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
  wealth_entry_weight = respondent_weight,
  wealth_entry_household_id = household_id,
  wealth_entry_design_stratum = design_stratum,
  wealth_entry_design_psu = design_psu
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

wealth_entry <- wealth_cutpoints[wealth_entry, on = .(cohort, country, wealth_entry_wave), nomatch = 0]
wealth_entry[, wealth3 := fifelse(
  wealth_entry_value <= q33, "low",
  fifelse(wealth_entry_value <= q67, "middle", "high")
)]

person_ses <- merge(
  education_person,
  wealth_entry[, .(
    cohort, country, person_id, wealth3, wealth_entry_wave, wealth_entry_time,
    wealth_entry_age, wealth_entry_value, wealth_entry_weight,
    wealth_entry_household_id, wealth_entry_design_stratum, wealth_entry_design_psu,
    cut_method, q33, q67
  )],
  by = c("cohort", "country", "person_id"), all.x = TRUE, sort = FALSE
)

# Entry and history-depth audit.
entry_function <- panel[
  observed_state == TRUE & !is.na(age) & age >= 60 & age <= 95
][order(cohort, person_id, interview_time, wave), .SD[1L], by = .(cohort, person_id)]
first_difficulty <- panel[
  observed_state == TRUE & difficulty == 1 & !is.na(age) & age >= 60 & age <= 95,
  .(first_observed_difficulty_time = min(interview_time)), by = .(cohort, person_id)
]
entry_audit_person <- merge(
  entry_function[, .(
    cohort, country, person_id, entry_wave = wave, entry_time = interview_time,
    entry_age = age, entry_difficulty = difficulty, entry_proxy = proxy,
    entry_weight = respondent_weight, entry_female = female
  )],
  person_ses,
  by = c("cohort", "country", "person_id"), all.x = TRUE, sort = FALSE
)
entry_audit_person <- merge(
  entry_audit_person, first_difficulty,
  by = c("cohort", "person_id"), all.x = TRUE, sort = FALSE
)
entry_audit_person[, `:=`(
  wealth_delay_years = wealth_entry_time - entry_time,
  wealth_after_entry = !is.na(wealth_entry_time) & wealth_entry_time > entry_time + 1e-8,
  wealth_after_first_observed_difficulty = !is.na(wealth_entry_time) &
    !is.na(first_observed_difficulty_time) & wealth_entry_time > first_observed_difficulty_time + 1e-8,
  entry_age_60_64 = entry_age >= 60 & entry_age < 65
)]

# Rebuild the formal risk sets from the locked functional histories, replacing
# only the corrected SES/household fields. No state transition is redefined.
drop_ses <- c(
  "education3_fixed", "wealth3", "wealth_entry_wave", "wealth_entry_time",
  "wealth_entry_age", "wealth_entry_value", "wealth_entry_weight",
  "household_id", "design_stratum", "design_psu"
)
intervals[, (intersect(drop_ses, names(intervals))) := NULL]
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

mortality[, (intersect(drop_ses, names(mortality))) := NULL]
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

# Initial-state data use the first functional observation at or after the
# exposure becomes observable. This avoids assigning a future wealth category
# to an earlier functional state.
panel_ses <- merge(
  panel,
  person_ses[, .(
    cohort, person_id, education3_fixed, wealth3, wealth_entry_time,
    wealth_entry_wave, wealth_entry_household_id
  )],
  by = c("cohort", "person_id"), all.x = TRUE, sort = FALSE
)
entry_education <- panel_ses[
  observed_state == TRUE & age >= 60 & age <= 95 & !is.na(education3_fixed)
][order(cohort, person_id, interview_time, wave), .SD[1L], by = .(cohort, person_id)]
entry_wealth_state <- panel_ses[
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
  entry_wealth_state[, .(
    cohort, country, person_id, exposure = "wealth", level = wealth3,
    entry_wave = wave, entry_time = interview_time, entry_age = age,
    entry_difficulty = difficulty, female, proxy, entry_weight = respondent_weight,
    household_id = wealth_entry_household_id,
    design_stratum, design_psu
  )]
), use.names = TRUE, fill = TRUE)

# Audits saved as aggregate files only.
wealth_timing_summary <- entry_audit_person[, .(
  entrants = .N,
  wealth_observed_n = sum(!is.na(wealth3)),
  wealth_missing_n = sum(is.na(wealth3)),
  wealth_missing_percent = 100 * mean(is.na(wealth3)),
  wealth_after_entry_n = sum(wealth_after_entry, na.rm = TRUE),
  wealth_after_entry_percent_all = 100 * mean(wealth_after_entry, na.rm = TRUE),
  wealth_after_first_difficulty_n = sum(wealth_after_first_observed_difficulty, na.rm = TRUE),
  wealth_after_first_difficulty_percent_all = 100 * mean(wealth_after_first_observed_difficulty, na.rm = TRUE),
  median_delay_years_observed = median(wealth_delay_years[!is.na(wealth_delay_years)], na.rm = TRUE),
  p25_delay_years_observed = quantile(wealth_delay_years[!is.na(wealth_delay_years)], 0.25, names = FALSE, na.rm = TRUE),
  p75_delay_years_observed = quantile(wealth_delay_years[!is.na(wealth_delay_years)], 0.75, names = FALSE, na.rm = TRUE)
), by = .(cohort, country)]

household_summary <- person_ses[, .(
  people = uniqueN(person_id),
  wealth_observed_people = uniqueN(person_id[!is.na(wealth3)]),
  household_clusters = uniqueN(wealth_entry_household_id[!is.na(wealth3)]),
  singleton_fallback_people = sum(grepl("::PERSON::", wealth_entry_household_id, fixed = TRUE), na.rm = TRUE),
  design_stratum_observed_people = sum(!is.na(wealth_entry_design_stratum)),
  design_psu_observed_people = sum(!is.na(wealth_entry_design_psu))
), by = .(cohort, country)]

entry_state_summary <- initial_state_data[level %in% c("low", "high"), .(
  people = uniqueN(person_id),
  difficulty_n = sum(entry_difficulty == 1, na.rm = TRUE),
  difficulty_percent_unweighted = 100 * mean(entry_difficulty == 1, na.rm = TRUE),
  entry_age_median = median(entry_age, na.rm = TRUE),
  entry_age_60_64_n = sum(entry_age >= 60 & entry_age < 65, na.rm = TRUE)
), by = .(cohort, exposure, level)]

interval_timing_summary <- function(z, time_col, exposure_name) {
  eligible_col <- if (exposure_name == "wealth") "wealth_time_eligible" else "education_eligible"
  level_col <- if (exposure_name == "wealth") "wealth3" else "education3_fixed"
  z[get(eligible_col) == TRUE & get(level_col) %in% c("low", "high"), .(
    intervals = .N,
    people = uniqueN(person_id),
    mean_interval_years = mean(interval_years),
    median_interval_years = median(interval_years),
    p25_interval_years = quantile(interval_years, 0.25, names = FALSE),
    p75_interval_years = quantile(interval_years, 0.75, names = FALSE)
  ), by = .(cohort, level = get(level_col))][, exposure := exposure_name]
}
interval_summary <- rbindlist(list(
  interval_timing_summary(function_risk, "interview_time", "education"),
  interval_timing_summary(function_risk, "interview_time", "wealth")
), use.names = TRUE, fill = TRUE)
setcolorder(interval_summary, c("cohort", "exposure", "level", setdiff(names(interval_summary), c("cohort", "exposure", "level"))))

saveRDS(panel, file.path(derived_dir, "individual_wave_panel_revision_v1_1.rds"), compress = "xz")
saveRDS(person_ses, file.path(derived_dir, "person_fixed_ses_revision_v1_1.rds"), compress = "xz")
saveRDS(function_risk, file.path(derived_dir, "formal_function_transition_riskset_revision_v1_1.rds"), compress = "xz")
saveRDS(mortality, file.path(derived_dir, "formal_mortality_riskset_revision_v1_1.rds"), compress = "xz")
saveRDS(initial_state_data, file.path(derived_dir, "initial_state_data_revision_v1_1.rds"), compress = "xz")

fwrite(wealth_cutpoints, file.path(out_dir, "wealth_entry_wave_cutpoints_revision_v1_1.csv"))
fwrite(wealth_timing_summary, file.path(out_dir, "wealth_measurement_timing_audit.csv"))
fwrite(household_summary, file.path(out_dir, "household_and_design_support.csv"))
fwrite(entry_state_summary, file.path(out_dir, "entry_state_support.csv"))
fwrite(interval_summary, file.path(out_dir, "interval_distribution_by_ses.csv"))
fwrite(entry_audit_person[, .(
  cohort, person_id, entry_wave, entry_age, entry_difficulty, wealth_entry_wave,
  wealth_entry_age, wealth_delay_years, wealth_after_entry,
  wealth_after_first_observed_difficulty, entry_age_60_64
)], file.path(derived_dir, "restricted_person_timing_audit_internal.csv"))

cat("Revision risk-set build complete\n")
cat("Panel rows: ", nrow(panel), "\n", sep = "")
cat("Function intervals: ", nrow(function_risk), "\n", sep = "")
cat("Mortality intervals: ", nrow(mortality), "\n", sep = "")
print(wealth_timing_summary)
print(household_summary)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_38_revision_risksets.txt"))
