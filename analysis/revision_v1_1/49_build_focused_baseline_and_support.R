#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
derived_dir <- file.path(revision_root, "02_derived")
out_dir <- file.path(revision_root, "03_outputs", "11_focused_baseline_support")
log_dir <- file.path(revision_root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
sink(file.path(log_dir, paste0("49_focused_baseline_support_", stamp, ".log")), split = TRUE)
on.exit(sink(), add = TRUE)

person_ses <- as.data.table(readRDS(file.path(derived_dir, "person_fixed_ses_revision_v1_1.rds")))
initial <- as.data.table(readRDS(file.path(derived_dir, "initial_state_data_revision_v1_1.rds")))
function_risk <- as.data.table(readRDS(file.path(derived_dir, "formal_function_transition_riskset_revision_v1_1.rds")))
mortality_risk <- as.data.table(readRDS(file.path(derived_dir, "formal_mortality_riskset_revision_v1_1.rds")))

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_quantile <- function(x, w, p = 0.5) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  o <- order(x)
  x <- x[o]
  w <- w[o]
  x[which(cumsum(w) / sum(w) >= p)[1L]]
}

wealth_entry <- copy(initial[exposure == "wealth"])
wealth_entry <- person_ses[wealth_entry, on = .(cohort, person_id), nomatch = 0]
focused_entry <- wealth_entry[
  education3_fixed == "low" & wealth3 %chin% c("high", "low") &
    is.finite(entry_weight) & entry_weight > 0
]

focused_baseline <- focused_entry[, .(
  people = uniqueN(person_id),
  households = uniqueN(wealth_entry_household_id),
  weighted_mean_entry_age = weighted_mean_safe(entry_age, entry_weight),
  weighted_female_percent = 100 * weighted_mean_safe(female, entry_weight),
  sex_missing_n = sum(is.na(female)),
  weighted_entry_difficulty_percent = 100 * weighted_mean_safe(entry_difficulty, entry_weight),
  proxy_observed_n = sum(!is.na(proxy)),
  weighted_proxy_percent_among_observed = if (sum(!is.na(proxy)) >= 30L) {
    100 * weighted_mean_safe(proxy, entry_weight)
  } else {
    NA_real_
  },
  proxy_missing_n = sum(is.na(proxy)),
  weighted_median_wealth_local_currency = weighted_quantile(wealth_entry_value, entry_weight, 0.5)
), by = .(cohort, country, wealth3)]
focused_baseline[, wealth_order := match(wealth3, c("high", "low"))]
setorder(focused_baseline, cohort, wealth_order)
focused_baseline[, wealth_order := NULL]

function_specs <- list(
  onset = list(origins = "I0", event = "event_onset"),
  recovery = list(origins = c("D1", "D2"), event = "event_recovery"),
  relapse = list(origins = "R1", event = "event_relapse")
)
mortality_specs <- list(
  death_pre = list(origins = "I0", event = "event_death"),
  death_post = list(origins = c("D1", "R1", "D2"), event = "event_death")
)

support_parts <- list()
for (cohort_name in unique(focused_entry$cohort)) {
  for (wealth_name in c("high", "low")) {
    for (process_name in names(function_specs)) {
      spec <- function_specs[[process_name]]
      z <- function_risk[
        cohort == cohort_name & wealth_time_eligible == TRUE &
          education3_fixed == "low" & wealth3 == wealth_name &
          origin_state %chin% spec$origins
      ]
      support_parts[[length(support_parts) + 1L]] <- data.table(
        cohort = cohort_name, wealth3 = wealth_name, process = process_name,
        intervals = nrow(z), people = uniqueN(z$person_id), households = uniqueN(z$household_id),
        events = sum(z[[spec$event]], na.rm = TRUE), person_years = sum(z$interval_years, na.rm = TRUE),
        proxy_observed_intervals = sum(!is.na(z$proxy)),
        proxy_intervals = sum(z$proxy == 1, na.rm = TRUE)
      )
    }
    for (process_name in names(mortality_specs)) {
      spec <- mortality_specs[[process_name]]
      z <- mortality_risk[
        cohort == cohort_name & primary_mortality_interval == TRUE &
          wealth_time_eligible == TRUE & education3_fixed == "low" &
          wealth3 == wealth_name & origin_state %chin% spec$origins
      ]
      support_parts[[length(support_parts) + 1L]] <- data.table(
        cohort = cohort_name, wealth3 = wealth_name, process = process_name,
        intervals = nrow(z), people = uniqueN(z$person_id), households = uniqueN(z$household_id),
        events = sum(z[[spec$event]], na.rm = TRUE), person_years = sum(z$interval_years, na.rm = TRUE),
        proxy_observed_intervals = sum(!is.na(z$proxy)),
        proxy_intervals = sum(z$proxy == 1, na.rm = TRUE)
      )
    }
  }
}
focused_support <- rbindlist(support_parts, use.names = TRUE, fill = TRUE)
focused_support[, proxy_percent_among_observed := fifelse(
  proxy_observed_intervals > 0,
  100 * proxy_intervals / proxy_observed_intervals,
  NA_real_
)]

ses_base <- initial[
  exposure %chin% c("education", "wealth") & !is.na(level) &
    is.finite(entry_weight) & entry_weight > 0
]
ses_totals <- ses_base[, .(weighted_total = sum(entry_weight)), by = .(cohort, exposure)]
ses_distribution <- ses_base[
  ,
  .(
    people = uniqueN(person_id),
    weighted_people = sum(entry_weight)
  ),
  by = .(cohort, country, exposure, level)
]
ses_distribution <- ses_totals[ses_distribution, on = .(cohort, exposure)]
ses_distribution[, weighted_percent := 100 * weighted_people / weighted_total]

focused_flow <- rbindlist(lapply(unique(focused_entry$cohort), function(cohort_name) {
  ids_entry <- focused_entry[cohort == cohort_name, unique(person_id)]
  ids_function <- function_risk[
    cohort == cohort_name & wealth_time_eligible == TRUE & education3_fixed == "low" &
      wealth3 %chin% c("high", "low"), unique(person_id)
  ]
  ids_mortality <- mortality_risk[
    cohort == cohort_name & primary_mortality_interval == TRUE &
      wealth_time_eligible == TRUE & education3_fixed == "low" &
      wealth3 %chin% c("high", "low"), unique(person_id)
  ]
  data.table(
    cohort = cohort_name,
    stage = c(
      "Low education with high or low wealth at first valid wealth wave",
      "Contributed at least one eligible functional transition interval",
      "Contributed at least one supported mortality interval"
    ),
    people = c(length(ids_entry), length(ids_function), length(ids_mortality))
  )
}))

fwrite(focused_baseline, file.path(out_dir, "focused_loweducation_wealth_baseline_characteristics.csv"))
fwrite(focused_support, file.path(out_dir, "focused_loweducation_wealth_transition_support.csv"))
fwrite(ses_distribution, file.path(out_dir, "primary_education_wealth_distribution.csv"))
fwrite(focused_flow, file.path(out_dir, "focused_loweducation_wealth_flow.csv"))

cat("\nFocused baseline characteristics:\n")
print(focused_baseline)
cat("\nFocused transition support:\n")
print(focused_support)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_49_focused_baseline_support.txt"))
