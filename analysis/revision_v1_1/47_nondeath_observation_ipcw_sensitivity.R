#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
  library(sandwich)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 47_nondeath_observation_ipcw_sensitivity.R CHARLS|ELSA|HRS|MHAS")
}
cohort_name <- toupper(args[[1L]])
stopifnot(cohort_name %in% c("CHARLS", "ELSA", "HRS", "MHAS"))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
derived_dir <- file.path(revision_root, "02_derived")
out_dir <- file.path(revision_root, "03_outputs", "09_nondeath_observation_ipcw")
log_dir <- file.path(revision_root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(
  log_dir,
  paste0("47_nondeath_ipcw_", tolower(cohort_name), "_", stamp, ".log")
)
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

panel <- as.data.table(readRDS(file.path(derived_dir, "individual_wave_panel_revision_v1_1.rds")))
person_ses <- as.data.table(readRDS(file.path(derived_dir, "person_fixed_ses_revision_v1_1.rds")))

ses <- person_ses[
  cohort == cohort_name & education3_fixed == "low" & wealth3 %chin% c("high", "low"),
  .(
    person_id,
    wealth3,
    wealth_entry_wave,
    wealth_entry_time,
    wealth_entry_household_id
  )
]
stopifnot(nrow(ses) > 0L, !anyDuplicated(ses$person_id))

origin <- panel[
  cohort == cohort_name & person_id %chin% ses$person_id & in_wave == 1 &
    observed_state == TRUE & age >= 60 & history_state %chin% c("I0", "D1", "R1", "D2") &
    !is.na(next_scheduled_wave) & !is.na(next_scheduled_time),
  .(
    person_id,
    origin_wave = wave,
    origin_time = interview_time,
    origin_age = age,
    female,
    proxy,
    respondent_weight,
    origin_state = history_state,
    next_scheduled_wave,
    next_scheduled_time,
    death_time
  )
]
origin <- ses[origin, on = "person_id", nomatch = 0]
origin <- origin[origin_time >= wealth_entry_time]
origin[, scheduled_interval := next_scheduled_time - origin_time]
origin <- origin[scheduled_interval >= 1 & scheduled_interval <= 4]

destination <- panel[
  cohort == cohort_name,
  .(
    person_id,
    destination_wave = wave,
    destination_in_wave = in_wave,
    destination_iwstat = iwstat,
    destination_observed_state = observed_state,
    destination_state = history_state,
    destination_difficulty = difficulty,
    destination_time = interview_time
  )
]
setkey(destination, person_id, destination_wave)
intervals <- destination[
  origin,
  on = .(person_id, destination_wave = next_scheduled_wave),
  nomatch = 0
]
intervals[, death_before_destination := !is.na(death_time) & death_time > origin_time & death_time <= next_scheduled_time]
intervals[, observed_next_function := as.integer(
  !death_before_destination & destination_in_wave == 1 & destination_observed_state == TRUE &
    destination_state %chin% c("I0", "D1", "R1", "D2")
)]
intervals[, nondeath_risk := !death_before_destination]
intervals[, `:=`(
  ses = factor(wealth3, levels = c("high", "low")),
  origin_factor = factor(origin_state, levels = c("I0", "D1", "R1", "D2")),
  sex_factor = factor(fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown"))),
  proxy_factor = factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self"))),
  period_factor = factor(origin_wave),
  age_model = pmin(pmax(origin_age, 60), 95),
  cluster_id = factor(wealth_entry_household_id)
)]
intervals <- intervals[
  !is.na(ses) & is.finite(respondent_weight) & respondent_weight > 0 &
    !is.na(cluster_id) & is.finite(scheduled_interval) & scheduled_interval > 0
]
for (v in c("ses", "origin_factor", "sex_factor", "proxy_factor", "period_factor")) {
  set(intervals, j = v, value = droplevels(intervals[[v]]))
}
intervals[, survey_weight := respondent_weight / mean(respondent_weight)]

at_risk <- intervals[nondeath_risk == TRUE]
age_range <- range(at_risk$age_model, na.rm = TRUE)
knots <- c(70, 80, 90)
knots <- knots[knots > age_range[[1L]] & knots < age_range[[2L]]]
age_term <- if (length(knots)) {
  paste0(
    "ns(age_model, knots=c(", paste(knots, collapse = ","),
    "), Boundary.knots=c(60,95))"
  )
} else {
  "age_model"
}

den_terms <- c("ses", age_term, "scheduled_interval")
for (v in c("origin_factor", "sex_factor", "period_factor", "proxy_factor")) {
  if (nlevels(at_risk[[v]]) > 1L) den_terms <- c(den_terms, v)
}
num_terms <- c("ses")
if (nlevels(at_risk$period_factor) > 1L) num_terms <- c(num_terms, "period_factor")

den_fit <- suppressWarnings(glm(
  as.formula(paste("observed_next_function ~", paste(den_terms, collapse = " + "))),
  data = at_risk,
  family = binomial(),
  weights = survey_weight,
  control = glm.control(maxit = 100, epsilon = 1e-9)
))
num_fit <- suppressWarnings(glm(
  as.formula(paste("observed_next_function ~", paste(num_terms, collapse = " + "))),
  data = at_risk,
  family = binomial(),
  weights = survey_weight,
  control = glm.control(maxit = 100, epsilon = 1e-9)
))
if (!isTRUE(den_fit$converged) || !isTRUE(num_fit$converged)) stop("Observation models did not converge")

at_risk[, `:=`(
  p_observed_denominator = pmin(pmax(as.numeric(predict(den_fit, type = "response")), 0.01), 0.995),
  p_observed_numerator = pmin(pmax(as.numeric(predict(num_fit, type = "response")), 0.01), 0.995)
)]
at_risk[, interval_ipcw_raw := p_observed_numerator / p_observed_denominator]
observed_weight_limits <- quantile(
  at_risk[observed_next_function == 1, interval_ipcw_raw],
  c(0.01, 0.99),
  na.rm = TRUE,
  names = FALSE,
  type = 6
)
at_risk[, interval_ipcw := pmin(pmax(interval_ipcw_raw, observed_weight_limits[[1L]]), observed_weight_limits[[2L]])]

retention_summary <- at_risk[, .(
  nondeath_intervals = .N,
  observed_next_function_n = sum(observed_next_function),
  observed_next_function_percent = 100 * mean(observed_next_function),
  interval_mean = mean(scheduled_interval),
  interval_sd = sd(scheduled_interval),
  interval_q25 = quantile(scheduled_interval, 0.25),
  interval_median = median(scheduled_interval),
  interval_q75 = quantile(scheduled_interval, 0.75),
  raw_ipcw_mean_observed = mean(interval_ipcw_raw[observed_next_function == 1]),
  truncated_ipcw_mean_observed = mean(interval_ipcw[observed_next_function == 1]),
  truncated_ipcw_max_observed = max(interval_ipcw[observed_next_function == 1])
), by = wealth3]
retention_summary[, `:=`(
  cohort = cohort_name,
  ipcw_lower_truncation = observed_weight_limits[[1L]],
  ipcw_upper_truncation = observed_weight_limits[[2L]],
  unknown_vital_status_treated_as_death = FALSE
)]
fwrite(
  retention_summary,
  file.path(out_dir, paste0("retention_interval_summary_", tolower(cohort_name), ".csv"))
)

observed <- at_risk[observed_next_function == 1]
observed[, `:=`(
  event_onset = as.integer(origin_state == "I0" & destination_state == "D1"),
  event_recovery = as.integer(origin_state %chin% c("D1", "D2") & destination_state == "R1"),
  event_relapse = as.integer(origin_state == "R1" & destination_state == "D2")
)]

process_specs <- list(
  onset = list(origins = "I0", event = "event_onset"),
  recovery = list(origins = c("D1", "D2"), event = "event_recovery"),
  relapse = list(origins = "R1", event = "event_relapse")
)

transition_formula <- function(z) {
  terms <- c("ses", age_term)
  for (v in c("origin_factor", "sex_factor", "period_factor", "proxy_factor")) {
    if (nlevels(z[[v]]) > 1L) terms <- c(terms, v)
  }
  as.formula(paste("event ~", paste(terms, collapse = " + "), "+ offset(log(scheduled_interval))"))
}

fit_transition <- function(z, use_ipcw) {
  z <- copy(z)
  z[, event := as.integer(get(process_specs[[process_name]]$event))]
  z[, fit_weight_raw := survey_weight * if (use_ipcw) interval_ipcw else 1]
  z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
  for (v in c("ses", "origin_factor", "sex_factor", "period_factor", "proxy_factor")) {
    set(z, j = v, value = droplevels(z[[v]]))
  }
  fit <- suppressWarnings(glm(
    transition_formula(z),
    data = z,
    family = binomial(link = "cloglog"),
    weights = fit_weight,
    control = glm.control(maxit = 100, epsilon = 1e-9)
  ))
  if (!isTRUE(fit$converged) || !"seslow" %in% names(coef(fit))) stop("IPCW transition fit failed")
  vc <- vcovCL(fit, cluster = z$cluster_id, type = "HC0")
  beta <- coef(fit)[["seslow"]]
  se <- sqrt(vc["seslow", "seslow"])
  data.table(
    adjustment = if (use_ipcw) "interval_ipcw_1_99_truncated" else "adjacent_complete_case",
    intensity_ratio = exp(beta),
    ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se),
    robust_se = se,
    intervals = nrow(z),
    people = uniqueN(z$person_id),
    households = uniqueN(z$wealth_entry_household_id),
    events = sum(z$event),
    converged = fit$converged
  )
}

result_parts <- list()
for (process_name in names(process_specs)) {
  spec <- process_specs[[process_name]]
  z <- observed[origin_state %chin% spec$origins]
  z[, origin_factor := droplevels(factor(origin_state, levels = spec$origins))]
  ans <- rbindlist(list(fit_transition(z, FALSE), fit_transition(z, TRUE)))
  ans[, `:=`(
    cohort = cohort_name,
    process = process_name,
    contrast = "low_vs_high_wealth_within_low_education"
  )]
  result_parts[[process_name]] <- ans
}
results <- rbindlist(result_parts)
setcolorder(results, c(
  "cohort", "process", "contrast", "adjustment", "intensity_ratio",
  "ci_low", "ci_high", "robust_se", "intervals", "people", "households",
  "events", "converged"
))
fwrite(
  results,
  file.path(out_dir, paste0("ipcw_transition_results_", tolower(cohort_name), ".csv"))
)

observation_qc <- data.table(
  cohort = cohort_name,
  eligible_origin_intervals = nrow(intervals),
  verified_death_before_next = sum(intervals$death_before_destination),
  nondeath_observation_risk = nrow(at_risk),
  observed_next_function = sum(at_risk$observed_next_function),
  nondeath_not_observed = sum(at_risk$observed_next_function == 0),
  denominator_model_converged = den_fit$converged,
  numerator_model_converged = num_fit$converged,
  ipcw_lower_truncation = observed_weight_limits[[1L]],
  ipcw_upper_truncation = observed_weight_limits[[2L]]
)
fwrite(
  observation_qc,
  file.path(out_dir, paste0("observation_model_qc_", tolower(cohort_name), ".csv"))
)

writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, paste0("sessionInfo_47_ipcw_", tolower(cohort_name), ".txt"))
)

cat("Retention and interval summary:\n")
print(retention_summary)
cat("\nIPCW transition sensitivity:\n")
print(results)
cat("Completed: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
