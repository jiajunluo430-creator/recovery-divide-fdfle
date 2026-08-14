#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
derived_dir <- file.path(revision_root, "02_derived")
point_dir <- file.path(revision_root, "03_outputs", "02_revision_point")
out_dir <- file.path(revision_root, "03_outputs", "10_internal_survival_calibration")
log_dir <- file.path(revision_root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("48_internal_survival_calibration_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

mortality <- as.data.table(readRDS(file.path(
  derived_dir, "formal_mortality_riskset_revision_v1_1.rds"
)))
model_life <- fread(file.path(point_dir, "life_expectancy_point_estimates.csv"))
model_rates <- fread(file.path(point_dir, "standardised_annual_hazards.csv"))
model_initial <- fread(file.path(point_dir, "initial_difficulty_probability_age60.csv"))

cohort_order <- c("CHARLS", "ELSA", "HRS", "MHAS")
grid_step <- 0.01
age_grid <- seq(60, 100, by = grid_step)
living_states <- c("I0", "D1", "R1", "D2")
states <- c(living_states, "DEAD")
allowed <- list(
  I0 = c("D1", "DEAD"),
  D1 = c("R1", "DEAD"),
  R1 = c("D2", "DEAD"),
  D2 = c("R1", "DEAD")
)

model_survival_curve <- function(cohort_name, group_name) {
  z <- model_rates[
    cohort == cohort_name & module == "wealth_within_low_education" & ses == group_name
  ]
  lambda <- array(
    0, c(40L, length(living_states), length(states)),
    dimnames = list(age = as.character(60:99), origin = living_states, destination = states)
  )
  for (i in seq_len(nrow(z))) {
    lambda[as.character(z$age[[i]]), z$origin_state[[i]], z$destination[[i]]] <- z$annual_hazard[[i]]
  }
  p_diff <- model_initial[
    cohort == cohort_name & module == "wealth_within_low_education" & ses == group_name,
    difficulty_probability_age60
  ]
  if (length(p_diff) != 1L || !is.finite(p_diff)) stop("Missing model initial state")
  v <- setNames(c(1 - p_diff, p_diff, 0, 0, 0), states)
  ans <- data.table(age = 60L, survival = sum(v[living_states]))
  for (age_value in 60:99) {
    m <- matrix(0, length(states), length(states), dimnames = list(states, states))
    m["DEAD", "DEAD"] <- 1
    for (origin_name in names(allowed)) {
      destinations <- allowed[[origin_name]]
      hazard <- lambda[as.character(age_value), origin_name, destinations]
      total <- sum(hazard)
      move <- if (total > 0) 1 - exp(-total) else 0
      m[origin_name, origin_name] <- 1 - move
      if (total > 0) m[origin_name, destinations] <- move * hazard / total
    }
    v <- setNames(as.numeric(v %*% m), states)
    ans <- rbind(ans, data.table(age = age_value + 1L, survival = sum(v[living_states])))
  }
  ans
}

curve_rmst <- function(curve, end_age) {
  x <- curve[age <= end_age][order(age)]
  if (nrow(x) < 2L || tail(x$age, 1L) != end_age) stop("Incomplete model curve")
  sum((head(x$survival, -1L) + tail(x$survival, -1L)) / 2)
}

person_parts <- list()
curve_parts <- list()
calibration_parts <- list()

for (cohort_name in cohort_order) {
  z <- copy(mortality[
    cohort == cohort_name & primary_mortality_interval == TRUE &
      wealth_time_eligible == TRUE & education3_fixed == "low" &
      wealth3 %chin% c("high", "low") &
      is.finite(origin_age) & is.finite(interval_years) & interval_years > 0 &
      is.finite(origin_weight) & origin_weight > 0
  ])
  if (!nrow(z)) stop("No mortality calibration support for ", cohort_name)
  setorder(z, person_id, origin_time, endpoint_time)

  person <- z[, {
    death_rows <- which(event_death == 1L)
    exit_row <- if (length(death_rows)) death_rows[[1L]] else .N
    exit_age_raw <- origin_age[[exit_row]] + interval_years[[exit_row]]
    list(
      wealth3 = wealth3[[1L]],
      household_id = household_id[[1L]],
      entry_age = origin_age[[1L]],
      exit_age_raw = exit_age_raw,
      death_event_raw = as.integer(length(death_rows) > 0L),
      calibration_weight = origin_weight[[1L]],
      intervals = .N
    )
  }, by = person_id]
  person[, `:=`(
    exit_age = pmin(exit_age_raw, 100),
    death_event = as.integer(death_event_raw == 1L & exit_age_raw <= 100)
  )]
  person <- person[
    is.finite(entry_age) & is.finite(exit_age) & exit_age > entry_age & entry_age < 100 &
      is.finite(calibration_weight) & calibration_weight > 0
  ]
  if (!all(c("high", "low") %in% person$wealth3)) stop("Missing high/low calibration group")
  person[, cohort := cohort_name]
  person_parts[[cohort_name]] <- person

  fit <- survfit(
    Surv(entry_age, exit_age, death_event) ~ wealth3,
    data = person,
    weights = calibration_weight,
    conf.type = "log"
  )
  sm <- summary(fit, times = age_grid, extend = TRUE)
  curves <- data.table(
    cohort = cohort_name,
    wealth3 = sub("^wealth3=", "", as.character(sm$strata)),
    age = sm$time,
    survival = sm$surv,
    survival_ci_low = sm$lower,
    survival_ci_high = sm$upper,
    weighted_risk = sm$n.risk
  )
  curve_parts[[cohort_name]] <- curves

  for (group_name in c("high", "low")) {
    g <- curves[wealth3 == group_name][order(age)]
    if (nrow(g) != length(age_grid)) stop("Incomplete survival grid: ", cohort_name, " / ", group_name)
    km_tle_to_100 <- sum(head(g$survival, -1L)) * grid_step
    km_alive_at_100 <- tail(g$survival, 1L)
    g90 <- g[age <= 90]
    km_tle_to_90 <- sum(head(g90$survival, -1L)) * grid_step
    km_alive_at_90 <- g[age == 90, survival]
    fitted_curve <- model_survival_curve(cohort_name, group_name)
    model_tle_to_90 <- curve_rmst(fitted_curve, 90L)
    model_alive_at_90 <- fitted_curve[age == 90, survival]
    model_tle <- model_life[
      cohort == cohort_name & module == "wealth_within_low_education" &
        ses == group_name & estimand == "population_initialised" &
        metric == "total_life_expectancy", estimate
    ]
    model_alive <- model_life[
      cohort == cohort_name & module == "wealth_within_low_education" &
        ses == group_name & estimand == "population_initialised" &
        metric == "residual_alive_at_end", estimate
    ]
    if (length(model_tle) != 1L || length(model_alive) != 1L) stop("Missing model life metric")
    calibration_parts[[length(calibration_parts) + 1L]] <- data.table(
      cohort = cohort_name,
      wealth3 = group_name,
      people = person[wealth3 == group_name, .N],
      deaths_by_age_100 = person[wealth3 == group_name, sum(death_event)],
      people_at_risk_age_90 = person[wealth3 == group_name & entry_age <= 90 & exit_age >= 90, .N],
      people_at_risk_age_100 = person[wealth3 == group_name & entry_age <= 100 & exit_age >= 100, .N],
      weighted_delayed_entry_km_tle_to_age_90 = km_tle_to_90,
      multistate_model_tle_to_age_90 = model_tle_to_90,
      model_minus_km_years_to_age_90 = model_tle_to_90 - km_tle_to_90,
      km_alive_at_age_90 = km_alive_at_90,
      multistate_model_alive_at_age_90 = model_alive_at_90,
      weighted_delayed_entry_km_tle_to_age_100 = km_tle_to_100,
      multistate_model_tle_to_age_100 = model_tle,
      model_minus_km_years = model_tle - km_tle_to_100,
      km_alive_at_age_100 = km_alive_at_100,
      multistate_model_alive_at_age_100 = model_alive,
      model_minus_km_alive_probability = model_alive - km_alive_at_100
    )
    if (abs(curve_rmst(fitted_curve, 100L) - model_tle) > 1e-6 ||
        abs(fitted_curve[age == 100, survival] - model_alive) > 1e-6) {
      stop("Model survival reconstruction mismatch: ", cohort_name, " / ", group_name)
    }
  }
}

person_audit <- rbindlist(person_parts, use.names = TRUE, fill = TRUE)
curves <- rbindlist(curve_parts, use.names = TRUE, fill = TRUE)
calibration <- rbindlist(calibration_parts, use.names = TRUE, fill = TRUE)

fwrite(person_audit[, .(
  cohort, person_id, wealth3, household_id, entry_age, exit_age,
  death_event, calibration_weight, intervals
)], file.path(out_dir, "focused_lowedu_person_survival_audit_internal.csv"))
fwrite(curves, file.path(out_dir, "weighted_delayed_entry_survival_curves.csv"))
fwrite(calibration, file.path(out_dir, "model_vs_weighted_km_tle_calibration.csv"))
fwrite(data.table(
  calibration_estimand = "remaining total life years from age 60 through age 100",
  nonparametric_method = "survey-weighted delayed-entry Kaplan-Meier on attained-age scale",
  model_method = "population-initialised discrete-time multistate model",
  grid_step_years = grid_step,
  limitation = paste(
    "Internal cohort calibration, not an official national life-table comparison;",
    "wealth analyses begin at first valid wealth and therefore condition on survival to observed wealth entry."
  )
), file.path(out_dir, "calibration_contract.csv"))

cat("\nInternal TLE calibration through age 100:\n")
print(calibration)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_48_internal_survival_calibration.txt"))
