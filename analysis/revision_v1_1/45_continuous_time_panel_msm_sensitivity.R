#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
project_lib <- file.path(revision_root, "07_r_library")
.libPaths(c(
  project_lib,
  "${RECOVERY_DIVIDE_AUX_DATA_ROOT}/14_charls_hrs_pulmonary_muscle_reserve/_rlib",
  .libPaths()
))

suppressPackageStartupMessages({
  library(data.table)
  library(msm)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 45_continuous_time_panel_msm_sensitivity.R CHARLS|ELSA|HRS|MHAS")
}
cohort_name <- toupper(args[[1L]])
stopifnot(cohort_name %in% c("CHARLS", "ELSA", "HRS", "MHAS"))

derived_dir <- file.path(revision_root, "02_derived")
out_dir <- file.path(revision_root, "03_outputs", "07_continuous_time_panel_msm")
log_dir <- file.path(revision_root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(
  log_dir,
  paste0("45_continuous_time_msm_", tolower(cohort_name), "_", stamp, ".log")
)
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Started: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
cat("Cohort: ", cohort_name, "\n", sep = "")
cat("msm version: ", as.character(packageVersion("msm")), "\n", sep = "")

panel <- as.data.table(readRDS(file.path(derived_dir, "individual_wave_panel_revision_v1_1.rds")))
mortality <- as.data.table(readRDS(file.path(derived_dir, "formal_mortality_riskset_revision_v1_1.rds")))
person_ses <- as.data.table(readRDS(file.path(derived_dir, "person_fixed_ses_revision_v1_1.rds")))

ses <- person_ses[
  cohort == cohort_name & education3_fixed == "low" & wealth3 %chin% c("high", "low"),
  .(
    person_id,
    wealth3,
    wealth_entry_wave,
    wealth_entry_time,
    wealth_entry_age,
    wealth_entry_household_id
  )
]
stopifnot(!anyDuplicated(ses$person_id), nrow(ses) > 0L)

state_codes <- c(I0 = 1L, D1 = 2L, R1 = 3L, D2 = 4L, DEAD = 5L)
state_labels <- names(state_codes)

living <- panel[
  cohort == cohort_name & person_id %chin% ses$person_id & in_wave == 1 &
    observed_state == TRUE & age >= 60 & history_state %chin% names(state_codes)[1:4],
  .(
    person_id,
    wave,
    observation_time = interview_time,
    age,
    female,
    history_state
  )
]
living <- ses[living, on = "person_id", nomatch = 0]
living <- living[observation_time >= wealth_entry_time]
living <- living[!is.na(female) & female %in% 0:1]
living[, `:=`(
  wealth_low = as.integer(wealth3 == "low"),
  age_center = age - 70,
  state = unname(state_codes[history_state]),
  obstype = 1L
)]
setorder(living, person_id, observation_time, wave)
living <- unique(living, by = c("person_id", "observation_time"))

# The primary discrete-time risk sets contain living observation gaps of at
# most four years. Split longer panel gaps into separate likelihood episodes
# so the continuous-time sensitivity uses the same observation support.
living[, prior_time := shift(observation_time), by = person_id]
living[, new_episode := as.integer(is.na(prior_time) | observation_time - prior_time > 4.000001), by = person_id]
living[, episode_number := cumsum(new_episode), by = person_id]
living[, episode_id := paste(cohort_name, person_id, episode_number, sep = "::")]
living[, c("prior_time", "new_episode") := NULL]

death <- mortality[
  cohort == cohort_name & primary_mortality_interval == TRUE & event_death == 1 &
    wealth_time_eligible == TRUE & education3_fixed == "low" & wealth3 %chin% c("high", "low"),
  .(
    person_id,
    origin_wave,
    death_time,
    origin_age,
    female,
    wealth3
  )
]
death <- death[!is.na(female) & female %in% 0:1]
death <- living[
  death,
  on = .(person_id, wave = origin_wave),
  nomatch = 0,
  .(
    person_id = x.person_id,
    wave = NA_integer_,
    observation_time = i.death_time,
    age = i.origin_age + (i.death_time - x.observation_time),
    female = i.female,
    history_state = "DEAD",
    wealth3 = i.wealth3,
    wealth_entry_wave = x.wealth_entry_wave,
    wealth_entry_time = x.wealth_entry_time,
    wealth_entry_age = x.wealth_entry_age,
    wealth_entry_household_id = x.wealth_entry_household_id,
    wealth_low = as.integer(i.wealth3 == "low"),
    age_center = i.origin_age + (i.death_time - x.observation_time) - 70,
    state = unname(state_codes[["DEAD"]]),
    obstype = 3L,
    episode_number = x.episode_number,
    episode_id = x.episode_id
  )
]

z <- rbindlist(list(living, death), use.names = TRUE, fill = TRUE)
setorder(z, episode_id, observation_time, obstype)
z <- unique(z, by = c("episode_id", "observation_time", "state"))
z[, records_in_episode := .N, by = episode_id]
z <- z[records_in_episode >= 2]
z[, analysis_time := observation_time - min(observation_time), by = episode_id]
z[, records_in_episode := NULL]

# A death is terminal. Any accidental later record would violate the source
# death audit and is blocked here rather than silently repaired.
terminal_qc <- z[, .(
  death_records = sum(state == state_codes[["DEAD"]]),
  last_is_death = if (any(state == state_codes[["DEAD"]])) state[.N] == state_codes[["DEAD"]] else TRUE
), by = episode_id]
if (any(terminal_qc$death_records > 1L) || any(!terminal_qc$last_is_death)) {
  stop("Death-terminal QC failed")
}

episode_flow <- data.table(
  cohort = cohort_name,
  eligible_low_education_high_low_wealth_people = nrow(ses),
  model_people = uniqueN(z$person_id),
  model_episodes = uniqueN(z$episode_id),
  model_records = nrow(z),
  living_records = sum(z$state != state_codes[["DEAD"]]),
  exact_death_records = sum(z$state == state_codes[["DEAD"]]),
  excluded_missing_sex_people = uniqueN(
    panel[cohort == cohort_name & person_id %chin% ses$person_id & is.na(female), person_id]
  ),
  long_gap_splits = sum(living$episode_number > 1L)
)
fwrite(
  episode_flow,
  file.path(out_dir, paste0("panel_flow_", tolower(cohort_name), ".csv"))
)

q_init <- matrix(0, 5, 5, dimnames = list(state_labels, state_labels))
q_init["I0", c("D1", "DEAD")] <- c(0.08, 0.02)
q_init["D1", c("R1", "DEAD")] <- c(0.20, 0.04)
q_init["R1", c("D2", "DEAD")] <- c(0.15, 0.03)
q_init["D2", c("R1", "DEAD")] <- c(0.15, 0.05)

cat("Fitting covariate-free initial model...\n")
fit0 <- msm(
  state ~ analysis_time,
  subject = episode_id,
  data = z,
  qmatrix = q_init,
  deathexact = state_codes[["DEAD"]],
  center = FALSE,
  control = list(maxit = 10000, reltol = 1e-8, fnscale = 4000)
)
if (fit0$opt$convergence != 0L) stop("Covariate-free msm did not converge")

q_start <- qmatrix.msm(fit0, covariates = "mean", ci = "none")
q_start <- as.matrix(q_start)

# Allowed intensity order, read across Q by row:
# I0-D1, I0-death, D1-R1, D1-death, R1-D2, R1-death,
# D2-R1, D2-death. Constraints reproduce the five process blocks used in
# the main model while retaining distinct baseline intensities by origin.
process_constraint <- c(1L, 2L, 3L, 4L, 5L, 4L, 3L, 4L)
cat("Fitting adjusted continuous-time panel model...\n")
fit <- msm(
  state ~ analysis_time,
  subject = episode_id,
  data = z,
  qmatrix = q_start,
  deathexact = state_codes[["DEAD"]],
  covariates = ~ wealth_low + age_center + female,
  constraint = list(
    wealth_low = process_constraint,
    age_center = process_constraint,
    female = process_constraint
  ),
  center = FALSE,
  control = list(maxit = 20000, reltol = 1e-8, fnscale = 4000)
)

if (fit$opt$convergence != 0L || !isTRUE(fit$foundse)) {
  stop("Adjusted continuous-time msm convergence/Hessian failure")
}

transition_order <- data.table(
  transition_index = 1:8,
  transition = c(
    "I0->D1", "I0->DEAD", "D1->R1", "D1->DEAD",
    "R1->D2", "R1->DEAD", "D2->R1", "D2->DEAD"
  ),
  process = c(
    "onset", "death_pre", "recovery", "death_post",
    "relapse", "death_post", "recovery", "death_post"
  )
)

haz <- as.data.table(hazard.msm(fit)$wealth_low, keep.rownames = "msm_transition")
setnames(haz, c("HR", "L", "U"), c("intensity_ratio", "ci_low", "ci_high"))
haz[, transition_index := .I]
haz <- transition_order[haz, on = "transition_index"]
haz[, `:=`(
  cohort = cohort_name,
  contrast = "low_vs_high_wealth_within_low_education",
  model = "continuous_time_panel_markov_age_sex_adjusted"
)]
setcolorder(haz, c(
  "cohort", "model", "contrast", "process", "transition",
  "intensity_ratio", "ci_low", "ci_high", "msm_transition"
))
fwrite(
  haz,
  file.path(out_dir, paste0("transition_intensity_ratios_", tolower(cohort_name), ".csv"))
)

process_haz <- unique(haz[, .(
  cohort,
  model,
  contrast,
  process,
  intensity_ratio,
  ci_low,
  ci_high
)], by = c("cohort", "process"))
fwrite(
  process_haz,
  file.path(out_dir, paste0("process_intensity_ratios_", tolower(cohort_name), ".csv"))
)

q_high <- as.matrix(qmatrix.msm(
  fit,
  covariates = list(wealth_low = 0, age_center = 0, female = mean(z$female)),
  ci = "none"
))
q_low <- as.matrix(qmatrix.msm(
  fit,
  covariates = list(wealth_low = 1, age_center = 0, female = mean(z$female)),
  ci = "none"
))
q_rows <- rbindlist(lapply(c("high", "low"), function(group_name) {
  q <- if (group_name == "high") q_high else q_low
  as.data.table(as.table(q))[, .(
    cohort = cohort_name,
    wealth = group_name,
    age = 70,
    female_standardisation = mean(z$female),
    origin = as.character(V1),
    destination = as.character(V2),
    intensity = as.numeric(N)
  )]
}))
fwrite(
  q_rows,
  file.path(out_dir, paste0("age70_intensity_matrices_", tolower(cohort_name), ".csv"))
)

# Endpoint-count calibration. The model probability matrix integrates over
# all permitted paths between observations; grouping only reduces repeated
# matrix-exponential calls and does not change observed counts.
pairs <- z[, .(
  origin = state,
  destination = shift(state, type = "lead"),
  interval = shift(analysis_time, type = "lead") - analysis_time,
  wealth_low,
  age_center,
  female
), by = episode_id][!is.na(destination) & interval > 0]
pairs[, `:=`(
  interval_group = round(interval, 6),
  age_center_group = round(age_center, 3)
)]
cal_groups <- pairs[, .N, by = .(
  wealth_low, female, age_center_group, interval_group, origin, destination
)]
prediction_groups <- unique(cal_groups[, .(
  wealth_low, female, age_center_group, interval_group, origin
)])

expected_parts <- vector("list", nrow(prediction_groups))
for (i in seq_len(nrow(prediction_groups))) {
  g <- prediction_groups[i]
  n_origin <- cal_groups[
    wealth_low == g$wealth_low & female == g$female &
      age_center_group == g$age_center_group & interval_group == g$interval_group &
      origin == g$origin,
    sum(N)
  ]
  p <- as.matrix(pmatrix.msm(
    fit,
    t = g$interval_group,
    covariates = list(
      wealth_low = g$wealth_low,
      age_center = g$age_center_group,
      female = g$female
    ),
    ci = "none"
  ))
  expected_parts[[i]] <- data.table(
    wealth_low = g$wealth_low,
    origin = g$origin,
    destination = 1:5,
    expected = n_origin * p[g$origin, ]
  )
}
expected <- rbindlist(expected_parts)[, .(expected = sum(expected)), by = .(
  wealth_low, origin, destination
)]
observed <- pairs[, .(observed = .N), by = .(wealth_low, origin, destination)]
calibration <- merge(
  CJ(wealth_low = 0:1, origin = 1:4, destination = 1:5),
  observed,
  by = c("wealth_low", "origin", "destination"),
  all.x = TRUE
)
calibration <- merge(
  calibration,
  expected,
  by = c("wealth_low", "origin", "destination"),
  all.x = TRUE
)
calibration[is.na(observed), observed := 0]
calibration[is.na(expected), expected := 0]
calibration[, `:=`(
  cohort = cohort_name,
  wealth = fifelse(wealth_low == 1, "low", "high"),
  origin_state = state_labels[origin],
  destination_state = state_labels[destination],
  observed_minus_expected = observed - expected,
  pearson_component = fifelse(expected > 0, (observed - expected)^2 / expected, NA_real_)
)]
setcolorder(calibration, c(
  "cohort", "wealth", "origin_state", "destination_state",
  "observed", "expected", "observed_minus_expected", "pearson_component"
))
fwrite(
  calibration,
  file.path(out_dir, paste0("endpoint_count_calibration_", tolower(cohort_name), ".csv"))
)

fit_qc <- data.table(
  cohort = cohort_name,
  model = "continuous_time_panel_markov_age_sex_adjusted",
  observations = nrow(z),
  transition_intervals = nrow(pairs),
  people = uniqueN(z$person_id),
  episodes = uniqueN(z$episode_id),
  exact_deaths = sum(z$state == state_codes[["DEAD"]]),
  parameters = length(fit$estimates),
  minus2loglik = fit$minus2loglik,
  optimizer_convergence = fit$opt$convergence,
  optimizer_message = if (is.null(fit$opt$message)) "" else fit$opt$message,
  hessian_standard_errors_found = fit$foundse,
  max_absolute_score = max(abs(fit$deriv)),
  calibration_pearson_sum = sum(calibration$pearson_component, na.rm = TRUE),
  calibration_absolute_count_error = sum(abs(calibration$observed_minus_expected)),
  calibration_total_observed = sum(calibration$observed),
  calibration_total_expected = sum(calibration$expected)
)
fwrite(
  fit_qc,
  file.path(out_dir, paste0("model_qc_", tolower(cohort_name), ".csv"))
)

writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, paste0("sessionInfo_45_continuous_time_", tolower(cohort_name), ".txt"))
)

cat("\nProcess intensity ratios (low vs high wealth):\n")
print(process_haz)
cat("\nModel QC:\n")
print(fit_qc)
cat("Completed: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
