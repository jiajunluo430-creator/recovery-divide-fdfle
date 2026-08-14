#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: 40_bootstrap_low_education_household.R COHORT REPLICATES PHASE [START] [END]")
}
cohort_name <- toupper(args[[1L]])
replicates <- as.integer(args[[2L]])
phase <- tolower(args[[3L]])
replicate_start <- if (length(args) >= 4L) as.integer(args[[4L]]) else 1L
replicate_end <- if (length(args) >= 5L) as.integer(args[[5L]]) else replicates
stopifnot(cohort_name %in% c("CHARLS", "ELSA", "HRS", "MHAS"))
stopifnot(phase %in% c("test", "preview", "final"))
stopifnot(replicates >= 1L, replicate_start >= 1L, replicate_end <= replicates, replicate_start <= replicate_end)

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
base_derived_dir <- file.path(revision_root, "02_derived")
sensitivity <- tolower(Sys.getenv("D4_REV_SENSITIVITY", unset = "primary"))
if (!sensitivity %in% c("primary", "adl_only", "at_least_two")) {
  stop("Unsupported D4_REV_SENSITIVITY for bootstrap: ", sensitivity)
}
analysis_module <- tolower(Sys.getenv(
  "D4_REV_ANALYSIS_MODULE",
  unset = "wealth_within_low_education"
))
resampling_unit <- tolower(Sys.getenv("D4_REV_RESAMPLING_UNIT", unset = "household"))
if (!resampling_unit %in% c("household", "stratified_psu")) {
  stop("Unsupported D4_REV_RESAMPLING_UNIT: ", resampling_unit)
}
allowed_modules <- c(
  "primary_education", "primary_wealth", "wealth_within_low_education"
)
if (!analysis_module %in% allowed_modules) {
  stop("Unsupported D4_REV_ANALYSIS_MODULE: ", analysis_module)
}
if (sensitivity != "primary" && analysis_module != "wealth_within_low_education") {
  stop("Threshold bootstrap sensitivities are defined only for wealth_within_low_education")
}
if (resampling_unit == "stratified_psu" &&
    (sensitivity != "primary" || analysis_module != "wealth_within_low_education" ||
      !cohort_name %in% c("ELSA", "HRS"))) {
  stop("Stratified PSU bootstrap is restricted to the primary focused analysis in ELSA/HRS")
}
derived_dir <- if (sensitivity == "primary") {
  base_derived_dir
} else {
  file.path(base_derived_dir, "sensitivity", sensitivity)
}
output_tag <- if (sensitivity == "primary") phase else paste0(sensitivity, "_", phase)
output_prefix <- if (resampling_unit == "stratified_psu") {
  "03_lowedu_stratified_psu_bootstrap_"
} else {
  switch(
    analysis_module,
    primary_education = "03_primary_education_household_bootstrap_",
    primary_wealth = "03_primary_wealth_household_bootstrap_",
    wealth_within_low_education = "03_lowedu_household_bootstrap_"
  )
}
out_dir <- file.path(revision_root, "03_outputs", paste0(output_prefix, output_tag))
log_dir <- file.path(revision_root, "06_logs")
for (d in c(out_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

result_path <- file.path(out_dir, paste0("bootstrap_", tolower(cohort_name), "_metrics.csv"))
qc_path <- file.path(out_dir, paste0("bootstrap_", tolower(cohort_name), "_qc.csv"))
cut_path <- file.path(out_dir, paste0("bootstrap_", tolower(cohort_name), "_wealth_cutpoints.csv"))

function_all <- as.data.table(readRDS(file.path(derived_dir, "formal_function_transition_riskset_revision_v1_1.rds")))
mortality_all <- as.data.table(readRDS(file.path(derived_dir, "formal_mortality_riskset_revision_v1_1.rds")))
initial_all <- as.data.table(readRDS(file.path(derived_dir, "initial_state_data_revision_v1_1.rds")))
person_all <- as.data.table(readRDS(file.path(base_derived_dir, "person_fixed_ses_revision_v1_1.rds")))

if (analysis_module == "primary_education") {
  function_base <- function_all[cohort == cohort_name & education_eligible == TRUE]
  mortality_base <- mortality_all[
    cohort == cohort_name & primary_mortality_interval == TRUE & education_eligible == TRUE
  ]
  initial_base <- initial_all[cohort == cohort_name & exposure == "education"]
  person_ses <- unique(initial_base[, .(
      person_id,
      cluster_household_id = household_id,
      cluster_design_stratum = design_stratum,
      cluster_design_psu = design_psu,
      fixed_ses = level
  )], by = "person_id")
} else {
  low_education_only <- analysis_module == "wealth_within_low_education"
  function_base <- function_all[cohort == cohort_name & wealth_time_eligible == TRUE]
  mortality_base <- mortality_all[
    cohort == cohort_name & primary_mortality_interval == TRUE & wealth_time_eligible == TRUE
  ]
  if (low_education_only) {
    function_base <- function_base[education3_fixed == "low"]
    mortality_base <- mortality_base[education3_fixed == "low"]
  }
  initial_base <- initial_all[cohort == cohort_name & exposure == "wealth"]
  if (low_education_only) {
    low_ids <- unique(c(function_base$person_id, mortality_base$person_id))
    initial_base <- initial_base[person_id %in% low_ids]
  }
  person_ses <- person_all[
    cohort == cohort_name & !is.na(wealth_entry_value) & !is.na(wealth_entry_household_id),
    .(
      person_id,
      cluster_household_id = wealth_entry_household_id,
      cluster_design_stratum = wealth_entry_design_stratum,
      cluster_design_psu = wealth_entry_design_psu,
      wealth_entry_wave,
      wealth_entry_value,
      wealth_entry_weight
    )
  ]
}

design_missing_people_excluded <- 0L
if (resampling_unit == "stratified_psu") {
  design_missing_people_excluded <- person_ses[
    is.na(cluster_design_stratum) | is.na(cluster_design_psu), .N
  ]
  person_ses <- person_ses[!is.na(cluster_design_stratum) & !is.na(cluster_design_psu)]
  complete_design_ids <- person_ses$person_id
  function_base <- function_base[person_id %in% complete_design_ids]
  mortality_base <- mortality_base[person_id %in% complete_design_ids]
  initial_base <- initial_base[person_id %in% complete_design_ids]
}

stopifnot(nrow(function_base) > 0L, nrow(mortality_base) > 0L, nrow(initial_base) > 0L, nrow(person_ses) > 0L)
stopifnot(!anyDuplicated(person_ses$person_id))

if (resampling_unit == "household") {
  person_ses[, cluster_id := cluster_household_id]
  cluster_frame <- unique(person_ses[, .(cluster_id, design_stratum = NA_character_)])
} else {
  person_ses[, cluster_id := paste(cluster_design_stratum, cluster_design_psu, sep = "::")]
  cluster_frame <- unique(person_ses[, .(
    cluster_id,
    design_stratum = cluster_design_stratum
  )])
}
setorder(cluster_frame, design_stratum, cluster_id)
cluster_frame[, cluster_index := .I]
person_ses[, cluster_index := cluster_frame$cluster_index[match(cluster_id, cluster_frame$cluster_id)]]
person_index <- person_ses[, .(person_id, cluster_index)]
function_base[, cluster_index := person_index$cluster_index[match(person_id, person_index$person_id)]]
mortality_base[, cluster_index := person_index$cluster_index[match(person_id, person_index$person_id)]]
initial_base[, cluster_index := person_index$cluster_index[match(person_id, person_index$person_id)]]
stopifnot(!anyNA(function_base$cluster_index), !anyNA(mortality_base$cluster_index), !anyNA(initial_base$cluster_index))

ses_levels <- c("high", "middle", "low")
ages <- 60:99
living_states <- c("I0", "D1", "R1", "D2")
states <- c(living_states, "DEAD")
allowed <- list(
  I0 = c("D1", "DEAD"),
  D1 = c("R1", "DEAD"),
  R1 = c("D2", "DEAD"),
  D2 = c("R1", "DEAD")
)
process_specs <- list(
  onset = list(source = "function", origins = "I0", destination = "D1", event = "event_onset"),
  recovery = list(source = "function", origins = c("D1", "D2"), destination = "R1", event = "event_recovery"),
  relapse = list(source = "function", origins = "R1", destination = "D2", event = "event_relapse"),
  death_pre = list(source = "mortality", origins = "I0", destination = "DEAD", event = "event_death"),
  death_post = list(source = "mortality", origins = c("D1", "R1", "D2"), destination = "DEAD", event = "event_death")
)
block_pairs <- list(
  onset = list(c("I0", "D1")),
  recovery = list(c("D1", "R1"), c("D2", "R1")),
  relapse = list(c("R1", "D2")),
  post_difficulty_mortality = list(c("D1", "DEAD"), c("R1", "DEAD"), c("D2", "DEAD")),
  pre_difficulty_mortality = list(c("I0", "DEAD"))
)

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

weighted_mean_finite <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

bootstrap_ses <- function(cluster_mult, replicate_id) {
  z <- copy(person_ses)
  z[, boot_mult := cluster_mult[cluster_index]]
  entry <- z[boot_mult > 0]
  if (analysis_module == "primary_education") {
    cuts <- data.table(
      cohort = cohort_name,
      replicate = replicate_id,
      wealth_entry_wave = NA_integer_,
      q33 = NA_real_,
      q67 = NA_real_,
      resampled_household_frequency = sum(cluster_mult),
      sampled_people_frequency = sum(entry$boot_mult),
      unique_sampled_people = nrow(entry),
      note = "education categories fixed; no replicate-specific wealth cutpoints"
    )
    return(list(map = entry[, .(person_id, ses_boot = fixed_ses)], cuts = cuts))
  }
  cuts <- entry[, {
    weighted_ok <- any(is.finite(wealth_entry_weight) & wealth_entry_weight > 0)
    ww <- if (weighted_ok) wealth_entry_weight * boot_mult else boot_mult
    list(
      q33 = weighted_cut(wealth_entry_value, ww, 1 / 3),
      q67 = weighted_cut(wealth_entry_value, ww, 2 / 3),
      resampled_household_frequency = sum(unique(boot_mult), na.rm = TRUE),
      sampled_people_frequency = sum(boot_mult),
      unique_sampled_people = .N
    )
  }, by = wealth_entry_wave]
  entry <- cuts[entry, on = "wealth_entry_wave"]
  entry[, ses_boot := fifelse(
    wealth_entry_value <= q33, "low",
    fifelse(wealth_entry_value <= q67, "middle", "high")
  )]
  cuts[, `:=`(cohort = cohort_name, replicate = replicate_id)]
  list(map = entry[, .(person_id, ses_boot)], cuts = cuts)
}

make_transition_formula <- function(z) {
  terms <- c("ses", "ns(age_model, knots = c(70, 80, 90), Boundary.knots = c(60, 95))")
  if (nlevels(z$origin_factor) > 1L) terms <- c(terms, "origin_factor")
  if (nlevels(z$sex_factor) > 1L) terms <- c(terms, "sex_factor")
  if (nlevels(z$period_factor) > 1L) terms <- c(terms, "period_factor")
  if (nlevels(z$proxy_factor) > 1L) terms <- c(terms, "proxy_factor")
  as.formula(paste("event ~", paste(terms, collapse = " + "), "+ offset(log(interval_years))"))
}

prepare_transition <- function(fboot, mboot, process_name) {
  spec <- process_specs[[process_name]]
  z <- if (spec$source == "function") copy(fboot[origin_state %in% spec$origins]) else copy(mboot[origin_state %in% spec$origins])
  z[, `:=`(
    event = as.integer(get(spec$event)),
    ses = factor(ses_boot, levels = ses_levels),
    sex_factor = factor(fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown"))),
    proxy_factor = factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self"))),
    period_factor = factor(origin_wave),
    origin_factor = factor(origin_state, levels = spec$origins),
    age_model = pmin(pmax(origin_age, 60), 95)
  )]
  z <- z[
    boot_mult > 0 & !is.na(ses) & is.finite(interval_years) & interval_years > 0 &
      is.finite(base_weight) & base_weight > 0 & !is.na(age_model)
  ]
  z[, fit_weight_raw := base_weight * boot_mult]
  z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
  for (v in c("ses", "sex_factor", "proxy_factor", "period_factor", "origin_factor")) set(z, j = v, value = droplevels(z[[v]]))
  z
}

fit_transition_rates <- function(fboot, mboot, process_name) {
  spec <- process_specs[[process_name]]
  z <- prepare_transition(fboot, mboot, process_name)
  if (!all(c("high", "low") %in% levels(z$ses))) stop("high/low wealth absent")
  fit <- suppressWarnings(glm(
    make_transition_formula(z), data = z, family = binomial(link = "cloglog"),
    weights = fit_weight, control = glm.control(maxit = 100, epsilon = 1e-9)
  ))
  if (!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) stop("transition model nonconvergent/nonfinite")

  low_term <- "seslow"
  if (!low_term %in% names(coef(fit))) stop("seslow coefficient absent")
  hr <- exp(coef(fit)[[low_term]])
  rate_parts <- list()
  for (origin_name in spec$origins) {
    ref <- z[origin_state == origin_name, .(fit_weight = sum(fit_weight)), by = .(sex_factor, period_factor, proxy_factor, origin_factor)]
    if (!nrow(ref)) stop("origin standardisation set empty")
    ref[, `:=`(
      sex_factor = factor(sex_factor, levels = levels(z$sex_factor)),
      period_factor = factor(period_factor, levels = levels(z$period_factor)),
      proxy_factor = factor(proxy_factor, levels = levels(z$proxy_factor)),
      origin_factor = factor(origin_factor, levels = levels(z$origin_factor))
    )]
    for (ses_name in c("high", "low")) {
      for (age_value in ages) {
        nd <- copy(ref)
        nd[, `:=`(ses = factor(ses_name, levels = levels(z$ses)), age_model = min(age_value, 95), interval_years = 1)]
        hazard <- exp(as.numeric(predict(fit, newdata = nd, type = "link")))
        rate_parts[[length(rate_parts) + 1L]] <- data.table(
          ses = ses_name, age = age_value, process = process_name,
          origin_state = origin_name, destination = spec$destination,
          annual_hazard = weighted_mean_finite(hazard, nd$fit_weight)
        )
      }
    }
  }
  list(rates = rbindlist(rate_parts), hazard_ratio_low_vs_high = hr)
}

make_initial_formula <- function(z, include_period = TRUE, include_proxy = TRUE) {
  terms <- c("ses", "ns(age_model, knots = c(70, 80, 90), Boundary.knots = c(60, 95))")
  if (nlevels(z$sex_factor) > 1L) terms <- c(terms, "sex_factor")
  if (include_period && nlevels(z$period_factor) > 1L) terms <- c(terms, "period_factor")
  if (include_proxy && nlevels(z$proxy_factor) > 1L) terms <- c(terms, "proxy_factor")
  as.formula(paste("entry_difficulty ~", paste(terms, collapse = " + ")))
}

fit_initial_state <- function(iboot) {
  z <- copy(iboot)
  z[, `:=`(
    ses = factor(ses_boot, levels = ses_levels),
    sex_factor = factor(fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown"))),
    proxy_factor = factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self"))),
    period_factor = factor(entry_wave),
    age_model = pmin(pmax(entry_age, 60), 95),
    base_weight = entry_weight
  )]
  z <- z[boot_mult > 0 & !is.na(ses) & !is.na(entry_difficulty) & is.finite(base_weight) & base_weight > 0]
  z[, fit_weight_raw := base_weight * boot_mult]
  z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
  for (v in c("ses", "sex_factor", "proxy_factor", "period_factor")) set(z, j = v, value = droplevels(z[[v]]))
  if (!all(c("high", "low") %in% levels(z$ses))) stop("initial high/low wealth absent")
  include_period <- TRUE
  include_proxy <- TRUE
  fit <- suppressWarnings(glm(make_initial_formula(z, include_period, include_proxy), data = z, family = binomial(), weights = fit_weight))
  if (any(!is.finite(coef(fit))) && nlevels(z$proxy_factor) > 1L) {
    include_proxy <- FALSE
    fit <- suppressWarnings(glm(make_initial_formula(z, include_period, include_proxy), data = z, family = binomial(), weights = fit_weight))
  }
  if (any(!is.finite(coef(fit))) && nlevels(z$period_factor) > 1L) {
    include_period <- FALSE
    fit <- suppressWarnings(glm(make_initial_formula(z, include_period, include_proxy), data = z, family = binomial(), weights = fit_weight))
  }
  if (!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) stop("initial model nonconvergent/nonfinite")
  ref <- z[, .(fit_weight = sum(fit_weight)), by = .(sex_factor, period_factor, proxy_factor)]
  ref[, `:=`(
    sex_factor = factor(sex_factor, levels = levels(z$sex_factor)),
    period_factor = factor(period_factor, levels = levels(z$period_factor)),
    proxy_factor = factor(proxy_factor, levels = levels(z$proxy_factor))
  )]
  out <- setNames(numeric(2L), c("high", "low"))
  for (ses_name in names(out)) {
    nd <- copy(ref)
    nd[, `:=`(ses = factor(ses_name, levels = levels(z$ses)), age_model = 60)]
    out[[ses_name]] <- weighted_mean_finite(as.numeric(predict(fit, newdata = nd, type = "response")), nd$fit_weight)
  }
  out
}

lambda_array <- function(z) {
  a <- array(NA_real_, c(length(ages), length(living_states), length(states)), dimnames = list(as.character(ages), living_states, states))
  for (i in seq_len(nrow(z))) a[as.character(z$age[[i]]), z$origin_state[[i]], z$destination[[i]]] <- z$annual_hazard[[i]]
  for (origin_name in names(allowed)) for (destination_name in allowed[[origin_name]]) {
    if (any(!is.finite(a[, origin_name, destination_name]))) stop("missing lambda")
  }
  a
}

transition_matrix <- function(lambda_age) {
  m <- matrix(0, length(states), length(states), dimnames = list(states, states))
  m["DEAD", "DEAD"] <- 1
  for (origin_name in names(allowed)) {
    destinations <- allowed[[origin_name]]
    lambda <- lambda_age[origin_name, destinations]
    total <- sum(lambda)
    move <- if (total > 0) 1 - exp(-total) else 0
    m[origin_name, origin_name] <- 1 - move
    if (total > 0) m[origin_name, destinations] <- move * lambda / total
  }
  m
}

life_metrics <- function(lambda, initial_vector) {
  v <- initial_vector[states]
  occupancy <- setNames(rep(0, length(states)), states)
  for (age_value in ages) {
    next_v <- as.numeric(v %*% transition_matrix(lambda[as.character(age_value), , , drop = TRUE]))
    names(next_v) <- states
    occupancy <- occupancy + (v + next_v) / 2
    v <- next_v
  }
  tle <- sum(occupancy[living_states])
  fdfle <- sum(occupancy[c("I0", "R1")])
  difficulty <- sum(occupancy[c("D1", "D2")])
  c(tle = tle, fdfle = fdfle, difficulty_years = difficulty, fdfle_percent = 100 * fdfle / tle)
}

replace_hazard_block <- function(base, donor, block_name) {
  out <- base
  for (pair in block_pairs[[block_name]]) out[, pair[[1]], pair[[2]]] <- donor[, pair[[1]], pair[[2]]]
  out
}

all_subsets <- function(x) {
  unlist(lapply(0:length(x), function(k) if (k == 0L) list(character()) else combn(x, k, simplify = FALSE)), recursive = FALSE)
}

subset_key <- function(x) if (!length(x)) "<empty>" else paste(sort(x), collapse = "+")

shapley_decompose <- function(high_lambda, low_lambda, high_initial, low_initial, include_initial) {
  hazard_blocks <- names(block_pairs)
  blocks <- c(if (include_initial) "initial_state", hazard_blocks)
  n_blocks <- length(blocks)
  values <- numeric()
  for (s in all_subsets(blocks)) {
    lambda <- high_lambda
    for (block_name in intersect(s, hazard_blocks)) lambda <- replace_hazard_block(lambda, low_lambda, block_name)
    initial <- if (include_initial && "initial_state" %in% s) low_initial else high_initial
    values[[subset_key(s)]] <- life_metrics(lambda, initial)[["fdfle"]]
  }
  contribution <- setNames(numeric(n_blocks), blocks)
  for (block_name in blocks) {
    for (s in all_subsets(setdiff(blocks, block_name))) {
      k <- length(s)
      weight <- factorial(k) * factorial(n_blocks - k - 1L) / factorial(n_blocks)
      contribution[[block_name]] <- contribution[[block_name]] + weight *
        (values[[subset_key(c(s, block_name))]] - values[[subset_key(s)]])
    }
  }
  gap <- values[[subset_key(blocks)]] - values[[subset_key(character())]]
  if (abs(sum(contribution) - gap) > 0.01) stop("Shapley closure failure")
  list(gap = gap, contribution = contribution)
}

calculate_metrics <- function(fboot, mboot, iboot) {
  transition_results <- lapply(names(process_specs), function(process_name) fit_transition_rates(fboot, mboot, process_name))
  rates <- rbindlist(lapply(transition_results, `[[`, "rates"))
  names(transition_results) <- names(process_specs)
  initial_probability <- fit_initial_state(iboot)
  lambdas <- list(
    high = lambda_array(rates[ses == "high"]),
    low = lambda_array(rates[ses == "low"])
  )
  initials <- list(
    high = setNames(c(1 - initial_probability[["high"]], initial_probability[["high"]], 0, 0, 0), states),
    low = setNames(c(1 - initial_probability[["low"]], initial_probability[["low"]], 0, 0, 0), states)
  )
  all_i0 <- setNames(c(1, 0, 0, 0, 0), states)
  pop_dec <- shapley_decompose(lambdas$high, lambdas$low, initials$high, initials$low, TRUE)
  con_dec <- shapley_decompose(lambdas$high, lambdas$low, all_i0, all_i0, FALSE)
  pop_high <- life_metrics(lambdas$high, initials$high)
  pop_low <- life_metrics(lambdas$low, initials$low)
  con_high <- life_metrics(lambdas$high, all_i0)
  con_low <- life_metrics(lambdas$low, all_i0)

  values <- c(
    initial_difficulty_high = initial_probability[["high"]],
    initial_difficulty_low = initial_probability[["low"]],
    population_high_tle = pop_high[["tle"]], population_low_tle = pop_low[["tle"]],
    population_high_fdfle = pop_high[["fdfle"]], population_low_fdfle = pop_low[["fdfle"]],
    population_high_difficulty_years = pop_high[["difficulty_years"]], population_low_difficulty_years = pop_low[["difficulty_years"]],
    population_high_fdfle_percent = pop_high[["fdfle_percent"]], population_low_fdfle_percent = pop_low[["fdfle_percent"]],
    population_fdfle_gap = pop_dec$gap,
    conditional_high_fdfle = con_high[["fdfle"]], conditional_low_fdfle = con_low[["fdfle"]],
    conditional_fdfle_gap = con_dec$gap,
    population_recovery_relapse = pop_dec$contribution[["recovery"]] + pop_dec$contribution[["relapse"]],
    conditional_recovery_relapse = con_dec$contribution[["recovery"]] + con_dec$contribution[["relapse"]]
  )
  for (block_name in names(pop_dec$contribution)) values[[paste0("population_contribution_", block_name)]] <- pop_dec$contribution[[block_name]]
  for (block_name in names(con_dec$contribution)) values[[paste0("conditional_contribution_", block_name)]] <- con_dec$contribution[[block_name]]
  for (process_name in names(transition_results)) values[[paste0("hazard_ratio_", process_name)]] <- transition_results[[process_name]]$hazard_ratio_low_vs_high
  data.table(metric = names(values), estimate = as.numeric(values))
}

safe_fwrite <- function(x, path, attempts = 10L) {
  last_error <- NULL
  for (i in seq_len(attempts)) {
    ok <- tryCatch({ fwrite(x, path); TRUE }, error = function(e) { last_error <<- e; FALSE })
    if (ok) return(invisible(TRUE))
    Sys.sleep(1)
  }
  stop("Failed to write ", path, ": ", conditionMessage(last_error))
}

existing_results <- if (file.exists(result_path)) fread(result_path) else data.table()
existing_qc <- if (file.exists(qc_path)) fread(qc_path) else data.table()
existing_cuts <- if (file.exists(cut_path)) fread(cut_path) else data.table()
completed <- if (nrow(existing_qc)) unique(existing_qc[status == "valid", replicate]) else integer()
result_parts <- if (nrow(existing_results)) list(existing_results) else list()
qc_parts <- if (nrow(existing_qc)) list(existing_qc) else list()
cut_parts <- if (nrow(existing_cuts)) list(existing_cuts) else list()

seed_base <- switch(phase, test = 91000L, preview = 92000L, final = 93000L)
cohort_offset <- c(CHARLS = 100000L, ELSA = 200000L, HRS = 300000L, MHAS = 400000L)[[cohort_name]]
module_offset <- c(
  primary_education = 1000000L,
  primary_wealth = 2000000L,
  wealth_within_low_education = 3000000L
)[[analysis_module]]
resampling_offset <- if (resampling_unit == "stratified_psu") 5000000L else 0L
sample_cluster_multiplicity <- function() {
  if (resampling_unit == "household") {
    sampled <- sample.int(nrow(cluster_frame), nrow(cluster_frame), replace = TRUE)
    return(tabulate(sampled, nbins = nrow(cluster_frame)))
  }
  out <- integer(nrow(cluster_frame))
  strata <- split(cluster_frame$cluster_index, cluster_frame$design_stratum)
  for (ids in strata) {
    sampled <- sample(ids, length(ids), replace = TRUE)
    counts <- table(sampled)
    out[as.integer(names(counts))] <- as.integer(counts)
  }
  out
}
cat(
  "Cluster bootstrap: unit=", resampling_unit, " sensitivity=", sensitivity,
  " module=", analysis_module,
  " cohort=", cohort_name, " clusters=", nrow(cluster_frame),
  " missing_design_people_excluded=", design_missing_people_excluded,
  " range=", replicate_start, "-", replicate_end,
  " completed=", length(completed), "\n", sep = ""
)

for (b in seq.int(replicate_start, replicate_end)) {
  if (b %in% completed) next
  started <- proc.time()[["elapsed"]]
  set.seed(seed_base + module_offset + cohort_offset + resampling_offset + b)
  cluster_mult <- sample_cluster_multiplicity()
  ses_bootstrap <- bootstrap_ses(cluster_mult, b)
  cut_parts[[length(cut_parts) + 1L]] <- ses_bootstrap$cuts

  fboot <- function_base[cluster_mult[cluster_index] > 0]
  fboot[, `:=`(
    boot_mult = cluster_mult[cluster_index], base_weight = respondent_weight,
    ses_boot = ses_bootstrap$map$ses_boot[match(person_id, ses_bootstrap$map$person_id)]
  )]
  mboot <- mortality_base[cluster_mult[cluster_index] > 0]
  mboot[, `:=`(
    boot_mult = cluster_mult[cluster_index], base_weight = origin_weight,
    ses_boot = ses_bootstrap$map$ses_boot[match(person_id, ses_bootstrap$map$person_id)]
  )]
  iboot <- initial_base[cluster_mult[cluster_index] > 0]
  iboot[, `:=`(
    boot_mult = cluster_mult[cluster_index],
    ses_boot = ses_bootstrap$map$ses_boot[match(person_id, ses_bootstrap$map$person_id)]
  )]

  status <- "valid"
  error_message <- ""
  metrics <- NULL
  tryCatch({
    metrics <- calculate_metrics(fboot, mboot, iboot)
    metrics[, `:=`(cohort = cohort_name, phase = phase, replicate = b)]
    setcolorder(metrics, c("cohort", "phase", "replicate", "metric", "estimate"))
  }, error = function(e) {
    status <<- "failed"
    error_message <<- conditionMessage(e)
  })
  if (!is.null(metrics)) result_parts[[length(result_parts) + 1L]] <- metrics
  qc_parts[[length(qc_parts) + 1L]] <- data.table(
    cohort = cohort_name, phase = phase, replicate = b, status = status,
    error_message = error_message,
    resampling_unit = resampling_unit,
    sampled_clusters_unique = sum(cluster_mult > 0),
    sampled_cluster_frequency = sum(cluster_mult),
    design_missing_people_excluded = design_missing_people_excluded,
    sampled_households_unique = sum(cluster_mult > 0),
    sampled_household_frequency = sum(cluster_mult),
    elapsed_seconds = proc.time()[["elapsed"]] - started
  )
  if (b %% 5L == 0L || b == replicate_end) {
    safe_fwrite(rbindlist(result_parts, fill = TRUE), result_path)
    safe_fwrite(rbindlist(qc_parts, fill = TRUE), qc_path)
    safe_fwrite(rbindlist(cut_parts, fill = TRUE), cut_path)
    cat("replicate=", b, " status=", status, " elapsed=", round(proc.time()[["elapsed"]] - started, 2), "s\n", sep = "")
  }
  rm(fboot, mboot, iboot, metrics)
  invisible(gc())
}

writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, paste0("sessionInfo_40_", output_tag, "_", tolower(cohort_name), ".txt"))
)
cat(
  "Completed household bootstrap for ", sensitivity, " / ", analysis_module,
  " / ", cohort_name, "\n", sep = ""
)
