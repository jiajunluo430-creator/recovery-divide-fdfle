#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
base_derived_dir <- file.path(revision_root, "02_derived")
sensitivity <- tolower(Sys.getenv("D4_REV_SENSITIVITY", unset = "primary"))
allowed_sensitivities <- c(
  "primary", "adl_only", "at_least_two", "incident_inception",
  "entry_age_60_64", "pre_covid", "common_interval", "common_calendar",
  "two_year_interval", "piecewise_exponential", "age110_tail", "age100_hazard",
  "wealth_pre_difficulty", "unweighted", "self_report_only", "mortality_common_waves"
)
if (!sensitivity %in% allowed_sensitivities) stop("Unknown D4_REV_SENSITIVITY: ", sensitivity)
derived_dir <- if (sensitivity %in% c("adl_only", "at_least_two")) {
  file.path(base_derived_dir, "sensitivity", sensitivity)
} else {
  base_derived_dir
}
out_dir <- file.path(
  revision_root, "03_outputs",
  if (sensitivity == "primary") "02_revision_point" else paste0("02_revision_point_", sensitivity)
)
log_dir <- file.path(revision_root, "06_logs")
for (d in c(out_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("39_fit_revision_point_models_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

function_risk <- as.data.table(readRDS(file.path(derived_dir, "formal_function_transition_riskset_revision_v1_1.rds")))
mortality_risk <- as.data.table(readRDS(file.path(derived_dir, "formal_mortality_riskset_revision_v1_1.rds")))
initial_data <- as.data.table(readRDS(file.path(derived_dir, "initial_state_data_revision_v1_1.rds")))

cohort_order <- c("CHARLS", "ELSA", "HRS", "MHAS")
ses_levels <- c("high", "middle", "low")
hazard_age_cap <- if (sensitivity == "age100_hazard") 100L else 95L
integration_end_age <- if (sensitivity %in% c("age110_tail", "age100_hazard")) 110L else 100L
ages <- 60:(integration_end_age - 1L)
living_states <- c("I0", "D1", "R1", "D2")
states <- c(living_states, "DEAD")

modules <- list(
  primary_education = list(exposure = "education", low_education_only = FALSE),
  primary_wealth = list(exposure = "wealth", low_education_only = FALSE),
  wealth_within_low_education = list(exposure = "wealth", low_education_only = TRUE)
)
requested_modules <- trimws(strsplit(
  Sys.getenv("D4_REV_MODULES", unset = paste(names(modules), collapse = ",")),
  ",", fixed = TRUE
)[[1L]])
unknown_modules <- setdiff(requested_modules, names(modules))
if (length(unknown_modules)) stop("Unknown D4_REV_MODULES: ", paste(unknown_modules, collapse = ", "))
modules <- modules[requested_modules]
if (!length(modules)) stop("No revision modules selected")

sensitivity_ids <- NULL
entry_audit_path <- file.path(base_derived_dir, "restricted_person_timing_audit_internal.csv")
if (sensitivity %in% c("incident_inception", "entry_age_60_64", "wealth_pre_difficulty")) {
  entry_audit <- fread(entry_audit_path)
  sensitivity_ids <- if (sensitivity == "incident_inception") {
    entry_audit[entry_difficulty == 0, .(cohort, person_id)]
  } else if (sensitivity == "entry_age_60_64") {
    entry_audit[entry_age_60_64 == TRUE, .(cohort, person_id)]
  } else {
    entry_audit[
      !is.na(wealth_entry_wave) & wealth_after_first_observed_difficulty == FALSE,
      .(cohort, person_id)
    ]
  }
}

process_specs <- list(
  onset = list(source = "function", origins = "I0", destination = "D1", event = "event_onset"),
  recovery = list(source = "function", origins = c("D1", "D2"), destination = "R1", event = "event_recovery"),
  relapse = list(source = "function", origins = "R1", destination = "D2", event = "event_relapse"),
  death_pre = list(source = "mortality", origins = "I0", destination = "DEAD", event = "event_death"),
  death_post = list(source = "mortality", origins = c("D1", "R1", "D2"), destination = "DEAD", event = "event_death")
)
mortality_wave_cap <- c(CHARLS = 4L, ELSA = 5L, HRS = 15L, MHAS = 5L)

allowed <- list(
  I0 = c("D1", "DEAD"),
  D1 = c("R1", "DEAD"),
  R1 = c("D2", "DEAD"),
  D2 = c("R1", "DEAD")
)

block_pairs <- list(
  onset = list(c("I0", "D1")),
  recovery = list(c("D1", "R1"), c("D2", "R1")),
  relapse = list(c("R1", "D2")),
  post_difficulty_mortality = list(c("D1", "DEAD"), c("R1", "DEAD"), c("D2", "DEAD")),
  pre_difficulty_mortality = list(c("I0", "DEAD"))
)

weighted_mean_finite <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

make_transition_formula <- function(z, include_period = TRUE, include_proxy = TRUE) {
  age_range <- range(z$age_model, na.rm = TRUE)
  valid_knots <- c(70, 80, 90)
  valid_knots <- valid_knots[valid_knots > age_range[[1L]] & valid_knots < age_range[[2L]]]
  age_term <- if (length(valid_knots)) {
    paste0(
      "ns(age_model, knots = c(", paste(valid_knots, collapse = ","),
      "), Boundary.knots = c(60, ", hazard_age_cap, "))"
    )
  } else {
    "age_model"
  }
  terms <- c("ses", age_term)
  if (nlevels(z$origin_factor) > 1L) terms <- c(terms, "origin_factor")
  if (nlevels(z$sex_factor) > 1L) terms <- c(terms, "sex_factor")
  if (include_period && nlevels(z$period_factor) > 1L) terms <- c(terms, "period_factor")
  if (include_proxy && nlevels(z$proxy_factor) > 1L) terms <- c(terms, "proxy_factor")
  as.formula(paste("event ~", paste(terms, collapse = " + "), "+ offset(log(interval_years))"))
}

prepare_transition_data <- function(cohort_name, module_name, process_name) {
  module <- modules[[module_name]]
  spec <- process_specs[[process_name]]
  if (spec$source == "function") {
    z <- copy(function_risk[cohort == cohort_name & origin_state %in% spec$origins])
    z[, fit_weight_raw := respondent_weight]
  } else {
    z <- copy(mortality_risk[
      cohort == cohort_name & primary_mortality_interval == TRUE & origin_state %in% spec$origins
    ])
    z[, fit_weight_raw := origin_weight]
  }
  if (module$exposure == "wealth") {
    z <- z[wealth_time_eligible == TRUE]
    z[, ses_value := wealth3]
  } else {
    z <- z[education_eligible == TRUE]
    z[, ses_value := education3_fixed]
  }
  if (isTRUE(module$low_education_only)) z <- z[education3_fixed == "low"]
  if (!is.null(sensitivity_ids)) {
    keep_ids <- sensitivity_ids[cohort == cohort_name, person_id]
    z <- z[person_id %in% keep_ids]
  }
  z[, calendar_cutoff_censored := FALSE]
  if (sensitivity == "pre_covid") {
    if (spec$source == "function") {
      z <- z[destination_time < 2020]
    } else {
      # Functional status cannot be imputed at the calendar cutoff, but vital
      # status can be administratively censored. Dropping intervals that cross
      # 2020 while retaining deaths before 2020 creates death-only terminal
      # wave cells and complete separation.
      cutoff <- 2020
      z <- z[origin_time < cutoff]
      z[endpoint_time >= cutoff, `:=`(
        endpoint_time = cutoff,
        interval_years = pmax(cutoff - origin_time, 1 / 365.25),
        death_event = 0L,
        calendar_cutoff_censored = TRUE
      )]
    }
  }
  if (sensitivity == "common_interval") z <- z[interval_years >= 1.5 & interval_years <= 3.25]
  if (sensitivity == "two_year_interval") z <- z[interval_years >= 1.75 & interval_years <= 2.25]
  if (sensitivity == "common_calendar") {
    if (spec$source == "function") {
      z <- z[origin_year >= 2002 & origin_year <= 2018 & destination_time <= 2020]
    } else {
      cutoff <- 2020
      z <- z[origin_year >= 2002 & origin_year <= 2018 & origin_time < cutoff]
      z[endpoint_time >= cutoff, `:=`(
        endpoint_time = cutoff,
        interval_years = pmax(cutoff - origin_time, 1 / 365.25),
        death_event = 0L,
        calendar_cutoff_censored = TRUE
      )]
    }
  }
  if (sensitivity == "mortality_common_waves") {
    z <- z[origin_wave <= mortality_wave_cap[[cohort_name]]]
  }
  if (sensitivity == "self_report_only" && any(!is.na(z$proxy))) z <- z[proxy == 0]
  if (sensitivity == "unweighted") z[, fit_weight_raw := 1]
  z[, `:=`(
    event = as.integer(get(spec$event)),
    ses = factor(ses_value, levels = ses_levels),
    sex_factor = factor(fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown"))),
    proxy_factor = factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self"))),
    period_factor = factor(origin_wave),
    origin_factor = factor(origin_state, levels = spec$origins),
    age_model = pmin(pmax(origin_age, 60), hazard_age_cap)
  )]
  z <- z[
    !is.na(ses) & is.finite(interval_years) & interval_years > 0 &
      is.finite(fit_weight_raw) & fit_weight_raw > 0 & !is.na(age_model)
  ]
  z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
  for (v in c("ses", "sex_factor", "proxy_factor", "period_factor", "origin_factor")) {
    set(z, j = v, value = droplevels(z[[v]]))
  }
  z
}

fit_transition_process <- function(cohort_name, module_name, process_name) {
  spec <- process_specs[[process_name]]
  z <- prepare_transition_data(cohort_name, module_name, process_name)
  if (!all(c("high", "low") %in% levels(z$ses))) stop("Missing high/low SES: ", cohort_name, " / ", module_name, " / ", process_name)
  warnings <- character()
  include_period <- TRUE
  include_proxy <- TRUE
  transition_family <- if (sensitivity == "piecewise_exponential") {
    poisson(link = "log")
  } else {
    binomial(link = "cloglog")
  }
  fit_once <- function() withCallingHandlers(
      glm(
        make_transition_formula(z, include_period, include_proxy),
        data = z, family = transition_family,
        weights = fit_weight, control = glm.control(maxit = 100, epsilon = 1e-9)
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  fit <- fit_once()
  if ((!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) && nlevels(z$proxy_factor) > 1L) {
    include_proxy <- FALSE
    fit <- fit_once()
  }
  if ((!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) && nlevels(z$period_factor) > 1L) {
    include_period <- FALSE
    fit <- fit_once()
  }
  if (!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) stop("Nonconvergent/nonfinite model: ", cohort_name, " / ", module_name, " / ", process_name)

  coefficients <- data.table(
    cohort = cohort_name, module = module_name, process = process_name,
    term = names(coef(fit)), estimate = as.numeric(coef(fit)),
    model_se = sqrt(diag(vcov(fit)))
  )
  coefficients[, `:=`(
    hazard_ratio = exp(estimate),
    model_ci_low = exp(estimate - 1.96 * model_se),
    model_ci_high = exp(estimate + 1.96 * model_se)
  )]

  qc <- data.table(
    cohort = cohort_name, module = module_name, process = process_name,
    intervals = nrow(z), people = uniqueN(z$person_id), events = sum(z$event),
    person_years = sum(z$interval_years), converged = fit$converged,
    iterations = fit$iter, period_adjusted = include_period,
    proxy_adjusted = include_proxy,
    transition_model = if (sensitivity == "piecewise_exponential") "piecewise_exponential_poisson" else "cloglog_interval_probability",
    administrative_censoring_n = sum(z$calendar_cutoff_censored),
    warnings = paste(unique(warnings), collapse = " | ")
  )

  rate_parts <- list()
  for (origin_name in spec$origins) {
    ref <- z[origin_state == origin_name, .(
      fit_weight = sum(fit_weight), standardization_rows = .N
    ), by = .(sex_factor, period_factor, proxy_factor, origin_factor)]
    if (!nrow(ref)) stop("No standardization support for ", origin_name)
    ref[, `:=`(
      sex_factor = factor(sex_factor, levels = levels(z$sex_factor)),
      period_factor = factor(period_factor, levels = levels(z$period_factor)),
      proxy_factor = factor(proxy_factor, levels = levels(z$proxy_factor)),
      origin_factor = factor(origin_factor, levels = levels(z$origin_factor))
    )]
    for (ses_name in ses_levels) {
      if (!ses_name %in% levels(z$ses)) next
      for (age_value in ages) {
        nd <- copy(ref)
        nd[, `:=`(
          ses = factor(ses_name, levels = levels(z$ses)),
          age_model = min(age_value, hazard_age_cap), interval_years = 1
        )]
        hazard <- exp(as.numeric(predict(fit, newdata = nd, type = "link")))
        annual_hazard <- weighted_mean_finite(hazard, nd$fit_weight)
        if (!is.finite(annual_hazard) || annual_hazard < 0) stop("Invalid standardised hazard")
        rate_parts[[length(rate_parts) + 1L]] <- data.table(
          cohort = cohort_name, module = module_name, ses = ses_name,
          age = age_value, process = process_name, origin_state = origin_name,
          destination = spec$destination, annual_hazard,
          standardization_rows = sum(nd$standardization_rows)
        )
      }
    }
  }
  list(coefficients = coefficients, qc = qc, rates = rbindlist(rate_parts))
}

make_initial_formula <- function(z, include_period = TRUE, include_proxy = TRUE) {
  age_span <- diff(range(z$age_model, na.rm = TRUE))
  age_term <- if (is.finite(age_span) && age_span >= 15 && uniqueN(z$age_model) >= 10L) {
    paste0("ns(age_model, knots = c(70, 80, 90), Boundary.knots = c(60, ", hazard_age_cap, "))")
  } else {
    "age_model"
  }
  terms <- c("ses", age_term)
  if (nlevels(z$sex_factor) > 1L) terms <- c(terms, "sex_factor")
  if (include_period && nlevels(z$period_factor) > 1L) terms <- c(terms, "period_factor")
  if (include_proxy && nlevels(z$proxy_factor) > 1L) terms <- c(terms, "proxy_factor")
  as.formula(paste("entry_difficulty ~", paste(terms, collapse = " + ")))
}

fit_initial_state <- function(cohort_name, module_name) {
  module <- modules[[module_name]]
  z <- copy(initial_data[cohort == cohort_name & exposure == module$exposure])
  if (isTRUE(module$low_education_only)) {
    low_ids <- unique(function_risk[cohort == cohort_name & education3_fixed == "low", person_id])
    z <- z[person_id %in% low_ids]
  }
  if (!is.null(sensitivity_ids)) {
    keep_ids <- sensitivity_ids[cohort == cohort_name, person_id]
    z <- z[person_id %in% keep_ids]
  }
  if (sensitivity == "pre_covid") z <- z[entry_time < 2020]
  if (sensitivity == "common_calendar") z <- z[entry_time >= 2002 & entry_time <= 2018]
  if (sensitivity == "mortality_common_waves") z <- z[entry_wave <= mortality_wave_cap[[cohort_name]]]
  if (sensitivity == "self_report_only" && any(!is.na(z$proxy))) z <- z[proxy == 0]
  z[, `:=`(
    ses = factor(level, levels = ses_levels),
    sex_factor = factor(fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown"))),
    proxy_factor = factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self"))),
    period_factor = factor(entry_wave),
    age_model = pmin(pmax(entry_age, 60), hazard_age_cap),
    fit_weight_raw = entry_weight
  )]
  if (sensitivity == "unweighted") z[, fit_weight_raw := 1]
  z <- z[!is.na(ses) & !is.na(entry_difficulty) & is.finite(age_model) & is.finite(fit_weight_raw) & fit_weight_raw > 0]
  z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
  for (v in c("ses", "sex_factor", "proxy_factor", "period_factor")) set(z, j = v, value = droplevels(z[[v]]))
  if (!all(c("high", "low") %in% levels(z$ses))) stop("Initial-state high/low absent: ", cohort_name, " / ", module_name)
  if (sensitivity == "incident_inception") {
    predictions <- data.table(
      cohort = cohort_name, module = module_name,
      ses = ses_levels[ses_levels %in% levels(z$ses)],
      difficulty_probability_age60 = 0,
      standardization_rows = nrow(z)
    )
    qc <- data.table(
      cohort = cohort_name, module = module_name, people = uniqueN(z$person_id),
      difficulty_n = sum(z$entry_difficulty), converged = TRUE,
      iterations = 0L, period_adjusted = FALSE, proxy_adjusted = FALSE,
      degenerate_initial_state = TRUE
    )
    return(list(predictions = predictions, qc = qc))
  }
  observed_initial_values <- unique(z$entry_difficulty[!is.na(z$entry_difficulty)])
  if (length(observed_initial_values) == 1L) {
    constant_probability <- as.numeric(observed_initial_values[[1L]])
    predictions <- data.table(
      cohort = cohort_name, module = module_name,
      ses = ses_levels[ses_levels %in% levels(z$ses)],
      difficulty_probability_age60 = constant_probability,
      standardization_rows = nrow(z)
    )
    qc <- data.table(
      cohort = cohort_name, module = module_name, people = uniqueN(z$person_id),
      difficulty_n = sum(z$entry_difficulty), converged = TRUE,
      iterations = 0L, period_adjusted = FALSE, proxy_adjusted = FALSE,
      degenerate_initial_state = TRUE
    )
    return(list(predictions = predictions, qc = qc))
  }
  include_period <- TRUE
  include_proxy <- TRUE
  fit <- suppressWarnings(glm(
    make_initial_formula(z, include_period, include_proxy),
    data = z, family = binomial(), weights = fit_weight
  ))
  if (any(!is.finite(coef(fit))) && nlevels(z$proxy_factor) > 1L) {
    include_proxy <- FALSE
    fit <- suppressWarnings(glm(
      make_initial_formula(z, include_period, include_proxy),
      data = z, family = binomial(), weights = fit_weight
    ))
  }
  if (any(!is.finite(coef(fit))) && nlevels(z$period_factor) > 1L) {
    include_period <- FALSE
    fit <- suppressWarnings(glm(
      make_initial_formula(z, include_period, include_proxy),
      data = z, family = binomial(), weights = fit_weight
    ))
  }
  if (!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) {
    stop(
      "Initial-state model failure: cohort=", cohort_name,
      " module=", module_name,
      " converged=", fit$converged,
      " coefficients=", paste(names(coef(fit)), format(coef(fit), digits = 6), collapse = " | ")
    )
  }

  ref <- z[, .(fit_weight = sum(fit_weight), standardization_rows = .N), by = .(sex_factor, period_factor, proxy_factor)]
  ref[, `:=`(
    sex_factor = factor(sex_factor, levels = levels(z$sex_factor)),
    period_factor = factor(period_factor, levels = levels(z$period_factor)),
    proxy_factor = factor(proxy_factor, levels = levels(z$proxy_factor))
  )]
  predictions <- rbindlist(lapply(ses_levels[ses_levels %in% levels(z$ses)], function(ses_name) {
    nd <- copy(ref)
    nd[, `:=`(ses = factor(ses_name, levels = levels(z$ses)), age_model = 60)]
    p <- as.numeric(predict(fit, newdata = nd, type = "response"))
    data.table(
      cohort = cohort_name, module = module_name, ses = ses_name,
      difficulty_probability_age60 = weighted_mean_finite(p, nd$fit_weight),
      standardization_rows = sum(nd$standardization_rows)
    )
  }))
  qc <- data.table(
    cohort = cohort_name, module = module_name, people = uniqueN(z$person_id),
    difficulty_n = sum(z$entry_difficulty), converged = fit$converged,
    iterations = fit$iter, period_adjusted = include_period,
    proxy_adjusted = include_proxy, degenerate_initial_state = FALSE
  )
  list(predictions = predictions, qc = qc)
}

lambda_array <- function(z) {
  a <- array(
    NA_real_, c(length(ages), length(living_states), length(states)),
    dimnames = list(age = as.character(ages), origin = living_states, destination = states)
  )
  for (i in seq_len(nrow(z))) {
    a[as.character(z$age[[i]]), z$origin_state[[i]], z$destination[[i]]] <- z$annual_hazard[[i]]
  }
  for (origin_name in names(allowed)) {
    for (destination_name in allowed[[origin_name]]) {
      if (any(!is.finite(a[, origin_name, destination_name]))) stop("Missing lambda ", origin_name, "->", destination_name)
    }
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
  if (max(abs(rowSums(m) - 1)) > 1e-10 || any(m < -1e-12)) stop("Transition matrix failure")
  m
}

life_metrics <- function(lambda, initial_vector, end_age = integration_end_age) {
  active_ages <- ages[ages < end_age]
  v <- initial_vector[states]
  occupancy <- setNames(rep(0, length(states)), states)
  for (age_value in active_ages) {
    m <- transition_matrix(lambda[as.character(age_value), , , drop = TRUE])
    v_next <- as.numeric(v %*% m)
    names(v_next) <- states
    occupancy <- occupancy + (v + v_next) / 2
    v <- v_next
  }
  tle <- sum(occupancy[living_states])
  fdfle <- sum(occupancy[c("I0", "R1")])
  difficulty_years <- sum(occupancy[c("D1", "D2")])
  c(
    total_life_expectancy = tle,
    fdfle = fdfle,
    difficulty_years = difficulty_years,
    fdfle_percent_of_life = if (tle > 0) 100 * fdfle / tle else NA_real_,
    residual_alive_at_end = sum(v[living_states]),
    closure = tle - fdfle - difficulty_years
  )
}

replace_hazard_block <- function(base, donor, block_name) {
  out <- base
  for (pair in block_pairs[[block_name]]) out[, pair[[1]], pair[[2]]] <- donor[, pair[[1]], pair[[2]]]
  out
}

all_subsets <- function(x) {
  unlist(lapply(0:length(x), function(k) {
    if (k == 0L) list(character()) else combn(x, k, simplify = FALSE)
  }), recursive = FALSE)
}

subset_key <- function(x) if (!length(x)) "<empty>" else paste(sort(x), collapse = "+")

shapley_decompose <- function(high_lambda, low_lambda, high_initial, low_initial, include_initial = TRUE) {
  hazard_blocks <- names(block_pairs)
  blocks <- c(if (include_initial) "initial_state", hazard_blocks)
  n_blocks <- length(blocks)
  subsets <- all_subsets(blocks)
  values <- setNames(numeric(length(subsets)), vapply(subsets, subset_key, character(1)))
  for (s in subsets) {
    lambda <- high_lambda
    for (block_name in intersect(s, hazard_blocks)) lambda <- replace_hazard_block(lambda, low_lambda, block_name)
    initial <- if (include_initial && "initial_state" %in% s) low_initial else high_initial
    values[[subset_key(s)]] <- life_metrics(lambda, initial)[["fdfle"]]
  }
  contributions <- setNames(numeric(n_blocks), blocks)
  for (block_name in blocks) {
    others <- setdiff(blocks, block_name)
    for (s in all_subsets(others)) {
      k <- length(s)
      weight <- factorial(k) * factorial(n_blocks - k - 1L) / factorial(n_blocks)
      contributions[[block_name]] <- contributions[[block_name]] + weight *
        (values[[subset_key(c(s, block_name))]] - values[[subset_key(s)]])
    }
  }
  low_value <- values[[subset_key(blocks)]]
  high_value <- values[[subset_key(character())]]
  list(
    low = low_value, high = high_value, gap = low_value - high_value,
    contributions = contributions,
    closure = sum(contributions) - (low_value - high_value)
  )
}

coefficient_parts <- list()
model_qc_parts <- list()
rate_parts <- list()
initial_parts <- list()
initial_qc_parts <- list()

for (module_name in names(modules)) {
  for (cohort_name in cohort_order) {
    cat("Initial-state model: ", module_name, " / ", cohort_name, "\n", sep = "")
    init <- fit_initial_state(cohort_name, module_name)
    initial_parts[[length(initial_parts) + 1L]] <- init$predictions
    initial_qc_parts[[length(initial_qc_parts) + 1L]] <- init$qc
    for (process_name in names(process_specs)) {
      cat("Transition model: ", module_name, " / ", cohort_name, " / ", process_name, "\n", sep = "")
      ans <- fit_transition_process(cohort_name, module_name, process_name)
      coefficient_parts[[length(coefficient_parts) + 1L]] <- ans$coefficients
      model_qc_parts[[length(model_qc_parts) + 1L]] <- ans$qc
      rate_parts[[length(rate_parts) + 1L]] <- ans$rates
      rm(ans)
      invisible(gc())
    }
  }
}

coefficients <- rbindlist(coefficient_parts, fill = TRUE)
model_qc <- rbindlist(model_qc_parts, fill = TRUE)
rates <- rbindlist(rate_parts, fill = TRUE)
initial_predictions <- rbindlist(initial_parts, fill = TRUE)
initial_qc <- rbindlist(initial_qc_parts, fill = TRUE)

life_parts <- list()
gap_parts <- list()
shapley_parts <- list()
for (module_name in names(modules)) {
  for (cohort_name in cohort_order) {
    lambda_store <- list()
    init_store <- list()
    for (ses_name in c("high", "low")) {
      lambda_store[[ses_name]] <- lambda_array(rates[module == module_name & cohort == cohort_name & ses == ses_name])
      p_diff <- initial_predictions[module == module_name & cohort == cohort_name & ses == ses_name, difficulty_probability_age60]
      if (length(p_diff) != 1L || !is.finite(p_diff) || p_diff < 0 || p_diff > 1) stop("Invalid initial probability")
      init_store[[ses_name]] <- setNames(c(1 - p_diff, p_diff, 0, 0, 0), states)
      conditional_init <- setNames(c(1, 0, 0, 0, 0), states)
      for (estimand_name in c("population_initialised", "conditional_all_I0")) {
        initial <- if (estimand_name == "population_initialised") init_store[[ses_name]] else conditional_init
        metrics <- life_metrics(lambda_store[[ses_name]], initial)
        life_parts[[length(life_parts) + 1L]] <- data.table(
          cohort = cohort_name, module = module_name, ses = ses_name,
          estimand = estimand_name, metric = names(metrics), estimate = as.numeric(metrics)
        )
      }
    }

    pop_dec <- shapley_decompose(
      lambda_store[["high"]], lambda_store[["low"]],
      init_store[["high"]], init_store[["low"]], include_initial = TRUE
    )
    conditional_init <- setNames(c(1, 0, 0, 0, 0), states)
    con_dec <- shapley_decompose(
      lambda_store[["high"]], lambda_store[["low"]],
      conditional_init, conditional_init, include_initial = FALSE
    )
    decomposition_store <- list(population_initialised = pop_dec, conditional_all_I0 = con_dec)
    for (estimand_name in names(decomposition_store)) {
      dec <- decomposition_store[[estimand_name]]
      if (!is.finite(dec$closure) || abs(dec$closure) > 0.01) stop("Shapley closure failure")
      gap_parts[[length(gap_parts) + 1L]] <- data.table(
        cohort = cohort_name, module = module_name, estimand = estimand_name,
        high_fdfle = dec$high, low_fdfle = dec$low, low_minus_high_gap = dec$gap,
        shapley_closure_error = dec$closure
      )
      shapley_parts[[length(shapley_parts) + 1L]] <- data.table(
        cohort = cohort_name, module = module_name, estimand = estimand_name,
        block = names(dec$contributions), contribution_years = as.numeric(dec$contributions)
      )
    }
  }
}

life <- rbindlist(life_parts)
gaps <- rbindlist(gap_parts)
shapley <- rbindlist(shapley_parts)

fwrite(coefficients, file.path(out_dir, "transition_model_coefficients.csv"))
fwrite(model_qc, file.path(out_dir, "transition_model_qc.csv"))
fwrite(rates, file.path(out_dir, "standardised_annual_hazards.csv"))
fwrite(initial_predictions, file.path(out_dir, "initial_difficulty_probability_age60.csv"))
fwrite(initial_qc, file.path(out_dir, "initial_state_model_qc.csv"))
fwrite(life, file.path(out_dir, "life_expectancy_point_estimates.csv"))
fwrite(gaps, file.path(out_dir, "low_high_fdfle_gaps.csv"))
fwrite(shapley, file.path(out_dir, "shapley_decomposition_point.csv"))
fwrite(data.table(
  sensitivity = sensitivity,
  start_age = 60L,
  integration_end_age = integration_end_age,
  hazard_age_cap = hazard_age_cap,
  tail_rule = if (integration_end_age > hazard_age_cap) {
    paste0("annual hazards frozen at age ", hazard_age_cap, " through age ", integration_end_age)
  } else {
    "no frozen tail"
  }
), file.path(out_dir, "analysis_horizon_qc.csv"))

cat("\nPopulation-initialised gaps:\n")
print(gaps[estimand == "population_initialised"][order(module, cohort)])
cat("\nPopulation Shapley contributions:\n")
print(shapley[estimand == "population_initialised"][order(module, cohort, block)])
writeLines(capture.output(sessionInfo()), file.path(log_dir, paste0("sessionInfo_39_revision_point_models_", sensitivity, ".txt")))
