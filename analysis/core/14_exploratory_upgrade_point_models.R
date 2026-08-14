#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
derived_dir <- file.path(root, "02_derived")
sensitivity <- tolower(Sys.getenv("D4_V2_SENSITIVITY", unset = "primary"))
allowed_sensitivities <- c("primary", "unweighted", "interval_1_to_3_years", "exclude_explicit_proxy")
if (!sensitivity %in% allowed_sensitivities) {
  stop("Unknown D4_V2_SENSITIVITY: ", sensitivity)
}
out_dir <- file.path(
  root, "03_outputs",
  if (sensitivity == "primary") "11_exploratory_upgrade_point" else paste0("15_exploratory_sensitivity_", sensitivity)
)
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("14_exploratory_upgrade_point_models_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

function_risk <- as.data.table(readRDS(file.path(derived_dir, "formal_function_transition_riskset.rds")))
mortality_risk <- as.data.table(readRDS(file.path(derived_dir, "formal_mortality_riskset.rds")))
panel <- as.data.table(readRDS(file.path(derived_dir, "individual_wave_panel_gate0.rds")))

function_risk[, event_recovery := as.integer(origin_state %in% c("D1", "D2") & destination == "R1")]

cohort_order <- c("CHARLS", "HRS", "ELSA", "MHAS")
ses_levels <- c("high", "middle", "low")
ages <- 60:99
age_profile_values <- c(60L, 70L, 80L, 90L)

process_specs <- list(
  onset = list(source = "function", origins = "I0", destination = "D1", event = "event_onset"),
  recovery = list(source = "function", origins = c("D1", "D2"), destination = "R1", event = "event_recovery"),
  relapse = list(source = "function", origins = "R1", destination = "D2", event = "event_relapse"),
  death_pre = list(source = "mortality", origins = "I0", destination = "DEAD", event = "event_death"),
  death_post = list(source = "mortality", origins = c("D1", "R1", "D2"), destination = "DEAD", event = "event_death")
)

modules <- list(
  wealth_by_sex = list(mode = "sex_interaction", low_education_only = FALSE, strata = c("male", "female")),
  wealth_within_low_education = list(mode = "main", low_education_only = TRUE, strata = "all"),
  wealth_age_varying = list(mode = "age_interaction", low_education_only = FALSE, strata = "all")
)
requested_modules <- trimws(strsplit(
  Sys.getenv("D4_V2_MODULES", unset = paste(names(modules), collapse = ",")),
  ",", fixed = TRUE
)[[1L]])
unknown_modules <- setdiff(requested_modules, names(modules))
if (length(unknown_modules)) stop("Unknown D4_V2_MODULES: ", paste(unknown_modules, collapse = ", "))
modules <- modules[requested_modules]
if (!length(modules)) stop("No exploratory modules selected")
cat("Selected life-table modules: ", paste(names(modules), collapse = ", "), "\n", sep = "")
cat("Sensitivity implementation: ", sensitivity, "\n", sep = "")

prepare_model_data <- function(cohort_name, process_name, module_name) {
  spec <- process_specs[[process_name]]
  module <- modules[[module_name]]
  if (spec$source == "function") {
    z <- copy(function_risk[cohort == cohort_name & origin_state %in% spec$origins])
    z[, fit_weight_raw := respondent_weight]
  } else {
    z <- copy(mortality_risk[
      cohort == cohort_name & primary_mortality_interval == TRUE & origin_state %in% spec$origins
    ])
    z[, fit_weight_raw := origin_weight]
  }
  if (isTRUE(module$low_education_only)) z <- z[education3_fixed == "low"]
  if (module$mode == "sex_interaction") z <- z[female %in% c(0, 1)]
  if (sensitivity == "interval_1_to_3_years") z <- z[interval_years >= 1 & interval_years <= 3]
  if (sensitivity == "exclude_explicit_proxy") z <- z[is.na(proxy) | proxy != 1]
  if (sensitivity == "unweighted") z[, fit_weight_raw := 1]
  z[, event := as.integer(get(spec$event))]
  z[, ses := factor(wealth3, levels = ses_levels)]
  z[, sex_factor := factor(
    fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown")),
    levels = c("male", "female", "unknown")
  )]
  z[, proxy_factor := factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self")))]
  z[, period_factor := factor(origin_wave)]
  z[, origin_factor := factor(origin_state, levels = spec$origins)]
  z[, age_model := pmin(pmax(origin_age, 60), 95)]
  z[, age_interaction := (pmin(pmax(origin_age, 60), 90) - 70) / 10]
  z <- z[
    !is.na(ses) & is.finite(interval_years) & interval_years > 0 &
      is.finite(fit_weight_raw) & fit_weight_raw > 0 & !is.na(age_model)
  ]
  z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
  for (column_name in c("ses", "sex_factor", "proxy_factor", "period_factor", "origin_factor")) {
    set(z, j = column_name, value = droplevels(z[[column_name]]))
  }
  z
}

make_formula <- function(z, mode) {
  age_main <- "ns(age_model, knots = c(70, 80, 90), Boundary.knots = c(60, 95))"
  if (mode == "sex_interaction") {
    terms <- c("ses * sex_factor", age_main)
  } else if (mode == "age_interaction") {
    terms <- c("ses", age_main, "ses:age_interaction")
  } else {
    terms <- c("ses", age_main)
  }
  if (nlevels(z$origin_factor) > 1L) terms <- c(terms, "origin_factor")
  if (mode != "sex_interaction" && nlevels(z$sex_factor) > 1L) terms <- c(terms, "sex_factor")
  if (nlevels(z$period_factor) > 1L) terms <- c(terms, "period_factor")
  if (nlevels(z$proxy_factor) > 1L) terms <- c(terms, "proxy_factor")
  as.formula(paste("event ~", paste(terms, collapse = " + "), "+ offset(log(interval_years))"))
}

weighted_mean_finite <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

fit_one_process <- function(cohort_name, process_name, module_name) {
  spec <- process_specs[[process_name]]
  module <- modules[[module_name]]
  z <- prepare_model_data(cohort_name, process_name, module_name)
  if (!all(c("high", "low") %in% levels(z$ses))) {
    stop(module_name, " ", cohort_name, " ", process_name, ": low/high wealth absent")
  }
  if (module$mode == "sex_interaction" && !all(c("male", "female") %in% levels(z$sex_factor))) {
    stop(module_name, " ", cohort_name, " ", process_name, ": male/female level absent")
  }
  formula <- make_formula(z, module$mode)
  warnings <- character()
  fit <- withCallingHandlers(
    glm(
      formula, data = z, family = binomial(link = "cloglog"),
      weights = fit_weight, control = glm.control(maxit = 100, epsilon = 1e-9)
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) {
    stop(module_name, " ", cohort_name, " ", process_name, ": nonconvergent/nonfinite model")
  }

  coef_mat <- summary(fit)$coefficients
  coef_dt <- data.table(
    module = module_name,
    cohort = cohort_name,
    process = process_name,
    term = rownames(coef_mat),
    estimate = coef_mat[, 1],
    std_error_model = coef_mat[, 2],
    statistic = coef_mat[, 3],
    p_value_model = coef_mat[, 4]
  )
  coef_dt[, `:=`(
    hazard_ratio = exp(estimate),
    hr_ci_low_model = exp(estimate - 1.96 * std_error_model),
    hr_ci_high_model = exp(estimate + 1.96 * std_error_model)
  )]

  qc <- z[, .(
    rows = .N,
    persons = uniqueN(person_id),
    events = sum(event),
    low_rows = sum(ses == "low"),
    low_events = sum(event[ses == "low"]),
    high_rows = sum(ses == "high"),
    high_events = sum(event[ses == "high"])
  )]
  qc[, `:=`(
    module = module_name,
    cohort = cohort_name,
    process = process_name,
    converged = fit$converged,
    iterations = fit$iter,
    coefficient_abs_max = max(abs(coef(fit))),
    warning_n = length(unique(warnings)),
    warnings = paste(unique(warnings), collapse = " | "),
    formula = paste(deparse(formula), collapse = " ")
  )]

  rate_parts <- list()
  for (origin_name in spec$origins) {
    if (module$mode == "sex_interaction") {
      ref <- z[origin_state == origin_name, .(
        fit_weight = sum(fit_weight), standardization_rows = .N
      ), by = .(period_factor, proxy_factor, origin_factor)]
    } else {
      ref <- z[origin_state == origin_name, .(
        fit_weight = sum(fit_weight), standardization_rows = .N
      ), by = .(sex_factor, period_factor, proxy_factor, origin_factor)]
    }
    if (!nrow(ref)) stop("No standardization rows for ", module_name, " ", cohort_name, " ", process_name)
    ref[, `:=`(
      period_factor = factor(period_factor, levels = levels(z$period_factor)),
      proxy_factor = factor(proxy_factor, levels = levels(z$proxy_factor)),
      origin_factor = factor(origin_factor, levels = levels(z$origin_factor))
    )]
    if (module$mode != "sex_interaction") {
      ref[, sex_factor := factor(sex_factor, levels = levels(z$sex_factor))]
    }
    for (stratum_name in module$strata) {
      for (ses_name in ses_levels) {
        if (!ses_name %in% levels(z$ses)) next
        for (age_value in ages) {
          eval_age <- min(age_value, 95)
          nd <- copy(ref)
          nd[, `:=`(
            ses = factor(ses_name, levels = levels(z$ses)),
            age_model = eval_age,
            age_interaction = (min(max(age_value, 60), 90) - 70) / 10,
            interval_years = 1
          )]
          if (module$mode == "sex_interaction") {
            nd[, sex_factor := factor(stratum_name, levels = levels(z$sex_factor))]
          }
          link <- as.numeric(predict(fit, newdata = nd, type = "link"))
          hazard <- exp(link)
          annual_hazard <- weighted_mean_finite(hazard, nd$fit_weight)
          if (!is.finite(annual_hazard) || annual_hazard < 0) {
            stop("Invalid standardized hazard for ", module_name, " ", cohort_name, " ", process_name)
          }
          rate_parts[[length(rate_parts) + 1L]] <- data.table(
            module = module_name,
            stratum = stratum_name,
            cohort = cohort_name,
            ses = ses_name,
            age = age_value,
            process = process_name,
            origin_state = origin_name,
            destination = spec$destination,
            annual_hazard = annual_hazard,
            standardization_rows = sum(nd$standardization_rows)
          )
        }
      }
    }
  }

  list(coefficients = coef_dt, qc = qc, rates = rbindlist(rate_parts))
}

coefficient_parts <- list()
qc_parts <- list()
rate_parts <- list()
for (module_name in names(modules)) {
  for (cohort_name in cohort_order) {
    for (process_name in names(process_specs)) {
      cat("Fitting ", module_name, " / ", cohort_name, " / ", process_name, "...\n", sep = "")
      ans <- fit_one_process(cohort_name, process_name, module_name)
      coefficient_parts[[length(coefficient_parts) + 1L]] <- ans$coefficients
      qc_parts[[length(qc_parts) + 1L]] <- ans$qc
      rate_parts[[length(rate_parts) + 1L]] <- ans$rates
      rm(ans)
      invisible(gc())
    }
  }
}

coefficients <- rbindlist(coefficient_parts, fill = TRUE)
model_qc <- rbindlist(qc_parts, fill = TRUE)
rates <- rbindlist(rate_parts, fill = TRUE)

living_states <- c("I0", "D1", "R1", "D2")
states <- c(living_states, "DEAD")
allowed <- list(
  I0 = c("D1", "DEAD"),
  D1 = c("R1", "DEAD"),
  R1 = c("D2", "DEAD"),
  D2 = c("R1", "DEAD")
)

lambda_array <- function(z) {
  a <- array(
    NA_real_,
    dim = c(length(ages), length(living_states), length(states)),
    dimnames = list(age = as.character(ages), origin = living_states, destination = states)
  )
  for (i in seq_len(nrow(z))) {
    a[as.character(z$age[[i]]), z$origin_state[[i]], z$destination[[i]]] <- z$annual_hazard[[i]]
  }
  for (origin_name in names(allowed)) {
    for (destination_name in allowed[[origin_name]]) {
      if (any(!is.finite(a[, origin_name, destination_name]))) {
        stop("Missing hazard for ", origin_name, " -> ", destination_name)
      }
    }
  }
  a
}

transition_matrix <- function(lambda_age) {
  m <- matrix(0, nrow = length(states), ncol = length(states), dimnames = list(states, states))
  m["DEAD", "DEAD"] <- 1
  for (origin_name in names(allowed)) {
    destinations <- allowed[[origin_name]]
    lambda <- lambda_age[origin_name, destinations]
    total <- sum(lambda)
    move <- if (total > 0) 1 - exp(-total) else 0
    m[origin_name, origin_name] <- 1 - move
    if (total > 0) m[origin_name, destinations] <- move * lambda / total
  }
  if (max(abs(rowSums(m) - 1)) > 1e-12 || any(m < -1e-12)) stop("Transition matrix failure")
  m
}

life_metrics <- function(lambda) {
  v <- setNames(c(1, 0, 0, 0, 0), states)
  occupancy <- setNames(rep(0, length(states)), states)
  for (age_value in ages) {
    m <- transition_matrix(lambda[as.character(age_value), , , drop = TRUE])
    v_next <- as.numeric(v %*% m)
    names(v_next) <- states
    occupancy <- occupancy + (v + v_next) / 2
    v <- v_next
  }
  c(
    total_life_expectancy = sum(occupancy[living_states]),
    dfle = sum(occupancy[c("I0", "R1")]),
    disabled_years = sum(occupancy[c("D1", "D2")]),
    residual_alive_age100 = sum(v[living_states])
  )
}

block_pairs <- list(
  onset = list(c("I0", "D1")),
  recovery = list(c("D1", "R1"), c("D2", "R1")),
  relapse = list(c("R1", "D2")),
  post_disability_mortality = list(c("D1", "DEAD"), c("R1", "DEAD"), c("D2", "DEAD")),
  pre_disability_mortality = list(c("I0", "DEAD"))
)
blocks <- names(block_pairs)

permutations <- function(x) {
  if (length(x) == 1L) return(list(x))
  out <- list()
  for (i in seq_along(x)) {
    for (tail_value in permutations(x[-i])) out[[length(out) + 1L]] <- c(x[[i]], tail_value)
  }
  out
}
block_permutations <- permutations(blocks)

replace_block <- function(base, donor, block_name) {
  out <- base
  for (pair in block_pairs[[block_name]]) out[, pair[[1]], pair[[2]]] <- donor[, pair[[1]], pair[[2]]]
  out
}

array_store <- list()
life_parts <- list()
for (module_name in unique(rates$module)) {
  for (stratum_name in unique(rates[module == module_name, stratum])) {
    for (cohort_name in cohort_order) {
      for (ses_name in ses_levels) {
        z <- rates[module == module_name & stratum == stratum_name & cohort == cohort_name & ses == ses_name]
        if (!nrow(z)) next
        key <- paste(module_name, stratum_name, cohort_name, ses_name, sep = "__")
        lambda <- lambda_array(z)
        array_store[[key]] <- lambda
        metrics <- life_metrics(lambda)
        life_parts[[length(life_parts) + 1L]] <- data.table(
          module = module_name,
          stratum = stratum_name,
          cohort = cohort_name,
          ses = ses_name,
          total_life_expectancy = unname(metrics[["total_life_expectancy"]]),
          dfle = unname(metrics[["dfle"]]),
          disabled_years = unname(metrics[["disabled_years"]]),
          residual_alive_age100 = unname(metrics[["residual_alive_age100"]]),
          state_years_closure = unname(metrics[["total_life_expectancy"]] - metrics[["dfle"]] - metrics[["disabled_years"]])
        )
      }
    }
  }
}
life_expectancy <- rbindlist(life_parts)

metric_names <- c("dfle", "disabled_years", "total_life_expectancy")
decomposition_parts <- list()
gap_parts <- list()
for (module_name in unique(rates$module)) {
  for (stratum_name in unique(rates[module == module_name, stratum])) {
    for (cohort_name in cohort_order) {
      high_lambda <- array_store[[paste(module_name, stratum_name, cohort_name, "high", sep = "__")]]
      low_lambda <- array_store[[paste(module_name, stratum_name, cohort_name, "low", sep = "__")]]
      high_metrics <- life_metrics(high_lambda)
      low_metrics <- life_metrics(low_lambda)
      contributions <- matrix(0, nrow = length(blocks), ncol = length(metric_names), dimnames = list(blocks, metric_names))
      for (perm in block_permutations) {
        hybrid <- high_lambda
        previous <- high_metrics
        for (block_name in perm) {
          hybrid <- replace_block(hybrid, low_lambda, block_name)
          current <- life_metrics(hybrid)
          contributions[block_name, ] <- contributions[block_name, ] + current[metric_names] - previous[metric_names]
          previous <- current
        }
      }
      contributions <- contributions / length(block_permutations)
      gaps <- low_metrics[metric_names] - high_metrics[metric_names]
      for (metric_name in metric_names) {
        for (block_name in blocks) {
          value <- contributions[block_name, metric_name]
          decomposition_parts[[length(decomposition_parts) + 1L]] <- data.table(
            module = module_name,
            stratum = stratum_name,
            cohort = cohort_name,
            estimand = metric_name,
            block = block_name,
            contribution_years = value,
            total_low_minus_high_gap = unname(gaps[[metric_name]]),
            contribution_percent = if (abs(gaps[[metric_name]]) > 1e-10) 100 * value / gaps[[metric_name]] else NA_real_
          )
        }
        gap_parts[[length(gap_parts) + 1L]] <- data.table(
          module = module_name,
          stratum = stratum_name,
          cohort = cohort_name,
          estimand = metric_name,
          high = unname(high_metrics[[metric_name]]),
          low = unname(low_metrics[[metric_name]]),
          low_minus_high_gap = unname(gaps[[metric_name]]),
          shapley_sum = sum(contributions[, metric_name]),
          closure_error = sum(contributions[, metric_name]) - unname(gaps[[metric_name]])
        )
      }
    }
  }
}
decomposition <- rbindlist(decomposition_parts)
absolute_gaps <- rbindlist(gap_parts)

point_screen <- decomposition[
  estimand == "dfle" & block %in% c("recovery", "relapse"),
  .(
    recovery_relapse_contribution_years = sum(contribution_years),
    dfle_gap = unique(total_low_minus_high_gap)
  ), by = .(module, stratum, cohort)
]
point_screen[, recovery_relapse_percent := fifelse(
  abs(dfle_gap) > 1e-10,
  100 * recovery_relapse_contribution_years / dfle_gap,
  NA_real_
)]
point_screen[, point_threshold_met :=
  abs(recovery_relapse_contribution_years) >= 0.5 & abs(recovery_relapse_percent) >= 20]

if ("wealth_by_sex" %in% point_screen$module) {
  sex_wide <- dcast(
    point_screen[module == "wealth_by_sex"],
    cohort ~ stratum,
    value.var = c("dfle_gap", "recovery_relapse_contribution_years", "recovery_relapse_percent")
  )
  sex_wide[, `:=`(
    female_minus_male_recovery_relapse = recovery_relapse_contribution_years_female - recovery_relapse_contribution_years_male,
    female_minus_male_dfle_gap = dfle_gap_female - dfle_gap_male
  )]
  sex_wide[, point_sex_pilot_positive := abs(female_minus_male_recovery_relapse) >= 0.5]
} else {
  sex_wide <- data.table()
}

lowedu_wealth_screen <- copy(point_screen[module == "wealth_within_low_education"])
if (nrow(lowedu_wealth_screen)) {
  lowedu_wealth_screen[, point_compensation_pilot_positive :=
    abs(recovery_relapse_contribution_years) >= 0.5 & abs(recovery_relapse_percent) >= 20]
}

if ("wealth_age_varying" %in% rates$module) {
  age_rate_profile <- rates[
    module == "wealth_age_varying" & age %in% age_profile_values & ses %in% c("high", "low") &
      process %in% c("recovery", "relapse")
  ]
  age_rate_profile <- dcast(
    age_rate_profile,
    cohort + process + origin_state + age ~ ses,
    value.var = "annual_hazard"
  )
  age_rate_profile[, low_high_hazard_ratio := low / high]
  age_change <- dcast(
    age_rate_profile[age %in% c(60, 80)],
    cohort + process + origin_state ~ age,
    value.var = "low_high_hazard_ratio"
  )
  setnames(age_change, c("60", "80"), c("ratio_age60", "ratio_age80"))
  age_change[, `:=`(
    ratio80_div_ratio60 = ratio_age80 / ratio_age60,
    point_age_pilot_positive = abs(ratio_age80 / ratio_age60 - 1) >= 0.20
  )]
} else {
  age_rate_profile <- data.table()
  age_change <- data.table()
}

# Recovery phase model: wealth-specific relapse after first versus sustained recovery.
setorder(panel, cohort, person_id, wave)
panel[, `:=`(
  previous_wave = shift(wave),
  previous_observed = shift(observed_state),
  previous_history_state = shift(history_state)
), by = .(cohort, person_id)]
panel[, recovery_phase := fcase(
  observed_state == TRUE & history_state == "R1" & previous_wave == wave - 1L &
    previous_observed == TRUE & previous_history_state %in% c("D1", "D2"), "early_recovery",
  observed_state == TRUE & history_state == "R1" & previous_wave == wave - 1L &
    previous_observed == TRUE & previous_history_state == "R1", "sustained_recovery",
  observed_state == TRUE & history_state == "R1", "unclassified",
  default = NA_character_
)]
phase_map <- unique(panel[
  observed_state == TRUE & history_state == "R1",
  .(cohort, person_id, origin_wave = wave, recovery_phase)
])
phase_risk <- merge(
  function_risk, phase_map,
  by = c("cohort", "person_id", "origin_wave"), all.x = TRUE, sort = FALSE
)

linear_contrast <- function(fit, weights) {
  beta <- coef(fit)
  vv <- vcov(fit)
  w <- setNames(rep(0, length(beta)), names(beta))
  for (name in names(weights)) {
    if (!name %in% names(w)) stop("Missing coefficient term: ", name)
    w[[name]] <- weights[[name]]
  }
  estimate <- sum(w * beta)
  se <- sqrt(as.numeric(t(w) %*% vv %*% w))
  data.table(
    log_hazard_ratio = estimate,
    std_error_model = se,
    hazard_ratio = exp(estimate),
    hr_ci_low_model = exp(estimate - 1.96 * se),
    hr_ci_high_model = exp(estimate + 1.96 * se),
    p_value_model = 2 * pnorm(-abs(estimate / se))
  )
}

phase_contrast_parts <- list()
phase_rate_parts <- list()
phase_qc_parts <- list()
for (cohort_name in cohort_order) {
  z <- copy(phase_risk[
    cohort == cohort_name & origin_state == "R1" &
      recovery_phase %in% c("early_recovery", "sustained_recovery") &
      wealth3 %in% ses_levels
  ])
  z[, `:=`(
    event = as.integer(event_relapse),
    ses = factor(wealth3, levels = ses_levels),
    recovery_phase = factor(recovery_phase, levels = c("early_recovery", "sustained_recovery")),
    sex_factor = factor(fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown"))),
    proxy_factor = factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self"))),
    period_factor = factor(origin_wave),
    age_model = pmin(pmax(origin_age, 60), 95),
    fit_weight_raw = respondent_weight
  )]
  if (sensitivity == "interval_1_to_3_years") z <- z[interval_years >= 1 & interval_years <= 3]
  if (sensitivity == "exclude_explicit_proxy") z <- z[is.na(proxy) | proxy != 1]
  if (sensitivity == "unweighted") z[, fit_weight_raw := 1]
  z <- z[
    is.finite(interval_years) & interval_years > 0 & is.finite(fit_weight_raw) &
      fit_weight_raw > 0 & !is.na(age_model)
  ]
  z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
  for (column_name in c("ses", "recovery_phase", "sex_factor", "proxy_factor", "period_factor")) {
    set(z, j = column_name, value = droplevels(z[[column_name]]))
  }
  phase_terms <- c(
    "ses * recovery_phase",
    "ns(age_model, knots = c(70, 80, 90), Boundary.knots = c(60, 95))"
  )
  if (nlevels(z$sex_factor) > 1L) phase_terms <- c(phase_terms, "sex_factor")
  if (nlevels(z$period_factor) > 1L) phase_terms <- c(phase_terms, "period_factor")
  if (nlevels(z$proxy_factor) > 1L) phase_terms <- c(phase_terms, "proxy_factor")
  phase_formula <- as.formula(paste(
    "event ~", paste(phase_terms, collapse = " + "), "+ offset(log(interval_years))"
  ))
  fit <- glm(
    phase_formula, data = z, family = binomial(link = "cloglog"),
    weights = fit_weight, control = glm.control(maxit = 100, epsilon = 1e-9)
  )
  if (!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) stop("Recovery-phase model failed: ", cohort_name)
  interaction_term <- "seslow:recovery_phasesustained_recovery"
  early <- linear_contrast(fit, c(seslow = 1))
  early[, `:=`(cohort = cohort_name, contrast = "low_vs_high_early_recovery")]
  sustained <- linear_contrast(fit, setNames(c(1, 1), c("seslow", interaction_term)))
  sustained[, `:=`(cohort = cohort_name, contrast = "low_vs_high_sustained_recovery")]
  modification <- linear_contrast(fit, setNames(1, interaction_term))
  modification[, `:=`(cohort = cohort_name, contrast = "sustained_vs_early_modification")]
  phase_contrast_parts[[length(phase_contrast_parts) + 1L]] <- rbindlist(list(early, sustained, modification), fill = TRUE)

  phase_qc_one <- z[, .(
    intervals = .N,
    persons = uniqueN(person_id),
    relapses = sum(event)
  ), by = .(recovery_phase, ses)]
  phase_qc_one[, cohort := cohort_name]
  setcolorder(phase_qc_one, c("cohort", "recovery_phase", "ses", "intervals", "persons", "relapses"))
  phase_qc_parts[[length(phase_qc_parts) + 1L]] <- phase_qc_one

  ref <- z[, .(fit_weight = sum(fit_weight)), by = .(sex_factor, period_factor, proxy_factor)]
  ref[, `:=`(
    sex_factor = factor(sex_factor, levels = levels(z$sex_factor)),
    period_factor = factor(period_factor, levels = levels(z$period_factor)),
    proxy_factor = factor(proxy_factor, levels = levels(z$proxy_factor))
  )]
  for (phase_name in c("early_recovery", "sustained_recovery")) {
    for (ses_name in c("high", "low")) {
      for (age_value in age_profile_values) {
        nd <- copy(ref)
        nd[, `:=`(
          ses = factor(ses_name, levels = levels(z$ses)),
          recovery_phase = factor(phase_name, levels = levels(z$recovery_phase)),
          age_model = age_value,
          interval_years = 1
        )]
        hazard <- exp(as.numeric(predict(fit, newdata = nd, type = "link")))
        phase_rate_parts[[length(phase_rate_parts) + 1L]] <- data.table(
          cohort = cohort_name,
          recovery_phase = phase_name,
          ses = ses_name,
          age = age_value,
          annual_relapse_hazard = weighted_mean_finite(hazard, nd$fit_weight)
        )
      }
    }
  }
}

phase_contrasts <- rbindlist(phase_contrast_parts, fill = TRUE)
phase_contrasts[, point_phase_pilot_positive := fifelse(
  contrast %in% c("low_vs_high_early_recovery", "low_vs_high_sustained_recovery"),
  hazard_ratio >= 1.15,
  abs(hazard_ratio - 1) >= 0.20
)]
phase_rates <- rbindlist(phase_rate_parts)
phase_rate_ratios <- dcast(
  phase_rates,
  cohort + recovery_phase + age ~ ses,
  value.var = "annual_relapse_hazard"
)
phase_rate_ratios[, low_high_hazard_ratio := low / high]
phase_qc <- rbindlist(phase_qc_parts, fill = TRUE)

if (max(abs(absolute_gaps$closure_error)) > 0.01) stop("Exploratory Shapley closure exceeds 0.01 years")
if (max(abs(life_expectancy$state_years_closure)) > 1e-8) stop("Exploratory life-table closure failure")

fwrite(model_qc, file.path(out_dir, "model_convergence_qc.csv"))
fwrite(coefficients, file.path(out_dir, "model_coefficients_model_based.csv"))
fwrite(rates, file.path(out_dir, "standardized_annual_hazards.csv"))
fwrite(life_expectancy, file.path(out_dir, "life_expectancy_point_estimates.csv"))
fwrite(absolute_gaps, file.path(out_dir, "low_high_absolute_gaps.csv"))
fwrite(decomposition, file.path(out_dir, "shapley_decomposition_point.csv"))
fwrite(point_screen, file.path(out_dir, "module_point_promotion_screen.csv"))
if (nrow(sex_wide)) {
  fwrite(sex_wide, file.path(out_dir, "sex_difference_point_screen.csv"))
}
fwrite(lowedu_wealth_screen, file.path(out_dir, "low_education_wealth_compensation_point_screen.csv"))
if (nrow(age_rate_profile)) {
  fwrite(age_rate_profile, file.path(out_dir, "age_varying_wealth_rate_profile.csv"))
}
if (nrow(age_change)) {
  fwrite(age_change, file.path(out_dir, "age_varying_wealth_change_screen.csv"))
}
fwrite(phase_contrasts, file.path(out_dir, "recovery_phase_model_contrasts.csv"))
fwrite(phase_rate_ratios, file.path(out_dir, "recovery_phase_standardized_rate_ratios.csv"))
fwrite(phase_qc, file.path(out_dir, "recovery_phase_model_qc.csv"))

cat("Exploratory module point screens:\n")
print(point_screen[order(module, cohort, stratum)])
cat("\nSex difference screen:\n")
if (nrow(sex_wide)) {
  print(sex_wide[order(cohort)])
} else {
  cat("Not selected in this run.\n")
}
cat("\nLow-education wealth compensation screen:\n")
print(lowedu_wealth_screen[order(cohort)])
cat("\nAge-varying wealth screen:\n")
if (nrow(age_change)) {
  print(age_change[order(process, origin_state, cohort)])
} else {
  cat("Not selected in this run.\n")
}
cat("\nRecovery-phase contrasts:\n")
print(phase_contrasts[order(contrast, cohort)])

writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, paste0("sessionInfo_exploratory_upgrade_point_models_", sensitivity, ".txt"))
)
cat("Completed: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
