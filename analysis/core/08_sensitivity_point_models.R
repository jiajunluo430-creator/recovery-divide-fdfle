#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
derived_dir <- file.path(root, "02_derived")
point_dir <- file.path(root, "03_outputs", "04_formal_point")
out_dir <- file.path(root, "03_outputs", "07_sensitivity_point")
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("08_sensitivity_point_models_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

function_base <- as.data.table(readRDS(file.path(derived_dir, "formal_function_transition_riskset.rds")))
mortality_base <- as.data.table(readRDS(file.path(derived_dir, "formal_mortality_riskset.rds")))[
  primary_mortality_interval == TRUE
]
function_base[, event_recovery := as.integer(origin_state %in% c("D1", "D2") & destination == "R1")]

states <- c("I0", "D1", "R1", "D2", "DEAD")
living_states <- states[1:4]
ses_levels <- c("high", "middle", "low")
ages <- 60:99
cohort_order <- c("CHARLS", "ELSA", "HRS", "MHAS")
exposure_order <- c("education", "wealth")

set_contract_order <- function(z, extra_columns = character()) {
  z[, `:=`(
    .sensitivity_order = match(sensitivity, variants$sensitivity),
    .cohort_order = match(cohort, cohort_order),
    .exposure_order = match(exposure, exposure_order)
  )]
  setorderv(z, c(".sensitivity_order", ".cohort_order", ".exposure_order", extra_columns))
  z[, c(".sensitivity_order", ".cohort_order", ".exposure_order") := NULL]
  invisible(z)
}

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

# One-feature-at-a-time, contract-preserving point sensitivities. The primary
# reproduction is mandatory and must match the frozen production output.
variants <- data.table(
  sensitivity = c(
    "primary_reproduction",
    "unweighted",
    "interval_1_to_3_years",
    "exclude_explicit_proxy",
    "exclude_pandemic_crossing",
    "pre_pandemic_end_before_2020"
  ),
  weighted = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
  max_interval = c(Inf, Inf, 3, Inf, Inf, Inf),
  exclude_proxy = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
  exclude_pandemic_crossing = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
  pre_pandemic = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
)

apply_variant <- function(z, source_name, variant) {
  out <- copy(z)
  if (is.finite(variant$max_interval)) {
    out <- out[interval_years >= 1 & interval_years <= variant$max_interval]
  }
  if (isTRUE(variant$exclude_proxy)) {
    # CHARLS pre-2020 has no proxy indicator; those rows remain explicitly
    # "unavailable". This check removes only interviews known to be proxy.
    out <- out[is.na(proxy) | proxy != 1]
  }
  if (isTRUE(variant$exclude_pandemic_crossing)) {
    end_time <- if (source_name == "function") out$destination_time else out$endpoint_time
    crosses <- is.finite(out$origin_year) & is.finite(end_time) & out$origin_year < 2020 & end_time >= 2020
    out <- out[!crosses]
  }
  if (isTRUE(variant$pre_pandemic)) {
    end_time <- if (source_name == "function") out$destination_time else out$endpoint_time
    # A path-consistent pandemic sensitivity cannot remove only the bridge
    # interval and then retain later R1/D2 histories that were constructed
    # through that omitted observation. Restrict the entire risk set to
    # endpoints before 2020 instead.
    out <- out[is.finite(end_time) & end_time < 2020]
  }
  out
}

make_formula <- function(z) {
  terms <- c("ses", "ns(age_model, knots = c(70, 80, 90), Boundary.knots = c(60, 95))")
  if (nlevels(z$origin_factor) > 1L) terms <- c(terms, "origin_factor")
  if (nlevels(z$sex_factor) > 1L) terms <- c(terms, "sex_factor")
  if (nlevels(z$period_factor) > 1L) terms <- c(terms, "period_factor")
  if (nlevels(z$proxy_factor) > 1L) terms <- c(terms, "proxy_factor")
  as.formula(paste("event ~", paste(terms, collapse = " + "), "+ offset(log(interval_years))"))
}

prepare_model_data <- function(cohort_name, exposure_name, spec, variant) {
  if (spec$source == "function") {
    z <- apply_variant(function_base[cohort == cohort_name & origin_state %in% spec$origins], "function", variant)
    z[, fit_weight_raw := respondent_weight]
  } else {
    z <- apply_variant(mortality_base[cohort == cohort_name & origin_state %in% spec$origins], "mortality", variant)
    z[, fit_weight_raw := origin_weight]
  }
  z[, event := as.integer(get(spec$event))]
  z[, ses_value := if (exposure_name == "education") education3_fixed else wealth3]
  z[, ses := factor(ses_value, levels = ses_levels)]
  z[, sex_factor := factor(fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown")))]
  z[, proxy_factor := factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self")))]
  z[, period_factor := factor(origin_wave)]
  z[, origin_factor := factor(origin_state, levels = spec$origins)]
  z[, age_model := pmin(pmax(origin_age, 60), 95)]
  z <- z[
    !is.na(ses) & is.finite(interval_years) & interval_years > 0 &
      is.finite(fit_weight_raw) & fit_weight_raw > 0 & !is.na(age_model)
  ]
  if (isTRUE(variant$weighted)) {
    z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
  } else {
    z[, fit_weight := 1]
  }
  z[, `:=`(
    ses = droplevels(ses),
    sex_factor = droplevels(sex_factor),
    proxy_factor = droplevels(proxy_factor),
    period_factor = droplevels(period_factor),
    origin_factor = droplevels(origin_factor)
  )]
  z
}

fit_process_rates <- function(cohort_name, exposure_name, process_name, spec, variant) {
  z <- prepare_model_data(cohort_name, exposure_name, spec, variant)
  if (!all(c("high", "low") %in% levels(z$ses))) stop("low/high SES level absent")
  formula <- make_formula(z)
  warnings <- character()
  fit <- withCallingHandlers(
    glm(
      formula, data = z, family = binomial(link = "cloglog"), weights = fit_weight,
      control = glm.control(maxit = 100, epsilon = 1e-9)
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) stop("nonconvergent/nonfinite model")

  rate_parts <- list()
  for (origin_name in spec$origins) {
    ref <- z[origin_state == origin_name, .(
      fit_weight = sum(fit_weight), standardization_rows = .N
    ), by = .(sex_factor, period_factor, proxy_factor, origin_factor)]
    if (!nrow(ref)) stop("empty origin standardization set")
    ref[, `:=`(
      sex_factor = factor(sex_factor, levels = levels(z$sex_factor)),
      period_factor = factor(period_factor, levels = levels(z$period_factor)),
      proxy_factor = factor(proxy_factor, levels = levels(z$proxy_factor)),
      origin_factor = factor(origin_factor, levels = levels(z$origin_factor))
    )]
    for (ses_name in c("high", "low")) {
      for (age_value in ages) {
        nd <- copy(ref)
        nd[, `:=`(
          ses = factor(ses_name, levels = levels(z$ses)),
          age_model = min(age_value, 95), interval_years = 1
        )]
        hazard <- exp(as.numeric(predict(fit, newdata = nd, type = "link")))
        annual_hazard <- sum(hazard * nd$fit_weight) / sum(nd$fit_weight)
        if (!is.finite(annual_hazard) || annual_hazard < 0) stop("invalid standardized hazard")
        rate_parts[[length(rate_parts) + 1L]] <- data.table(
          cohort = cohort_name, exposure = exposure_name, ses = ses_name,
          age = age_value, process = process_name, origin_state = origin_name,
          destination = spec$destination, annual_hazard = annual_hazard
        )
      }
    }
  }
  qc <- data.table(
    cohort = cohort_name, exposure = exposure_name, process = process_name,
    rows = nrow(z), persons = uniqueN(z$person_id), events = sum(z$event),
    low_rows = sum(z$ses == "low"), low_events = sum(z$event[z$ses == "low"]),
    high_rows = sum(z$ses == "high"), high_events = sum(z$event[z$ses == "high"]),
    converged = fit$converged, iterations = fit$iter,
    coefficient_abs_max = max(abs(coef(fit))), warning_n = uniqueN(warnings),
    warnings = paste(unique(warnings), collapse = " | ")
  )
  coefficients <- data.table(
    cohort = cohort_name, exposure = exposure_name, process = process_name,
    term = names(coef(fit)), estimate = as.numeric(coef(fit)),
    hazard_ratio = exp(as.numeric(coef(fit)))
  )
  list(rates = rbindlist(rate_parts), qc = qc, coefficients = coefficients)
}

lambda_array <- function(z) {
  a <- array(
    NA_real_, c(length(ages), length(living_states), length(states)),
    dimnames = list(as.character(ages), living_states, states)
  )
  for (i in seq_len(nrow(z))) {
    a[as.character(z$age[[i]]), z$origin_state[[i]], z$destination[[i]]] <- z$annual_hazard[[i]]
  }
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
  if (max(abs(rowSums(m) - 1)) > 1e-10) stop("matrix row-sum failure")
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
  for (i in seq_along(x)) for (tail in permutations(x[-i])) out[[length(out) + 1L]] <- c(x[[i]], tail)
  out
}
block_permutations <- permutations(blocks)

replace_block <- function(base, donor, block_name) {
  out <- base
  for (pair in block_pairs[[block_name]]) out[, pair[[1]], pair[[2]]] <- donor[, pair[[1]], pair[[2]]]
  out
}

decompose_dfle <- function(high_lambda, low_lambda) {
  high <- life_metrics(high_lambda)
  low <- life_metrics(low_lambda)
  contribution <- setNames(rep(0, length(blocks)), blocks)
  for (perm in block_permutations) {
    hybrid <- high_lambda
    previous <- high[["dfle"]]
    for (block_name in perm) {
      hybrid <- replace_block(hybrid, low_lambda, block_name)
      current <- life_metrics(hybrid)[["dfle"]]
      contribution[[block_name]] <- contribution[[block_name]] + current - previous
      previous <- current
    }
  }
  contribution <- contribution / length(block_permutations)
  gap <- low[["dfle"]] - high[["dfle"]]
  list(high = high, low = low, gap = gap, contribution = contribution, closure = sum(contribution) - gap)
}

rate_parts <- list()
qc_parts <- list()
coefficient_parts <- list()
for (variant_i in seq_len(nrow(variants))) {
  variant <- variants[variant_i]
  cat("\nSensitivity: ", variant$sensitivity, "\n", sep = "")
  for (cohort_name in cohort_order) for (exposure_name in exposure_order) for (process_name in names(process_specs)) {
    cat("  ", cohort_name, " / ", exposure_name, " / ", process_name, "\n", sep = "")
    ans <- fit_process_rates(cohort_name, exposure_name, process_name, process_specs[[process_name]], variant)
    ans$rates[, sensitivity := variant$sensitivity]
    ans$qc[, sensitivity := variant$sensitivity]
    ans$coefficients[, sensitivity := variant$sensitivity]
    rate_parts[[length(rate_parts) + 1L]] <- ans$rates
    qc_parts[[length(qc_parts) + 1L]] <- ans$qc
    coefficient_parts[[length(coefficient_parts) + 1L]] <- ans$coefficients
    rm(ans)
    invisible(gc())
  }
}
rates <- rbindlist(rate_parts)
model_qc <- rbindlist(qc_parts)
coefficients <- rbindlist(coefficient_parts)
rate_qc <- rates[, .(
  finite_hazards = sum(is.finite(annual_hazard)),
  max_annual_hazard = max(annual_hazard, na.rm = TRUE),
  p99_annual_hazard = as.numeric(quantile(annual_hazard, 0.99, na.rm = TRUE, names = FALSE)),
  hazards_gt_10 = sum(annual_hazard > 10, na.rm = TRUE),
  hazards_gt_100 = sum(annual_hazard > 100, na.rm = TRUE)
), by = .(sensitivity, cohort, exposure, process)]

life_parts <- list()
decomposition_parts <- list()
for (sensitivity_name in variants$sensitivity) for (cohort_name in cohort_order) for (exposure_name in exposure_order) {
  z <- rates[sensitivity == sensitivity_name & cohort == cohort_name & exposure == exposure_name]
  high_lambda <- lambda_array(z[ses == "high"])
  low_lambda <- lambda_array(z[ses == "low"])
  ans <- decompose_dfle(high_lambda, low_lambda)
  life_parts[[length(life_parts) + 1L]] <- rbind(
    data.table(sensitivity = sensitivity_name, cohort = cohort_name, exposure = exposure_name, ses = "high", as.list(ans$high)),
    data.table(sensitivity = sensitivity_name, cohort = cohort_name, exposure = exposure_name, ses = "low", as.list(ans$low))
  )
  decomposition_parts[[length(decomposition_parts) + 1L]] <- data.table(
    sensitivity = sensitivity_name, cohort = cohort_name, exposure = exposure_name,
    block = names(ans$contribution), contribution_years = as.numeric(ans$contribution),
    dfle_gap_low_minus_high = ans$gap, closure_error = ans$closure
  )
}
life <- rbindlist(life_parts)
decomposition <- rbindlist(decomposition_parts)
summary <- decomposition[, .(
  dfle_gap_low_minus_high = unique(dfle_gap_low_minus_high),
  onset_contribution = contribution_years[block == "onset"],
  recovery_contribution = contribution_years[block == "recovery"],
  relapse_contribution = contribution_years[block == "relapse"],
  recovery_relapse_contribution = sum(contribution_years[block %in% c("recovery", "relapse")]),
  recovery_relapse_percent = 100 * sum(contribution_years[block %in% c("recovery", "relapse")]) / unique(dfle_gap_low_minus_high),
  max_abs_closure_error = max(abs(closure_error))
), by = .(sensitivity, cohort, exposure)]

# Exact primary reproduction audit against the frozen production outputs.
frozen_gap <- fread(file.path(point_dir, "low_high_absolute_gaps.csv"))[estimand == "dfle", .(
  cohort, exposure, frozen_dfle_gap = low_minus_high_gap
)]
frozen_decomp <- fread(file.path(point_dir, "shapley_decomposition_point.csv"))[estimand == "dfle", .(
  cohort, exposure, block, frozen_contribution = contribution_years
)]
repro_gap <- merge(
  summary[sensitivity == "primary_reproduction", .(cohort, exposure, reproduced_dfle_gap = dfle_gap_low_minus_high)],
  frozen_gap, by = c("cohort", "exposure")
)
repro_gap[, abs_difference := abs(reproduced_dfle_gap - frozen_dfle_gap)]
repro_decomp <- merge(
  decomposition[sensitivity == "primary_reproduction", .(cohort, exposure, block, reproduced_contribution = contribution_years)],
  frozen_decomp, by = c("cohort", "exposure", "block")
)
repro_decomp[, abs_difference := abs(reproduced_contribution - frozen_contribution)]
reproduction_qc <- data.table(
  component = c("dfle_gap", "shapley_contribution"),
  max_abs_difference = c(max(repro_gap$abs_difference), max(repro_decomp$abs_difference)),
  tolerance = 1e-8
)
reproduction_qc[, passed := max_abs_difference <= tolerance]
if (!all(reproduction_qc$passed)) stop("Primary reproduction check failed")

primary_summary <- summary[sensitivity == "primary_reproduction", .(
  cohort, exposure, primary_dfle_gap = dfle_gap_low_minus_high,
  primary_rr = recovery_relapse_contribution
)]
concordance <- merge(summary, primary_summary, by = c("cohort", "exposure"), all.x = TRUE)
concordance[, `:=`(
  dfle_gap_direction_concordant = sign(dfle_gap_low_minus_high) == sign(primary_dfle_gap),
  rr_direction_concordant = sign(recovery_relapse_contribution) == sign(primary_rr),
  dfle_gap_ratio_to_primary = dfle_gap_low_minus_high / primary_dfle_gap,
  rr_ratio_to_primary = recovery_relapse_contribution / primary_rr
)]

set_contract_order(model_qc, "process")
set_contract_order(coefficients, c("process", "term"))
set_contract_order(rate_qc, "process")
set_contract_order(life, "ses")
set_contract_order(decomposition, "block")
set_contract_order(concordance)

fwrite(variants, file.path(out_dir, "sensitivity_specifications.csv"))
fwrite(model_qc, file.path(out_dir, "sensitivity_model_qc.csv"))
fwrite(coefficients, file.path(out_dir, "sensitivity_model_coefficients.csv"))
fwrite(rates, file.path(out_dir, "sensitivity_standardized_annual_hazards.csv"))
fwrite(rate_qc, file.path(out_dir, "sensitivity_rate_support_qc.csv"))
fwrite(life, file.path(out_dir, "sensitivity_life_expectancy.csv"))
fwrite(decomposition, file.path(out_dir, "sensitivity_shapley_decomposition.csv"))
fwrite(summary, file.path(out_dir, "sensitivity_summary.csv"))
fwrite(repro_gap, file.path(out_dir, "primary_reproduction_gap_detail.csv"))
fwrite(repro_decomp, file.path(out_dir, "primary_reproduction_decomposition_detail.csv"))
fwrite(reproduction_qc, file.path(out_dir, "primary_reproduction_qc.csv"))
fwrite(concordance, file.path(out_dir, "sensitivity_direction_and_magnitude.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_sensitivity_point.txt"))

cat("\nPrimary reproduction QC:\n")
print(reproduction_qc)
cat("\nSensitivity summary:\n")
print(concordance)
cat("\nCompleted: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
