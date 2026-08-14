#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
analysis_tag <- Sys.getenv("D4_ANALYSIS_TAG", unset = "primary")
derived_dir <- Sys.getenv("D4_DERIVED_DIR", unset = file.path(root, "02_derived"))
out_dir <- Sys.getenv("D4_POINT_OUT_DIR", unset = file.path(root, "03_outputs", "04_formal_point"))
log_dir <- file.path(root, "06_logs")
for (d in c(out_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("04_fit_formal_models_and_decompose_", analysis_tag, "_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

function_risk <- as.data.table(readRDS(file.path(derived_dir, "formal_function_transition_riskset.rds")))
mortality_risk <- as.data.table(readRDS(file.path(derived_dir, "formal_mortality_riskset.rds")))

states <- c("I0", "D1", "R1", "D2", "DEAD")
ses_levels <- c("high", "middle", "low")
ages <- 60:99

process_specs <- list(
  onset = list(source = "function", origins = "I0", destination = "D1", event = "event_onset"),
  recovery = list(source = "function", origins = c("D1", "D2"), destination = "R1", event = "event_recovery"),
  relapse = list(source = "function", origins = "R1", destination = "D2", event = "event_relapse"),
  death_pre = list(source = "mortality", origins = "I0", destination = "DEAD", event = "event_death"),
  death_post = list(source = "mortality", origins = c("D1", "R1", "D2"), destination = "DEAD", event = "event_death")
)

function_risk[, event_recovery := as.integer(origin_state %in% c("D1", "D2") & destination == "R1")]

prepare_model_data <- function(cohort_name, exposure_name, spec) {
  level_col <- if (exposure_name == "education") "education3_fixed" else "wealth3"
  if (spec$source == "function") {
    z <- copy(function_risk[cohort == cohort_name & origin_state %in% spec$origins])
    z[, fit_weight_raw := respondent_weight]
  } else {
    z <- copy(mortality_risk[
      cohort == cohort_name & primary_mortality_interval == TRUE & origin_state %in% spec$origins
    ])
    z[, fit_weight_raw := origin_weight]
  }
  z[, event := as.integer(get(spec$event))]
  z[, ses := factor(get(level_col), levels = ses_levels)]
  z[, sex_factor := factor(fifelse(female == 1, "female",
    fifelse(female == 0, "male", "unknown")))]
  z[, proxy_factor := factor(fifelse(is.na(proxy), "unavailable",
    fifelse(proxy == 1, "proxy", "self")))]
  z[, period_factor := factor(origin_wave)]
  z[, origin_factor := factor(origin_state, levels = spec$origins)]
  z[, age_model := pmin(pmax(origin_age, 60), 95)]
  z <- z[
    !is.na(ses) & is.finite(interval_years) & interval_years > 0 &
    is.finite(fit_weight_raw) & fit_weight_raw > 0 & !is.na(age_model)
  ]
  z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
  z[, `:=`(
    sex_factor = droplevels(sex_factor),
    proxy_factor = droplevels(proxy_factor),
    period_factor = droplevels(period_factor),
    origin_factor = droplevels(origin_factor),
    ses = droplevels(ses)
  )]
  z
}

make_formula <- function(z) {
  terms <- c(
    "ses",
    "ns(age_model, knots = c(70, 80, 90), Boundary.knots = c(60, 95))"
  )
  if (nlevels(z$origin_factor) > 1L) terms <- c(terms, "origin_factor")
  if (nlevels(z$sex_factor) > 1L) terms <- c(terms, "sex_factor")
  if (nlevels(z$period_factor) > 1L) terms <- c(terms, "period_factor")
  if (nlevels(z$proxy_factor) > 1L) terms <- c(terms, "proxy_factor")
  as.formula(paste("event ~", paste(terms, collapse = " + "), "+ offset(log(interval_years))"))
}

weighted_mean_finite <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

fit_one_process <- function(cohort_name, exposure_name, process_name, spec) {
  z <- prepare_model_data(cohort_name, exposure_name, spec)
  if (!all(c("high", "low") %in% levels(z$ses))) {
    stop(cohort_name, " ", exposure_name, " ", process_name, ": low/high SES level absent")
  }
  event_by_level <- z[, .(rows = .N, events = sum(event)), by = ses]
  if (any(event_by_level[ses %in% c("high", "low"), rows == 0])) {
    stop(cohort_name, " ", exposure_name, " ", process_name, ": empty low/high risk set")
  }
  formula <- make_formula(z)
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
    stop(cohort_name, " ", exposure_name, " ", process_name, ": model did not converge with finite coefficients")
  }

  coef_mat <- summary(fit)$coefficients
  coef_dt <- data.table(
    cohort = cohort_name,
    exposure = exposure_name,
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

  qc <- data.table(
    cohort = cohort_name,
    exposure = exposure_name,
    process = process_name,
    rows = nrow(z),
    persons = uniqueN(z$person_id),
    events = sum(z$event),
    low_rows = sum(z$ses == "low"),
    low_events = sum(z$event[z$ses == "low"]),
    high_rows = sum(z$ses == "high"),
    high_events = sum(z$event[z$ses == "high"]),
    converged = fit$converged,
    iterations = fit$iter,
    coefficient_abs_max = max(abs(coef(fit))),
    warning_n = length(unique(warnings)),
    warnings = paste(unique(warnings), collapse = " | "),
    formula = paste(deparse(formula), collapse = " ")
  )

  rate_parts <- list()
  for (origin_name in spec$origins) {
    ref <- z[origin_state == origin_name, .(
      fit_weight = sum(fit_weight),
      standardization_rows = .N
    ), by = .(sex_factor, period_factor, proxy_factor, origin_factor)]
    if (!nrow(ref)) stop("No reference rows for ", cohort_name, " ", process_name, " ", origin_name)
    ref[, `:=`(
      sex_factor = factor(sex_factor, levels = levels(z$sex_factor)),
      period_factor = factor(period_factor, levels = levels(z$period_factor)),
      proxy_factor = factor(proxy_factor, levels = levels(z$proxy_factor)),
      origin_factor = factor(origin_factor, levels = levels(z$origin_factor))
    )]
    for (ses_name in ses_levels) {
      if (!ses_name %in% levels(z$ses)) next
      for (age in ages) {
        eval_age <- min(age, 95)
        nd <- copy(ref)
        nd[, `:=`(
          ses = factor(ses_name, levels = levels(z$ses)),
          age_model = eval_age,
          interval_years = 1
        )]
        link <- as.numeric(predict(fit, newdata = nd, type = "link"))
        hazard <- exp(link)
        rate_parts[[length(rate_parts) + 1L]] <- data.table(
          cohort = cohort_name,
          exposure = exposure_name,
          ses = ses_name,
          age = age,
          age_rate_evaluated = eval_age,
          tail_rate_held_at_95 = age > 95,
          process = process_name,
          origin_state = origin_name,
          destination = spec$destination,
          annual_hazard = weighted_mean_finite(hazard, nd$fit_weight),
          standardization_rows = sum(nd$standardization_rows)
        )
      }
    }
  }

  list(
    fit = fit,
    coefficients = coef_dt,
    qc = qc,
    rates = rbindlist(rate_parts)
  )
}

model_store <- list()
coef_parts <- list()
qc_parts <- list()
rate_parts <- list()

for (cohort_name in c("CHARLS", "HRS", "ELSA", "MHAS")) {
  for (exposure_name in c("education", "wealth")) {
    for (process_name in names(process_specs)) {
      cat("Fitting ", cohort_name, " / ", exposure_name, " / ", process_name, "...\n", sep = "")
      ans <- fit_one_process(cohort_name, exposure_name, process_name, process_specs[[process_name]])
      key <- paste(cohort_name, exposure_name, process_name, sep = "__")
      model_store[[key]] <- list(
        coefficients = coef(ans$fit),
        formula = formula(ans$fit),
        xlevels = ans$fit$xlevels,
        contrasts = ans$fit$contrasts,
        converged = ans$fit$converged,
        iterations = ans$fit$iter
      )
      coef_parts[[length(coef_parts) + 1L]] <- ans$coefficients
      qc_parts[[length(qc_parts) + 1L]] <- ans$qc
      rate_parts[[length(rate_parts) + 1L]] <- ans$rates
      cat("Completed ", cohort_name, " / ", exposure_name, " / ", process_name, ".\n", sep = "")
      rm(ans)
      invisible(gc())
    }
  }
}

coefficients <- rbindlist(coef_parts, fill = TRUE)
model_qc <- rbindlist(qc_parts, fill = TRUE)
rates <- rbindlist(rate_parts, fill = TRUE)
cat("All 40 models and standardized hazards completed.\n")
fwrite(model_qc, file.path(out_dir, "model_convergence_qc_checkpoint.csv"))
fwrite(coefficients, file.path(out_dir, "model_coefficients_checkpoint.csv"))
fwrite(rates, file.path(out_dir, "standardized_annual_hazards_checkpoint.csv"))

# The array-based implementation in 05_life_tables_shapley_from_checkpoint.R is
# the audited production postprocessor. Keep the older row-filter implementation
# below only as an explicitly requested legacy cross-check; it is much slower and
# is not part of the formal execution path.
run_legacy_postprocess <- "--run-legacy-postprocess" %in% commandArgs(trailingOnly = TRUE)
if (!run_legacy_postprocess) {
  saveRDS(model_store, file.path(derived_dir, "formal_point_models_internal.rds"), compress = "xz")
  writeLines(capture.output(sessionInfo()), file.path(log_dir, paste0("sessionInfo_formal_model_checkpoint_", analysis_tag, ".txt")))
  cat("Model checkpoints written. Continue with 05_life_tables_shapley_from_checkpoint.R.\n")
  quit(save = "no", status = 0L)
}

transition_matrix_from_rates <- function(rate_slice) {
  m <- matrix(0, nrow = length(states), ncol = length(states), dimnames = list(states, states))
  m["DEAD", "DEAD"] <- 1
  allowed <- list(
    I0 = c("D1", "DEAD"),
    D1 = c("R1", "DEAD"),
    R1 = c("D2", "DEAD"),
    D2 = c("R1", "DEAD")
  )
  for (origin_name in names(allowed)) {
    destinations <- allowed[[origin_name]]
    lambda <- vapply(destinations, function(destination_name) {
      z <- rate_slice[origin_state == origin_name & destination == destination_name, annual_hazard]
      if (length(z) != 1L || !is.finite(z)) {
        stop("Missing/nonunique hazard for ", origin_name, " -> ", destination_name)
      }
      z
    }, numeric(1))
    total <- sum(lambda)
    move <- if (total > 0) 1 - exp(-total) else 0
    m[origin_name, origin_name] <- 1 - move
    if (total > 0) {
      for (j in seq_along(destinations)) {
        m[origin_name, destinations[[j]]] <- move * lambda[[j]] / total
      }
    }
  }
  if (max(abs(rowSums(m) - 1)) > 1e-10 || any(m < -1e-12)) {
    stop("Transition-matrix coherence failure")
  }
  m
}

life_table <- function(rate_set, return_trace = FALSE) {
  v <- setNames(c(1, 0, 0, 0, 0), states)
  occupancy <- setNames(rep(0, length(states)), states)
  trace_parts <- list()
  for (age_value in ages) {
    m <- transition_matrix_from_rates(rate_set[age == age_value])
    v_next <- as.numeric(v %*% m)
    names(v_next) <- states
    half <- (v + v_next) / 2
    occupancy <- occupancy + half
    if (return_trace) {
      trace_parts[[length(trace_parts) + 1L]] <- data.table(
        age = age_value,
        state = states,
        probability_start = as.numeric(v),
        probability_end = as.numeric(v_next),
        half_cycle_occupancy = as.numeric(half)
      )
    }
    v <- v_next
  }
  metrics <- c(
    total_life_expectancy = sum(occupancy[c("I0", "D1", "R1", "D2")]),
    dfle = sum(occupancy[c("I0", "R1")]),
    disabled_years = sum(occupancy[c("D1", "D2")]),
    residual_alive_age100 = sum(v[c("I0", "D1", "R1", "D2")])
  )
  list(metrics = metrics, trace = if (return_trace) rbindlist(trace_parts) else NULL)
}

life_parts <- list()
trace_parts <- list()
for (cohort_name in unique(rates$cohort)) {
  for (exposure_name in unique(rates$exposure)) {
    for (ses_name in ses_levels) {
      z <- rates[cohort == cohort_name & exposure == exposure_name & ses == ses_name]
      if (!nrow(z)) next
      lt <- life_table(z, return_trace = TRUE)
      life_parts[[length(life_parts) + 1L]] <- data.table(
        cohort = cohort_name,
        exposure = exposure_name,
        ses = ses_name,
        total_life_expectancy = unname(lt$metrics["total_life_expectancy"]),
        dfle = unname(lt$metrics["dfle"]),
        disabled_years = unname(lt$metrics["disabled_years"]),
        residual_alive_age100 = unname(lt$metrics["residual_alive_age100"]),
        state_years_closure = unname(lt$metrics["total_life_expectancy"] -
          lt$metrics["dfle"] - lt$metrics["disabled_years"])
      )
      lt$trace[, `:=`(cohort = cohort_name, exposure = exposure_name, ses = ses_name)]
      trace_parts[[length(trace_parts) + 1L]] <- lt$trace
    }
  }
}
life_expectancy <- rbindlist(life_parts)
occupancy_trace <- rbindlist(trace_parts)

block_map <- data.table(
  origin_state = c("I0", "D1", "D2", "R1", "D1", "R1", "D2", "I0"),
  destination = c("D1", "R1", "R1", "D2", "DEAD", "DEAD", "DEAD", "DEAD"),
  block = c("onset", "recovery", "recovery", "relapse",
    "post_disability_mortality", "post_disability_mortality", "post_disability_mortality",
    "pre_disability_mortality")
)
blocks <- c("onset", "recovery", "relapse", "post_disability_mortality", "pre_disability_mortality")

permutations <- function(x) {
  if (length(x) == 1L) return(list(x))
  out <- list()
  for (i in seq_along(x)) {
    rest <- x[-i]
    for (tail in permutations(rest)) out[[length(out) + 1L]] <- c(x[[i]], tail)
  }
  out
}
block_permutations <- permutations(blocks)

replace_block <- function(base, donor, block_name) {
  pairs <- block_map[block == block_name]
  out <- copy(base)
  for (i in seq_len(nrow(pairs))) {
    o <- pairs$origin_state[[i]]
    d <- pairs$destination[[i]]
    replacement <- donor[origin_state == o & destination == d, annual_hazard]
    if (length(replacement) != length(ages)) stop("Shapley donor rate length failure")
    out[origin_state == o & destination == d, annual_hazard := replacement]
  }
  out
}

decomposition_parts <- list()
gap_parts <- list()

for (cohort_name in unique(rates$cohort)) {
  for (exposure_name in unique(rates$exposure)) {
    high_rates <- copy(rates[cohort == cohort_name & exposure == exposure_name & ses == "high"])
    low_rates <- copy(rates[cohort == cohort_name & exposure == exposure_name & ses == "low"])
    high_metrics <- life_table(high_rates)$metrics
    low_metrics <- life_table(low_rates)$metrics
    metric_names <- c("dfle", "disabled_years", "total_life_expectancy")
    contribution <- matrix(0, nrow = length(blocks), ncol = length(metric_names),
      dimnames = list(blocks, metric_names))
    for (perm in block_permutations) {
      hybrid <- copy(high_rates)
      previous <- high_metrics
      for (block_name in perm) {
        hybrid <- replace_block(hybrid, low_rates, block_name)
        current <- life_table(hybrid)$metrics
        contribution[block_name, ] <- contribution[block_name, ] +
          current[metric_names] - previous[metric_names]
        previous <- current
      }
    }
    contribution <- contribution / length(block_permutations)
    gaps <- low_metrics[metric_names] - high_metrics[metric_names]
    for (metric_name in metric_names) {
      for (block_name in blocks) {
        value <- contribution[block_name, metric_name]
        decomposition_parts[[length(decomposition_parts) + 1L]] <- data.table(
          cohort = cohort_name,
          exposure = exposure_name,
          estimand = metric_name,
          block = block_name,
          contribution_years = value,
          total_low_minus_high_gap = unname(gaps[[metric_name]]),
          contribution_percent = if (abs(gaps[[metric_name]]) > 1e-10) 100 * value / gaps[[metric_name]] else NA_real_
        )
      }
      gap_parts[[length(gap_parts) + 1L]] <- data.table(
        cohort = cohort_name,
        exposure = exposure_name,
        estimand = metric_name,
        high = unname(high_metrics[[metric_name]]),
        low = unname(low_metrics[[metric_name]]),
        low_minus_high_gap = unname(gaps[[metric_name]]),
        shapley_sum = sum(contribution[, metric_name]),
        closure_error = sum(contribution[, metric_name]) - unname(gaps[[metric_name]])
      )
    }
  }
}

decomposition <- rbindlist(decomposition_parts)
gaps <- rbindlist(gap_parts)

point_promotion_screen <- decomposition[estimand == "dfle" & block %in% c("recovery", "relapse"), .(
  recovery_relapse_contribution_years = sum(contribution_years),
  recovery_relapse_percent = 100 * sum(contribution_years) / unique(total_low_minus_high_gap),
  dfle_gap = unique(total_low_minus_high_gap)
), by = .(cohort, exposure)]
point_promotion_screen[, point_threshold_met :=
  abs(recovery_relapse_contribution_years) >= 0.5 & abs(recovery_relapse_percent) >= 20]

if (max(abs(gaps$closure_error)) > 0.01) stop("Shapley closure exceeds 0.01 years")
if (max(abs(life_expectancy$state_years_closure)) > 1e-8) stop("Life-table state-year closure failure")

saveRDS(model_store, file.path(derived_dir, "formal_point_models_internal.rds"), compress = "xz")
saveRDS(rates, file.path(derived_dir, "formal_point_annual_hazards.rds"), compress = "xz")
fwrite(model_qc, file.path(out_dir, "model_convergence_qc.csv"))
fwrite(coefficients, file.path(out_dir, "model_coefficients_model_based.csv"))
fwrite(rates, file.path(out_dir, "standardized_annual_hazards.csv"))
fwrite(life_expectancy, file.path(out_dir, "life_expectancy_point_estimates.csv"))
fwrite(occupancy_trace, file.path(out_dir, "state_occupancy_trace.csv"))
fwrite(gaps, file.path(out_dir, "low_high_absolute_gaps.csv"))
fwrite(decomposition, file.path(out_dir, "shapley_decomposition_point.csv"))
fwrite(point_promotion_screen, file.path(out_dir, "point_promotion_screen.csv"))

cat("\nModel QC:\n")
print(model_qc[, .(cohort, exposure, process, rows, persons, events, low_events, high_events, converged, coefficient_abs_max, warning_n)])
cat("\nLife expectancy point estimates:\n")
print(life_expectancy[ses %in% c("high", "low")][order(cohort, exposure, ses)])
cat("\nDFLE decomposition:\n")
print(decomposition[estimand == "dfle"][order(cohort, exposure, block)])
cat("\nPoint promotion screen (uncertainty still required):\n")
print(point_promotion_screen)

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_formal_point.txt"))
cat("Completed: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
