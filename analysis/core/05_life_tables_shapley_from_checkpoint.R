#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
analysis_tag <- Sys.getenv("D4_ANALYSIS_TAG", unset = "primary")
derived_dir <- Sys.getenv("D4_DERIVED_DIR", unset = file.path(root, "02_derived"))
out_dir <- Sys.getenv("D4_POINT_OUT_DIR", unset = file.path(root, "03_outputs", "04_formal_point"))
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("05_life_tables_shapley_", analysis_tag, "_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

rates_path <- file.path(out_dir, "standardized_annual_hazards_checkpoint.csv")
coef_path <- file.path(out_dir, "model_coefficients_checkpoint.csv")
qc_path <- file.path(out_dir, "model_convergence_qc_checkpoint.csv")
stopifnot(file.exists(rates_path), file.exists(coef_path), file.exists(qc_path))

rates <- fread(rates_path)
coefficients <- fread(coef_path)
model_qc <- fread(qc_path)

ages <- 60:99
living_states <- c("I0", "D1", "R1", "D2")
states <- c(living_states, "DEAD")
ses_levels <- c("high", "middle", "low")
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
      v <- a[, origin_name, destination_name]
      if (any(!is.finite(v))) stop("Missing hazard for ", origin_name, " -> ", destination_name)
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
  if (max(abs(rowSums(m) - 1)) > 1e-12 || any(m < -1e-12)) {
    stop("Transition-matrix coherence failure")
  }
  m
}

life_table <- function(lambda, return_trace = FALSE, return_matrices = FALSE) {
  v <- setNames(c(1, 0, 0, 0, 0), states)
  occupancy <- setNames(rep(0, length(states)), states)
  trace_parts <- list()
  matrix_parts <- list()
  for (age_value in ages) {
    m <- transition_matrix(lambda[as.character(age_value), , , drop = TRUE])
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
    if (return_matrices) {
      matrix_long <- as.data.table(as.table(m))
      setnames(matrix_long, names(matrix_long)[1:3], c("origin_state", "destination", "annual_probability"))
      matrix_long[, `:=`(
        age = age_value,
        origin_state = as.character(origin_state),
        destination = as.character(destination)
      )]
      setcolorder(matrix_long, c("age", "origin_state", "destination", "annual_probability"))
      matrix_parts[[length(matrix_parts) + 1L]] <- matrix_long
    }
    v <- v_next
  }
  metrics <- c(
    total_life_expectancy = sum(occupancy[living_states]),
    dfle = sum(occupancy[c("I0", "R1")]),
    disabled_years = sum(occupancy[c("D1", "D2")]),
    residual_alive_age100 = sum(v[living_states])
  )
  list(
    metrics = metrics,
    trace = if (return_trace) rbindlist(trace_parts) else NULL,
    matrices = if (return_matrices) rbindlist(matrix_parts) else NULL
  )
}

array_store <- list()
life_parts <- list()
trace_parts <- list()
matrix_parts <- list()

for (cohort_name in unique(rates$cohort)) {
  for (exposure_name in unique(rates$exposure)) {
    for (ses_name in ses_levels) {
      z <- rates[cohort == cohort_name & exposure == exposure_name & ses == ses_name]
      if (!nrow(z)) next
      key <- paste(cohort_name, exposure_name, ses_name, sep = "__")
      lambda <- lambda_array(z)
      array_store[[key]] <- lambda
      lt <- life_table(lambda, return_trace = TRUE, return_matrices = TRUE)
      life_parts[[length(life_parts) + 1L]] <- data.table(
        cohort = cohort_name,
        exposure = exposure_name,
        ses = ses_name,
        total_life_expectancy = unname(lt$metrics["total_life_expectancy"]),
        dfle = unname(lt$metrics["dfle"]),
        disabled_years = unname(lt$metrics["disabled_years"]),
        residual_alive_age100 = unname(lt$metrics["residual_alive_age100"]),
        state_years_closure = unname(
          lt$metrics["total_life_expectancy"] - lt$metrics["dfle"] - lt$metrics["disabled_years"]
        )
      )
      lt$trace[, `:=`(cohort = cohort_name, exposure = exposure_name, ses = ses_name)]
      lt$matrices[, `:=`(cohort = cohort_name, exposure = exposure_name, ses = ses_name)]
      trace_parts[[length(trace_parts) + 1L]] <- lt$trace
      matrix_parts[[length(matrix_parts) + 1L]] <- lt$matrices
    }
  }
}

life_expectancy <- rbindlist(life_parts)
occupancy_trace <- rbindlist(trace_parts)
annual_probabilities <- rbindlist(matrix_parts)

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
    for (tail in permutations(x[-i])) out[[length(out) + 1L]] <- c(x[[i]], tail)
  }
  out
}
block_permutations <- permutations(blocks)

replace_block <- function(base, donor, block_name) {
  out <- base
  for (pair in block_pairs[[block_name]]) out[, pair[[1]], pair[[2]]] <- donor[, pair[[1]], pair[[2]]]
  out
}

decomposition_parts <- list()
gap_parts <- list()
metric_names <- c("dfle", "disabled_years", "total_life_expectancy")

for (cohort_name in unique(rates$cohort)) {
  for (exposure_name in unique(rates$exposure)) {
    high_lambda <- array_store[[paste(cohort_name, exposure_name, "high", sep = "__")]]
    low_lambda <- array_store[[paste(cohort_name, exposure_name, "low", sep = "__")]]
    high_metrics <- life_table(high_lambda)$metrics
    low_metrics <- life_table(low_lambda)$metrics
    contributions <- matrix(0, nrow = length(blocks), ncol = length(metric_names),
      dimnames = list(blocks, metric_names))
    for (perm in block_permutations) {
      hybrid <- high_lambda
      previous <- high_metrics
      for (block_name in perm) {
        hybrid <- replace_block(hybrid, low_lambda, block_name)
        current <- life_table(hybrid)$metrics
        contributions[block_name, ] <- contributions[block_name, ] +
          current[metric_names] - previous[metric_names]
        previous <- current
      }
    }
    contributions <- contributions / length(block_permutations)
    gaps <- low_metrics[metric_names] - high_metrics[metric_names]
    for (metric_name in metric_names) {
      for (block_name in blocks) {
        value <- contributions[block_name, metric_name]
        decomposition_parts[[length(decomposition_parts) + 1L]] <- data.table(
          cohort = cohort_name,
          exposure = exposure_name,
          estimand = metric_name,
          block = block_name,
          contribution_years = value,
          total_low_minus_high_gap = unname(gaps[[metric_name]]),
          contribution_percent = if (abs(gaps[[metric_name]]) > 1e-10) {
            100 * value / gaps[[metric_name]]
          } else {
            NA_real_
          }
        )
      }
      gap_parts[[length(gap_parts) + 1L]] <- data.table(
        cohort = cohort_name,
        exposure = exposure_name,
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

decomposition <- rbindlist(decomposition_parts)
gaps <- rbindlist(gap_parts)

point_promotion_screen <- decomposition[
  estimand == "dfle" & block %in% c("recovery", "relapse"),
  .(
    recovery_relapse_contribution_years = sum(contribution_years),
    dfle_gap = unique(total_low_minus_high_gap)
  ),
  by = .(cohort, exposure)
]
point_promotion_screen[, recovery_relapse_percent := fifelse(
  abs(dfle_gap) > 1e-10,
  100 * recovery_relapse_contribution_years / dfle_gap,
  NA_real_
)]
point_promotion_screen[, point_threshold_met :=
  abs(recovery_relapse_contribution_years) >= 0.5 & abs(recovery_relapse_percent) >= 20]

country_heterogeneity_screen <- decomposition[estimand == "dfle", .(
  min_contribution_years = min(contribution_years),
  max_contribution_years = max(contribution_years),
  country_range_years = max(contribution_years) - min(contribution_years)
), by = .(exposure, block)]
country_heterogeneity_screen[, point_range_ge_0_75y := country_range_years >= 0.75]

if (max(abs(gaps$closure_error)) > 0.01) stop("Shapley closure exceeds 0.01 years")
if (max(abs(life_expectancy$state_years_closure)) > 1e-8) stop("Life-table state-year closure failure")

saveRDS(rates, file.path(derived_dir, "formal_point_annual_hazards.rds"), compress = "xz")
fwrite(model_qc, file.path(out_dir, "model_convergence_qc.csv"))
fwrite(coefficients, file.path(out_dir, "model_coefficients_model_based.csv"))
fwrite(rates, file.path(out_dir, "standardized_annual_hazards.csv"))
fwrite(annual_probabilities, file.path(out_dir, "annual_transition_probabilities.csv"))
fwrite(life_expectancy, file.path(out_dir, "life_expectancy_point_estimates.csv"))
fwrite(occupancy_trace, file.path(out_dir, "state_occupancy_trace.csv"))
fwrite(gaps, file.path(out_dir, "low_high_absolute_gaps.csv"))
fwrite(decomposition, file.path(out_dir, "shapley_decomposition_point.csv"))
fwrite(point_promotion_screen, file.path(out_dir, "point_promotion_screen.csv"))
fwrite(country_heterogeneity_screen, file.path(out_dir, "country_heterogeneity_point_screen.csv"))

cat("Life expectancy point estimates:\n")
print(life_expectancy[ses %in% c("high", "low")][order(cohort, exposure, ses)])
cat("\nDFLE decomposition:\n")
print(decomposition[estimand == "dfle"][order(cohort, exposure, block)])
cat("\nPoint promotion screen (bootstrap uncertainty still required):\n")
print(point_promotion_screen)
cat("\nPoint country heterogeneity screen (bootstrap uncertainty still required):\n")
print(country_heterogeneity_screen)

writeLines(capture.output(sessionInfo()), file.path(log_dir, paste0("sessionInfo_formal_point_postprocess_", analysis_tag, ".txt")))
cat("Completed: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
