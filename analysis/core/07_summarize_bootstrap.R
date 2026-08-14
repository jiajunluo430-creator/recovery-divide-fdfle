#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) tolower(args[[1]]) else "preview"
stopifnot(phase %in% c("preview", "final"))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
point_dir <- file.path(root, "03_outputs", "04_formal_point")
boot_dir <- file.path(root, "03_outputs", paste0("05_bootstrap_", phase))
summary_dir <- file.path(root, "03_outputs", paste0("06_bootstrap_summary_", phase))
log_dir <- file.path(root, "06_logs")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

target_replicates <- if (phase == "final") 500L else 100L
minimum_valid <- if (phase == "final") 450L else 90L
cohort_order <- c("CHARLS", "ELSA", "HRS", "MHAS")
exposure_order <- c("education", "wealth")

set_contract_order <- function(z, extra_columns = character()) {
  z[, `:=`(
    .cohort_order = match(cohort, cohort_order),
    .exposure_order = match(exposure, exposure_order)
  )]
  setorderv(z, c(".cohort_order", ".exposure_order", extra_columns))
  z[, c(".cohort_order", ".exposure_order") := NULL]
  invisible(z)
}

read_family <- function(suffix) {
  paths <- file.path(boot_dir, paste0("bootstrap_", tolower(cohort_order), suffix))
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths)) stop("Missing bootstrap files: ", paste(missing_paths, collapse = "; "))
  rbindlist(lapply(paths, fread), fill = TRUE)
}

replicates <- read_family("_replicates.csv")
qc <- read_family("_qc.csv")
setorder(replicates, cohort, exposure, replicate)
setorder(qc, cohort, exposure, replicate)

qc_summary <- qc[, .(
  target_replicates = target_replicates,
  qc_replicates = uniqueN(replicate),
  valid_replicates = uniqueN(replicate[status == "valid"]),
  failed_replicates = uniqueN(replicate[status != "valid"]),
  failure_messages = paste(unique(error_message[nzchar(error_message)]), collapse = " | "),
  median_unique_sampled_people = median(unique_sampled_people, na.rm = TRUE),
  median_elapsed_seconds = median(elapsed_seconds_replicate_exposure, na.rm = TRUE)
), by = .(cohort, exposure)]
result_counts <- replicates[, .(result_replicates = uniqueN(replicate)), by = .(cohort, exposure)]
qc_summary <- merge(qc_summary, result_counts, by = c("cohort", "exposure"), all.x = TRUE)
qc_summary[, `:=`(
  minimum_valid_required = minimum_valid,
  validity_gate_pass = valid_replicates >= minimum_valid & result_replicates >= minimum_valid
)]
set_contract_order(qc_summary)
fwrite(qc_summary, file.path(summary_dir, "bootstrap_validity_qc.csv"))

metric_columns <- c(
  "high_tle", "low_tle", "gap_tle",
  "high_dfle", "low_dfle", "gap_dfle",
  "high_disabled", "low_disabled", "gap_disabled",
  "high_residual_age100", "low_residual_age100",
  "contribution_onset", "contribution_recovery", "contribution_relapse",
  "contribution_post_disability_mortality", "contribution_pre_disability_mortality",
  "recovery_relapse_contribution", "recovery_relapse_percent", "closure_error"
)
long <- melt(
  replicates,
  id.vars = c("cohort", "phase", "replicate", "seed", "exposure"),
  measure.vars = metric_columns,
  variable.name = "metric", value.name = "estimate"
)
ci_summary <- long[is.finite(estimate), .(
  valid_replicates = .N,
  bootstrap_mean = mean(estimate),
  bootstrap_sd = sd(estimate),
  ci_low = as.numeric(quantile(estimate, 0.025, names = FALSE, type = 6)),
  bootstrap_median = median(estimate),
  ci_high = as.numeric(quantile(estimate, 0.975, names = FALSE, type = 6))
), by = .(cohort, exposure, metric)]

life <- fread(file.path(point_dir, "life_expectancy_point_estimates.csv"))
life_long <- melt(
  life[ses %in% c("high", "low")],
  id.vars = c("cohort", "exposure", "ses"),
  measure.vars = c("total_life_expectancy", "dfle", "disabled_years", "residual_alive_age100"),
  variable.name = "estimand", value.name = "point_estimate"
)
life_long[, metric := paste0(
  ses, "_",
  fifelse(estimand == "total_life_expectancy", "tle",
    fifelse(estimand == "disabled_years", "disabled",
      fifelse(estimand == "residual_alive_age100", "residual_age100", as.character(estimand))))
)]
point_parts <- list(life_long[, .(cohort, exposure, metric, point_estimate)])

gaps <- fread(file.path(point_dir, "low_high_absolute_gaps.csv"))
gap_map <- c(dfle = "gap_dfle", disabled_years = "gap_disabled", total_life_expectancy = "gap_tle")
gaps <- gaps[estimand %in% names(gap_map)]
gaps[, metric := unname(gap_map[estimand])]
point_parts[[length(point_parts) + 1L]] <- gaps[, .(cohort, exposure, metric, point_estimate = low_minus_high_gap)]

shapley <- fread(file.path(point_dir, "shapley_decomposition_point.csv"))[estimand == "dfle"]
shapley[, metric := paste0("contribution_", block)]
point_parts[[length(point_parts) + 1L]] <- shapley[, .(cohort, exposure, metric, point_estimate = contribution_years)]

promotion_point <- fread(file.path(point_dir, "point_promotion_screen.csv"))
point_parts[[length(point_parts) + 1L]] <- promotion_point[, .(
  cohort, exposure, metric = "recovery_relapse_contribution",
  point_estimate = recovery_relapse_contribution_years
)]
point_parts[[length(point_parts) + 1L]] <- promotion_point[, .(
  cohort, exposure, metric = "recovery_relapse_percent",
  point_estimate = recovery_relapse_percent
)]
point_parts[[length(point_parts) + 1L]] <- data.table(
  cohort = cohort_order, exposure = exposure_order[[1]], metric = "closure_error", point_estimate = NA_real_
)[0]

point_long <- unique(rbindlist(point_parts, fill = TRUE), by = c("cohort", "exposure", "metric"))
ci_summary <- merge(ci_summary, point_long, by = c("cohort", "exposure", "metric"), all.x = TRUE)
setcolorder(ci_summary, c(
  "cohort", "exposure", "metric", "point_estimate", "valid_replicates",
  "bootstrap_mean", "bootstrap_sd", "ci_low", "bootstrap_median", "ci_high"
))
set_contract_order(ci_summary, "metric")
fwrite(ci_summary, file.path(summary_dir, "bootstrap_percentile_intervals.csv"))

rr_ci <- ci_summary[metric == "recovery_relapse_contribution", .(
  cohort, exposure, rr_point = point_estimate, rr_ci_low = ci_low, rr_ci_high = ci_high,
  rr_valid_replicates = valid_replicates
)]
rr_pct_ci <- ci_summary[metric == "recovery_relapse_percent", .(
  cohort, exposure, rr_percent_point = point_estimate,
  rr_percent_ci_low = ci_low, rr_percent_ci_high = ci_high
)]
promotion <- merge(promotion_point, rr_ci, by = c("cohort", "exposure"), all.x = TRUE)
promotion <- merge(promotion, rr_pct_ci, by = c("cohort", "exposure"), all.x = TRUE)
promotion <- merge(
  promotion,
  qc_summary[, .(cohort, exposure, validity_gate_pass)],
  by = c("cohort", "exposure"), all.x = TRUE
)
promotion[, `:=`(
  rr_ci_excludes_zero = rr_ci_high < 0 | rr_ci_low > 0,
  rr_direction = fifelse(rr_point < 0, "widens_low_high_dfle_deficit",
    fifelse(rr_point > 0, "narrows_low_high_dfle_deficit", "null"))
)]
promotion[, mechanism_trigger_cohort := validity_gate_pass & point_threshold_met & rr_ci_excludes_zero]
set_contract_order(promotion)
fwrite(promotion, file.path(summary_dir, "mechanism_promotion_by_cohort.csv"))

block_metric_map <- c(
  onset = "contribution_onset",
  recovery = "contribution_recovery",
  relapse = "contribution_relapse",
  post_disability_mortality = "contribution_post_disability_mortality",
  pre_disability_mortality = "contribution_pre_disability_mortality"
)
heterogeneity_parts <- list()
for (exposure_name in exposure_order) {
  for (block_name in names(block_metric_map)) {
    metric_name <- unname(block_metric_map[[block_name]])
    z <- dcast(
      long[exposure == exposure_name & metric == metric_name],
      replicate ~ cohort, value.var = "estimate"
    )
    point_z <- shapley[exposure == exposure_name & block == block_name]
    for (pair in combn(cohort_order, 2L, simplify = FALSE)) {
      if (!all(pair %in% names(z))) next
      diff_values <- z[[pair[[1]]]] - z[[pair[[2]]]]
      diff_values <- diff_values[is.finite(diff_values)]
      if (!length(diff_values)) next
      point_a <- point_z[cohort == pair[[1]], contribution_years]
      point_b <- point_z[cohort == pair[[2]], contribution_years]
      point_diff <- point_a - point_b
      ci <- as.numeric(quantile(diff_values, c(0.025, 0.975), names = FALSE, type = 6))
      heterogeneity_parts[[length(heterogeneity_parts) + 1L]] <- data.table(
        exposure = exposure_name, block = block_name,
        cohort_a = pair[[1]], cohort_b = pair[[2]],
        point_difference_a_minus_b = point_diff,
        valid_paired_replicates = length(diff_values),
        ci_low = ci[[1]], ci_high = ci[[2]],
        absolute_point_difference_ge_0_75y = abs(point_diff) >= 0.75,
        ci_excludes_zero = ci[[2]] < 0 | ci[[1]] > 0,
        heterogeneity_trigger_pair = abs(point_diff) >= 0.75 & (ci[[2]] < 0 | ci[[1]] > 0) &
          length(diff_values) >= minimum_valid
      )
    }
  }
}
heterogeneity <- rbindlist(heterogeneity_parts, fill = TRUE)
setorder(heterogeneity, exposure, block, cohort_a, cohort_b)
fwrite(heterogeneity, file.path(summary_dir, "pairwise_country_heterogeneity.csv"))

mechanism_trigger <- promotion[mechanism_trigger_cohort == TRUE, .N, by = .(exposure, rr_direction)][N >= 2L]
mechanism_trigger_met <- nrow(mechanism_trigger) > 0L
heterogeneity_trigger_met <- any(heterogeneity$heterogeneity_trigger_pair, na.rm = TRUE)
decision <- data.table(
  phase = phase,
  target_replicates = target_replicates,
  minimum_valid_required = minimum_valid,
  all_cohort_exposure_validity_gates_pass = all(qc_summary$validity_gate_pass),
  mechanism_same_direction_two_cohorts_met = mechanism_trigger_met,
  qualifying_mechanism_strata = if (mechanism_trigger_met) paste(
    mechanism_trigger[, paste(exposure, rr_direction, paste0("n=", N), sep = ":")], collapse = " | "
  ) else "",
  country_heterogeneity_trigger_met = heterogeneity_trigger_met,
  qualifying_heterogeneity_pairs = if (heterogeneity_trigger_met) paste(
    heterogeneity[heterogeneity_trigger_pair == TRUE,
      paste(exposure, block, paste0(cohort_a, "-", cohort_b), sep = ":")], collapse = " | "
  ) else "",
  top_journal_promotion_gate = all(qc_summary$validity_gate_pass) &
    (mechanism_trigger_met | heterogeneity_trigger_met),
  decision_scope = if (phase == "final") "binding_final" else "computational_preview_only"
)
fwrite(decision, file.path(summary_dir, "bootstrap_promotion_decision.csv"))

max_closure <- replicates[, max(abs(closure_error), na.rm = TRUE)]
writeLines(capture.output(sessionInfo()), file.path(log_dir, paste0("sessionInfo_bootstrap_summary_", phase, ".txt")))
cat("Bootstrap summary completed\n")
cat("phase=", phase, " target=", target_replicates, " minimum_valid=", minimum_valid, "\n", sep = "")
cat("all_validity_gates=", all(qc_summary$validity_gate_pass), "\n", sep = "")
cat("mechanism_trigger=", mechanism_trigger_met, " heterogeneity_trigger=", heterogeneity_trigger_met, "\n", sep = "")
cat("top_journal_promotion_gate=", decision$top_journal_promotion_gate, " scope=", decision$decision_scope, "\n", sep = "")
cat("max_abs_shapley_closure=", format(max_closure, scientific = TRUE), "\n", sep = "")
