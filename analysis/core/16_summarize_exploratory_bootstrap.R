#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) tolower(args[[1L]]) else "preview"
stopifnot(phase %in% c("test", "preview", "final"))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
point_dir <- file.path(root, "03_outputs", "11_exploratory_upgrade_point")
boot_dir <- file.path(root, "03_outputs", paste0("12_exploratory_bootstrap_", phase))
summary_dir <- file.path(root, "03_outputs", paste0("13_exploratory_bootstrap_summary_", phase))
log_dir <- file.path(root, "06_logs")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cohort_order <- c("CHARLS", "ELSA", "HRS", "MHAS")
module_order <- c("wealth_within_low_education", "recovery_phase", "wealth_by_sex")
target_replicates <- switch(phase, test = 2L, preview = 100L, final = 500L)
minimum_valid <- switch(phase, test = 2L, preview = 90L, final = 450L)
decision_scope <- switch(
  phase,
  test = "smoke_test_only_do_not_interpret_intervals",
  preview = "computational_preview_only",
  final = "binding_final"
)

set_order <- function(z, extra = character()) {
  if ("cohort" %in% names(z)) z[, .cohort_order := match(cohort, cohort_order)]
  if ("module" %in% names(z)) z[, .module_order := match(module, module_order)]
  ordering <- c(
    if (".module_order" %in% names(z)) ".module_order",
    if (".cohort_order" %in% names(z)) ".cohort_order",
    extra
  )
  if (length(ordering)) setorderv(z, ordering)
  drop <- intersect(c(".module_order", ".cohort_order"), names(z))
  if (length(drop)) z[, (drop) := NULL]
  invisible(z)
}

read_family <- function(suffix) {
  available <- list.files(boot_dir, full.names = TRUE)
  pieces <- lapply(cohort_order, function(cohort_name) {
    prefix <- paste0("bootstrap_", tolower(cohort_name))
    base_path <- file.path(boot_dir, paste0(prefix, suffix))
    shard_paths <- available[
      startsWith(basename(available), paste0(prefix, "_part_")) &
        endsWith(basename(available), suffix)
    ]
    paths <- unique(c(if (file.exists(base_path)) base_path else character(), shard_paths))
    if (!length(paths)) {
      stop("Missing bootstrap family for ", cohort_name, suffix)
    }
    rbindlist(lapply(paths, function(path) {
      value <- fread(path)
      value[, `:=`(
        source_file = basename(path),
        source_mtime = as.numeric(file.info(path)$mtime)
      )]
      value
    }), fill = TRUE)
  })
  rbindlist(pieces, fill = TRUE)
}

metrics <- read_family("_module_metrics.csv")
qc <- read_family("_qc.csv")
stopifnot(all(c("cohort", "replicate", "module", "stratum", "metric", "estimate") %in% names(metrics)))
stopifnot(all(c("cohort", "replicate", "module", "status") %in% names(qc)))

# A failed checkpoint can leave module results without the matching QC file. In
# addition, a later sharded resume can legitimately supersede an earlier row.
# Resolve duplicate keys by source modification time, then retain metrics only
# for module-replicates whose selected QC row is explicitly valid.
qc_key <- c("cohort", "replicate", "module")
metric_key <- c(qc_key, "stratum", "metric")

qc_overlap_audit <- qc[, .(
  source_rows = .N,
  distinct_status = uniqueN(status),
  distinct_seed = uniqueN(seed),
  source_files = paste(sort(unique(source_file)), collapse = " | ")
), by = qc_key][source_rows > 1L]
metric_overlap_audit <- metrics[, {
  finite_estimates <- estimate[is.finite(estimate)]
  list(
    source_rows = .N,
    estimate_min = if (length(finite_estimates)) min(finite_estimates) else NA_real_,
    estimate_max = if (length(finite_estimates)) max(finite_estimates) else NA_real_,
    nonfinite_rows = sum(!is.finite(estimate)),
    distinct_seed = uniqueN(seed),
    source_files = paste(sort(unique(source_file)), collapse = " | ")
  )
}, by = metric_key][source_rows > 1L]
metric_overlap_audit[, estimate_range := estimate_max - estimate_min]
metric_overlap_audit[, mixed_finite_nonfinite :=
  nonfinite_rows > 0L & nonfinite_rows < source_rows]

if (nrow(qc_overlap_audit) && any(
  qc_overlap_audit$distinct_status > 1L | qc_overlap_audit$distinct_seed > 1L
)) {
  stop("Conflicting QC rows across bootstrap sources")
}
if (nrow(metric_overlap_audit) && any(
  metric_overlap_audit$distinct_seed > 1L |
    metric_overlap_audit$mixed_finite_nonfinite |
    (!is.na(metric_overlap_audit$estimate_range) & metric_overlap_audit$estimate_range > 1e-10)
)) {
  stop("Conflicting metric rows across bootstrap sources")
}

setorderv(qc, c(qc_key, "source_mtime", "source_file"))
qc <- qc[, .SD[.N], by = qc_key]
setorderv(metrics, c(metric_key, "source_mtime", "source_file"))
metrics <- metrics[, .SD[.N], by = metric_key]
metrics_before_qc_filter <- nrow(metrics)
valid_qc_keys <- qc[status == "valid", ..qc_key]
metrics <- metrics[valid_qc_keys, on = qc_key, nomatch = 0L]
orphan_or_failed_metric_rows_removed <- metrics_before_qc_filter - nrow(metrics)
fwrite(qc_overlap_audit, file.path(summary_dir, "bootstrap_qc_source_overlap_audit.csv"))
fwrite(metric_overlap_audit, file.path(summary_dir, "bootstrap_metric_source_overlap_audit.csv"))
fwrite(data.table(
  metric_rows_before_qc_filter = metrics_before_qc_filter,
  metric_rows_retained = nrow(metrics),
  orphan_or_failed_metric_rows_removed = orphan_or_failed_metric_rows_removed
), file.path(summary_dir, "bootstrap_metric_qc_filter_audit.csv"))

# Public, non-disclosive consolidated replicate outputs. Source filenames and
# modification times are operational provenance only and are deliberately
# omitted; these rows contain cohort/module/replicate aggregate estimates and
# QC, never participant identifiers or resampled person frequencies.
qc_public <- copy(qc)
qc_public[, c("source_file", "source_mtime") := NULL]
metrics_public <- copy(metrics)
metrics_public[, c("source_file", "source_mtime") := NULL]
setorderv(qc_public, qc_key)
setorderv(metrics_public, metric_key)
fwrite(qc_public, file.path(summary_dir, "bootstrap_qc_consolidated.csv"))
fwrite(metrics_public, file.path(summary_dir, "bootstrap_valid_metrics_consolidated.csv"))
cat(
  "Bootstrap family merge: metric_rows_before_qc_filter=", metrics_before_qc_filter,
  " retained=", nrow(metrics),
  " removed=", orphan_or_failed_metric_rows_removed, "\n",
  sep = ""
)

qc_summary <- qc[, .(
  target_replicates = target_replicates,
  qc_replicates = uniqueN(replicate),
  valid_replicates = uniqueN(replicate[status == "valid"]),
  failed_replicates = uniqueN(replicate[status != "valid"]),
  failure_messages = paste(unique(na.omit(error_message[nzchar(error_message)])), collapse = " | "),
  median_unique_sampled_people = median(unique_sampled_people, na.rm = TRUE),
  median_elapsed_seconds = median(elapsed_seconds, na.rm = TRUE)
), by = .(cohort, module)]
result_counts <- metrics[, .(result_replicates = uniqueN(replicate)), by = .(cohort, module)]
qc_summary <- merge(qc_summary, result_counts, by = c("cohort", "module"), all.x = TRUE)
expected_qc_groups <- CJ(cohort = cohort_order, module = module_order, unique = TRUE)
expected_qc_keys <- expected_qc_groups[, paste(cohort, module, sep = "::")]
observed_qc_keys <- qc_summary[, paste(cohort, module, sep = "::")]
if (!setequal(expected_qc_keys, observed_qc_keys) || nrow(qc_summary) != nrow(expected_qc_groups)) {
  stop("Incomplete or unexpected cohort-module QC grid")
}
qc_summary[, `:=`(
  minimum_valid_required = minimum_valid,
  validity_gate_pass = qc_replicates == target_replicates &
    valid_replicates >= minimum_valid & result_replicates >= minimum_valid
)]
set_order(qc_summary)
fwrite(qc_summary, file.path(summary_dir, "bootstrap_validity_qc.csv"))

intervals <- metrics[is.finite(estimate), .(
  valid_replicates = .N,
  bootstrap_mean = mean(estimate),
  bootstrap_sd = sd(estimate),
  ci_low = as.numeric(quantile(estimate, 0.025, names = FALSE, type = 6)),
  bootstrap_median = median(estimate),
  ci_high = as.numeric(quantile(estimate, 0.975, names = FALSE, type = 6))
), by = .(cohort, module, stratum, metric)]

life <- fread(file.path(point_dir, "life_expectancy_point_estimates.csv"))
life_point <- melt(
  life[ses %in% c("high", "low")],
  id.vars = c("module", "stratum", "cohort", "ses"),
  measure.vars = "dfle", variable.name = "estimand", value.name = "point_estimate"
)
life_point[, metric := paste0(ses, "_dfle")]
life_point <- life_point[, .(cohort, module, stratum, metric, point_estimate)]

gaps <- fread(file.path(point_dir, "low_high_absolute_gaps.csv"))[estimand == "dfle"]
gap_point <- gaps[, .(
  cohort, module, stratum, metric = "gap_dfle", point_estimate = low_minus_high_gap
)]

shapley <- fread(file.path(point_dir, "shapley_decomposition_point.csv"))[estimand == "dfle"]
shapley_point <- shapley[, .(
  cohort, module, stratum, metric = paste0("contribution_", block),
  point_estimate = contribution_years
)]

promotion_point <- fread(file.path(point_dir, "module_point_promotion_screen.csv"))
rr_point <- rbind(
  promotion_point[, .(
    cohort, module, stratum, metric = "recovery_relapse_contribution",
    point_estimate = recovery_relapse_contribution_years
  )],
  promotion_point[, .(
    cohort, module, stratum, metric = "recovery_relapse_percent",
    point_estimate = recovery_relapse_percent
  )]
)

life_module_point <- unique(
  rbindlist(list(life_point, gap_point, shapley_point, rr_point), fill = TRUE),
  by = c("cohort", "module", "stratum", "metric")
)
sex_difference_point <- dcast(
  life_module_point[module == "wealth_by_sex" & stratum %in% c("male", "female")],
  cohort + module + metric ~ stratum,
  value.var = "point_estimate"
)
sex_difference_point[, `:=`(
  stratum = "female_minus_male",
  point_estimate = female - male
)]
sex_difference_point <- sex_difference_point[, .(cohort, module, stratum, metric, point_estimate)]

phase_point <- fread(file.path(point_dir, "recovery_phase_model_contrasts.csv"))
phase_strata <- c(
  low_vs_high_early_recovery = "early_recovery",
  low_vs_high_sustained_recovery = "sustained_recovery",
  sustained_vs_early_modification = "sustained_vs_early_modification"
)
phase_point[, `:=`(
  module = "recovery_phase",
  stratum = unname(phase_strata[contrast])
)]
phase_point <- melt(
  phase_point,
  id.vars = c("cohort", "module", "stratum"),
  measure.vars = c("log_hazard_ratio", "hazard_ratio", "std_error_model"),
  variable.name = "metric", value.name = "point_estimate"
)
phase_point[metric == "std_error_model", metric := "model_se"]

point_long <- unique(
  rbindlist(list(life_module_point, sex_difference_point, phase_point), fill = TRUE),
  by = c("cohort", "module", "stratum", "metric")
)
intervals <- merge(
  intervals, point_long,
  by = c("cohort", "module", "stratum", "metric"), all.x = TRUE
)
setcolorder(intervals, c(
  "cohort", "module", "stratum", "metric", "point_estimate", "valid_replicates",
  "bootstrap_mean", "bootstrap_sd", "ci_low", "bootstrap_median", "ci_high"
))
set_order(intervals, c("stratum", "metric"))
fwrite(intervals, file.path(summary_dir, "bootstrap_percentile_intervals.csv"))

validity_lookup <- qc_summary[, .(cohort, module, validity_gate_pass)]

lowedu_point <- fread(file.path(point_dir, "low_education_wealth_compensation_point_screen.csv"))[, .(
  cohort, module, stratum,
  recovery_relapse_contribution_years,
  dfle_gap,
  recovery_relapse_percent
)]
lowedu_ci <- dcast(
  intervals[module == "wealth_within_low_education" & stratum == "all" &
    metric %in% c("recovery_relapse_contribution", "gap_dfle")],
  cohort + module + stratum ~ metric,
  value.var = c("ci_low", "ci_high", "valid_replicates")
)
lowedu <- merge(lowedu_point, lowedu_ci, by = c("cohort", "module", "stratum"), all.x = TRUE)
lowedu <- merge(lowedu, validity_lookup, by = c("cohort", "module"), all.x = TRUE)
lowedu[, `:=`(
  metric_validity_gate_pass =
    valid_replicates_gap_dfle >= minimum_valid &
    valid_replicates_recovery_relapse_contribution >= minimum_valid,
  point_magnitude_adverse_ge_0_50y = recovery_relapse_contribution_years <= -0.50,
  bootstrap_interval_below_zero = ci_high_recovery_relapse_contribution < 0,
  cohort_promotion_gate = validity_gate_pass &
    valid_replicates_gap_dfle >= minimum_valid &
    valid_replicates_recovery_relapse_contribution >= minimum_valid &
    recovery_relapse_contribution_years <= -0.50 &
    ci_high_recovery_relapse_contribution < 0
)]
set_order(lowedu)
fwrite(lowedu, file.path(summary_dir, "low_education_wealth_promotion_by_cohort.csv"))
lowedu_same_direction_two <- lowedu[cohort_promotion_gate == TRUE, .N] >= 2L

phase_point_screen <- fread(file.path(point_dir, "recovery_phase_model_contrasts.csv"))
phase_point_screen[, `:=`(
  module = "recovery_phase",
  stratum = unname(phase_strata[contrast])
)]
phase_ci <- intervals[module == "recovery_phase" & metric == "hazard_ratio", .(
  cohort, module, stratum,
  bootstrap_hr_ci_low = ci_low,
  bootstrap_hr_ci_high = ci_high,
  valid_replicates
)]
phase_screen <- merge(
  phase_point_screen[, .(cohort, module, stratum, contrast, hazard_ratio)],
  phase_ci, by = c("cohort", "module", "stratum"), all.x = TRUE
)
phase_screen <- merge(phase_screen, validity_lookup, by = c("cohort", "module"), all.x = TRUE)
phase_screen[, `:=`(
  metric_validity_gate_pass = valid_replicates >= minimum_valid,
  point_threshold_met = fifelse(
    stratum == "sustained_vs_early_modification",
    hazard_ratio >= 1.20 | hazard_ratio <= 0.80,
    hazard_ratio >= 1.15
  ),
  bootstrap_interval_support = fifelse(
    stratum == "sustained_vs_early_modification",
    bootstrap_hr_ci_low > 1 | bootstrap_hr_ci_high < 1,
    bootstrap_hr_ci_low > 1
  )
)]
phase_screen[, cohort_phase_gate := validity_gate_pass & metric_validity_gate_pass &
  point_threshold_met & bootstrap_interval_support]
set_order(phase_screen, "stratum")
fwrite(phase_screen, file.path(summary_dir, "recovery_phase_promotion_by_cohort.csv"))
# The v2.1 addendum defines primary/secondary durability as replication of the
# low-versus-high wealth contrast in the same observed recovery phase.  The
# sustained-versus-early modification is a separate diagnostic and cannot by
# itself satisfy either same-phase replication rule.
phase_counts <- phase_screen[
  cohort_phase_gate == TRUE &
    stratum %in% c("early_recovery", "sustained_recovery"),
  .N, by = stratum
]
phase_same_phase_three <- any(phase_counts$N >= 3L)
phase_same_phase_two <- any(phase_counts$N >= 2L)

sex_point <- fread(file.path(point_dir, "sex_difference_point_screen.csv"))[, .(
  cohort,
  module = "wealth_by_sex",
  stratum = "female_minus_male",
  female_minus_male_recovery_relapse
)]
sex_ci <- intervals[
  module == "wealth_by_sex" & stratum == "female_minus_male" &
    metric == "recovery_relapse_contribution",
  .(
    cohort, module, stratum,
    bootstrap_ci_low = ci_low,
    bootstrap_ci_high = ci_high,
    valid_replicates
  )
]
sex_screen <- merge(sex_point, sex_ci, by = c("cohort", "module", "stratum"), all.x = TRUE)
sex_screen <- merge(sex_screen, validity_lookup, by = c("cohort", "module"), all.x = TRUE)
sex_screen[, `:=`(
  metric_validity_gate_pass = valid_replicates >= minimum_valid,
  direction = fifelse(
    female_minus_male_recovery_relapse < 0,
    "more_adverse_recovery_relapse_contribution_in_women",
    "more_adverse_recovery_relapse_contribution_in_men"
  ),
  magnitude_ge_0_50y = abs(female_minus_male_recovery_relapse) >= 0.50,
  magnitude_ge_0_75y = abs(female_minus_male_recovery_relapse) >= 0.75,
  bootstrap_interval_excludes_zero = bootstrap_ci_low > 0 | bootstrap_ci_high < 0
)]
sex_screen[, interval_supported := validity_gate_pass & metric_validity_gate_pass &
  bootstrap_interval_excludes_zero]
set_order(sex_screen)
fwrite(sex_screen, file.path(summary_dir, "sex_heterogeneity_promotion_by_cohort.csv"))

sex_route1 <- sex_screen[
  interval_supported & magnitude_ge_0_50y,
  .N, by = direction
][N >= 2L, .N] > 0L
sex_route2 <- FALSE
for (direction_name in unique(sex_screen$direction)) {
  anchor <- sex_screen[direction == direction_name & interval_supported & magnitude_ge_0_75y, .N]
  replication <- sex_screen[direction == direction_name & validity_gate_pass, .N]
  if (anchor >= 1L && replication >= 2L) sex_route2 <- TRUE
}
sex_sign_reversal_across_cohorts <- uniqueN(sex_screen$direction) > 1L
sex_promotion_gate <- sex_route1 | sex_route2

# Versioned v2.4 cross-country differences for wealth within low education.
# Cohort resamples are independent; matching replicate labels provides a
# reproducible Monte Carlo draw from the difference of independent estimates.
pair_metrics <- metrics[
  module == "wealth_within_low_education" & stratum == "all" &
    metric %in% c(
      "gap_dfle", "contribution_onset", "contribution_recovery",
      "contribution_relapse", "contribution_post_disability_mortality",
      "contribution_pre_disability_mortality", "recovery_relapse_contribution"
    ),
  .(cohort, replicate, module, stratum, metric, estimate)
]
pair_point <- point_long[
  module == "wealth_within_low_education" & stratum == "all" &
    metric %in% unique(pair_metrics$metric),
  .(cohort, module, stratum, metric, point_estimate)
]
country_pairs <- combn(cohort_order, 2L, simplify = FALSE)
pairwise_parts <- lapply(country_pairs, function(pair_names) {
  cohort_a <- pair_names[[1L]]
  cohort_b <- pair_names[[2L]]
  paired <- merge(
    pair_metrics[cohort == cohort_a],
    pair_metrics[cohort == cohort_b],
    by = c("replicate", "module", "stratum", "metric"),
    suffixes = c("_a", "_b")
  )
  paired[, difference := estimate_a - estimate_b]
  out <- paired[, .(
    paired_valid_replicates = .N,
    bootstrap_mean_difference = mean(difference),
    bootstrap_sd_difference = sd(difference),
    ci_low = as.numeric(quantile(difference, 0.025, names = FALSE, type = 6)),
    bootstrap_median_difference = median(difference),
    ci_high = as.numeric(quantile(difference, 0.975, names = FALSE, type = 6))
  ), by = .(module, stratum, metric)]
  point_a <- pair_point[cohort == cohort_a, .(module, stratum, metric, point_a = point_estimate)]
  point_b <- pair_point[cohort == cohort_b, .(module, stratum, metric, point_b = point_estimate)]
  out <- merge(out, point_a, by = c("module", "stratum", "metric"), all.x = TRUE)
  out <- merge(out, point_b, by = c("module", "stratum", "metric"), all.x = TRUE)
  out[, `:=`(
    cohort_a = cohort_a,
    cohort_b = cohort_b,
    contrast = paste0(cohort_a, "_minus_", cohort_b),
    point_difference = point_a - point_b,
    minimum_paired_valid_required = minimum_valid,
    paired_validity_gate = paired_valid_replicates >= minimum_valid,
    interval_excludes_zero = ci_low > 0 | ci_high < 0
  )]
  out
})
pairwise_country <- rbindlist(pairwise_parts, fill = TRUE)
setcolorder(pairwise_country, c(
  "cohort_a", "cohort_b", "contrast", "module", "stratum", "metric",
  "point_a", "point_b", "point_difference", "paired_valid_replicates",
  "minimum_paired_valid_required", "paired_validity_gate",
  "bootstrap_mean_difference", "bootstrap_sd_difference", "ci_low",
  "bootstrap_median_difference", "ci_high", "interval_excludes_zero"
))
fwrite(pairwise_country, file.path(summary_dir, "low_education_wealth_country_heterogeneity.csv"))

mhas_heterogeneity <- pairwise_country[
  cohort_b == "MHAS" & metric == "recovery_relapse_contribution"
]
mhas_heterogeneity[, supported_country_contrast :=
  paired_validity_gate & abs(point_difference) >= 0.75 & interval_excludes_zero]
fwrite(mhas_heterogeneity, file.path(summary_dir, "low_education_wealth_mhas_contrast_screen.csv"))
mhas_heterogeneity_same_direction_two <- mhas_heterogeneity[
  supported_country_contrast == TRUE,
  .N, by = sign(point_difference)
][N >= 2L, .N] > 0L

decision <- data.table(
  phase = phase,
  target_replicates = target_replicates,
  minimum_valid_required = minimum_valid,
  all_cohort_module_validity_gates_pass = all(qc_summary$validity_gate_pass),
  low_education_wealth_same_direction_two_cohorts = lowedu_same_direction_two,
  low_education_wealth_promoted = all(qc_summary[module == "wealth_within_low_education"]$validity_gate_pass) &
    lowedu_same_direction_two,
  recovery_phase_same_phase_two_cohorts_secondary = phase_same_phase_two,
  recovery_phase_same_phase_three_cohorts_primary = phase_same_phase_three,
  recovery_phase_promoted_primary = all(qc_summary[module == "recovery_phase"]$validity_gate_pass) &
    phase_same_phase_three,
  sex_route1_two_interval_supported_0_50y = sex_route1,
  sex_route2_one_interval_supported_0_75y_plus_directional_replication = sex_route2,
  sex_sign_reversal_across_cohorts = sex_sign_reversal_across_cohorts,
  sex_heterogeneity_promoted = all(qc_summary[module == "wealth_by_sex"]$validity_gate_pass) & sex_promotion_gate,
  universal_sex_vulnerability_claim_blocked = sex_sign_reversal_across_cohorts,
  low_education_wealth_mhas_contrast_two_cohorts = mhas_heterogeneity_same_direction_two,
  low_education_wealth_cross_country_heterogeneity_promoted =
    mhas_heterogeneity_same_direction_two,
  decision_scope = decision_scope
)
decision[, exploratory_upgrade_gate :=
  all_cohort_module_validity_gates_pass &
  (low_education_wealth_promoted | recovery_phase_promoted_primary |
    sex_heterogeneity_promoted | low_education_wealth_cross_country_heterogeneity_promoted)]
fwrite(decision, file.path(summary_dir, "exploratory_upgrade_decision.csv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, paste0("sessionInfo_exploratory_bootstrap_summary_", phase, ".txt"))
)
cat("Exploratory bootstrap summary completed\n")
cat("phase=", phase, " target=", target_replicates, " minimum_valid=", minimum_valid, "\n", sep = "")
cat("all_validity_gates=", decision$all_cohort_module_validity_gates_pass, "\n", sep = "")
cat("low_education_wealth_promoted=", decision$low_education_wealth_promoted, "\n", sep = "")
cat("recovery_phase_primary=", decision$recovery_phase_promoted_primary,
  " secondary_two_cohort=", decision$recovery_phase_same_phase_two_cohorts_secondary, "\n", sep = "")
cat("sex_heterogeneity_promoted=", decision$sex_heterogeneity_promoted,
  " sign_reversal=", decision$sex_sign_reversal_across_cohorts, "\n", sep = "")
cat("exploratory_upgrade_gate=", decision$exploratory_upgrade_gate,
  " scope=", decision$decision_scope, "\n", sep = "")
