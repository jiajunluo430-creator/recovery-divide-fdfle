options(stringsAsFactors = FALSE, width = 220, warn = 1)

if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
sensitivity <- tolower(Sys.getenv("D4_REV_SENSITIVITY", unset = "primary"))
if (!sensitivity %in% c("primary", "adl_only", "at_least_two")) {
  stop("Unsupported D4_REV_SENSITIVITY for bootstrap summary: ", sensitivity)
}
analysis_module <- tolower(Sys.getenv(
  "D4_REV_ANALYSIS_MODULE",
  unset = "wealth_within_low_education"
))
allowed_modules <- c(
  "primary_education", "primary_wealth", "wealth_within_low_education"
)
if (!analysis_module %in% allowed_modules) {
  stop("Unsupported D4_REV_ANALYSIS_MODULE for bootstrap summary: ", analysis_module)
}
if (sensitivity != "primary" && analysis_module != "wealth_within_low_education") {
  stop("Threshold bootstrap summaries are defined only for wealth_within_low_education")
}
boot_tag <- if (sensitivity == "primary") "final" else paste0(sensitivity, "_final")
point_tag <- if (sensitivity == "primary") {
  "02_revision_point"
} else {
  paste0("02_revision_point_", sensitivity)
}
boot_prefix <- switch(
  analysis_module,
  primary_education = "03_primary_education_household_bootstrap_",
  primary_wealth = "03_primary_wealth_household_bootstrap_",
  wealth_within_low_education = "03_lowedu_household_bootstrap_"
)
inference_tag <- switch(
  analysis_module,
  primary_education = "04_primary_education_inference",
  primary_wealth = "04_primary_wealth_inference",
  wealth_within_low_education = if (sensitivity == "primary") {
    "04_lowedu_inference"
  } else {
    paste0("04_lowedu_inference_", sensitivity)
  }
)
boot_dir <- file.path(revision_root, "03_outputs", paste0(boot_prefix, boot_tag))
point_dir <- file.path(revision_root, "03_outputs", point_tag)
out_dir <- file.path(revision_root, "03_outputs", inference_tag)
log_dir <- file.path(revision_root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- c("CHARLS", "ELSA", "HRS", "MHAS")
default_minimum <- as.integer(Sys.getenv("D4_REV_MIN_BOOTSTRAP", unset = "500"))
default_target <- as.integer(Sys.getenv("D4_REV_TARGET_BOOTSTRAP", unset = "500"))
minimum_by_cohort <- setNames(vapply(cohorts, function(cohort_name) {
  as.integer(Sys.getenv(
    paste0("D4_REV_MIN_", cohort_name),
    unset = as.character(default_minimum)
  ))
}, integer(1)), cohorts)
target_by_cohort <- setNames(vapply(cohorts, function(cohort_name) {
  as.integer(Sys.getenv(
    paste0("D4_REV_TARGET_", cohort_name),
    unset = as.character(default_target)
  ))
}, integer(1)), cohorts)

read_one <- function(cohort_name, suffix) {
  path <- file.path(boot_dir, sprintf("bootstrap_%s_%s.csv", tolower(cohort_name), suffix))
  if (!file.exists(path)) stop("Missing bootstrap file: ", path)
  fread(path)
}

qc <- rbindlist(lapply(cohorts, read_one, suffix = "qc"), fill = TRUE)
metrics <- rbindlist(lapply(cohorts, read_one, suffix = "metrics"), fill = TRUE)
qc[, status := tolower(status)]
valid <- qc[status == "valid", .(cohort, replicate)]
replicate_qc <- qc[, .(
  attempted = .N,
  valid = sum(status == "valid"),
  failed = sum(status != "valid"),
  median_seconds = median(elapsed_seconds, na.rm = TRUE),
  p95_seconds = quantile(elapsed_seconds, 0.95, na.rm = TRUE)
), by = cohort]
replicate_qc[, minimum_required := minimum_by_cohort[cohort]]
replicate_qc[, target_for_inference := target_by_cohort[cohort]]
if (any(replicate_qc$valid < replicate_qc$minimum_required)) {
  stop(
    "Bootstrap incomplete: ",
    paste(replicate_qc[, paste0(cohort, "=", valid)], collapse = ", "),
    "; required=", paste(replicate_qc[, paste0(cohort, "=", minimum_required)], collapse = ", ")
  )
}
setorder(valid, cohort, replicate)
valid[, inference_order := seq_len(.N), by = cohort]
valid[, target_for_inference := target_by_cohort[cohort]]
valid_inference <- valid[inference_order <= target_for_inference, .(cohort, replicate)]
replicate_qc[, used_in_inference := pmin(valid, target_for_inference)]
if (any(replicate_qc$used_in_inference < replicate_qc$target_for_inference)) {
  stop("Fewer valid replicates than D4_REV_TARGET_BOOTSTRAP")
}
metrics <- metrics[valid_inference, on = .(cohort, replicate), nomatch = 0]
setorder(metrics, cohort, replicate, metric)

qfun <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE, type = 6))
bootstrap_ci <- metrics[, .(
  valid_replicates = .N,
  bootstrap_mean = mean(estimate),
  bootstrap_sd = sd(estimate),
  ci_low = qfun(estimate, 0.025),
  median = qfun(estimate, 0.5),
  ci_high = qfun(estimate, 0.975)
), by = .(cohort, metric)]

# Convert the independently fitted point-model outputs to the metric names used
# in the household bootstrap. These are the reported estimates; bootstrap means
# are retained only as bias and Monte Carlo diagnostics.
module_name <- analysis_module
point_parts <- list()

init <- fread(file.path(point_dir, "initial_difficulty_probability_age60.csv"))[
  module == module_name & ses %chin% c("high", "low")
]
point_parts[[length(point_parts) + 1L]] <- init[, .(
  cohort,
  metric = paste0("initial_difficulty_", ses),
  point_estimate = difficulty_probability_age60
)]

life <- fread(file.path(point_dir, "life_expectancy_point_estimates.csv"))[
  module == module_name & ses %chin% c("high", "low")
]
metric_suffix <- c(
  total_life_expectancy = "tle",
  fdfle = "fdfle",
  difficulty_years = "difficulty_years",
  fdfle_percent_of_life = "fdfle_percent"
)
life <- life[metric %chin% names(metric_suffix)]
life[, prefix := fifelse(estimand == "population_initialised", "population", "conditional")]
point_parts[[length(point_parts) + 1L]] <- life[, .(
  cohort,
  metric = paste(prefix, ses, unname(metric_suffix[metric]), sep = "_"),
  point_estimate = estimate
)]

gaps <- fread(file.path(point_dir, "low_high_fdfle_gaps.csv"))[module == module_name]
point_parts[[length(point_parts) + 1L]] <- gaps[, .(
  cohort,
  metric = fifelse(estimand == "population_initialised", "population_fdfle_gap", "conditional_fdfle_gap"),
  point_estimate = low_minus_high_gap
)]

shapley <- fread(file.path(point_dir, "shapley_decomposition_point.csv"))[module == module_name]
shapley[, prefix := fifelse(estimand == "population_initialised", "population", "conditional")]
point_parts[[length(point_parts) + 1L]] <- shapley[, .(
  cohort,
  metric = paste0(prefix, "_contribution_", block),
  point_estimate = contribution_years
)]
point_parts[[length(point_parts) + 1L]] <- shapley[
  block %chin% c("recovery", "relapse"),
  .(point_estimate = sum(contribution_years)),
  by = .(cohort, estimand)
][, .(
  cohort,
  metric = fifelse(estimand == "population_initialised", "population_recovery_relapse", "conditional_recovery_relapse"),
  point_estimate
)]

coefs <- fread(file.path(point_dir, "transition_model_coefficients.csv"))[
  module == module_name & term == "seslow"
]
point_parts[[length(point_parts) + 1L]] <- coefs[, .(
  cohort,
  metric = paste0("hazard_ratio_", process),
  point_estimate = hazard_ratio
)]

point_metrics <- unique(rbindlist(point_parts, fill = TRUE), by = c("cohort", "metric"))
bootstrap_ci <- point_metrics[bootstrap_ci, on = .(cohort, metric)]
bootstrap_ci[, bootstrap_bias := bootstrap_mean - point_estimate]
bootstrap_ci[, null_value := fifelse(grepl("^hazard_ratio_", metric), 1, 0)]
sign_counts <- metrics[, .(
  lower_or_equal_null = sum(estimate <= fifelse(grepl("^hazard_ratio_", metric), 1, 0)),
  higher_or_equal_null = sum(estimate >= fifelse(grepl("^hazard_ratio_", metric), 1, 0)),
  n = .N
), by = .(cohort, metric)]
bootstrap_ci <- sign_counts[bootstrap_ci, on = .(cohort, metric)]
bootstrap_ci[, bootstrap_p_two_sided := pmin(1, 2 * pmin(
  (lower_or_equal_null + 1) / (n + 1),
  (higher_or_equal_null + 1) / (n + 1)
))]

# Percent contribution is interpretable only when the total population FDFLE
# gap is bounded away from zero. Ratios are calculated within replicate.
ratio_blocks <- c(
  "population_recovery_relapse",
  "population_contribution_initial_state",
  "population_contribution_onset",
  "population_contribution_recovery",
  "population_contribution_relapse",
  "population_contribution_post_difficulty_mortality",
  "population_contribution_pre_difficulty_mortality"
)
ratio_wide <- dcast(
  metrics[metric %chin% c("population_fdfle_gap", ratio_blocks)],
  cohort + replicate ~ metric,
  value.var = "estimate"
)
ratio_parts <- list()
for (block_name in ratio_blocks) {
  for (cohort_name in cohorts) {
    gap_ci <- bootstrap_ci[cohort == cohort_name & metric == "population_fdfle_gap"]
    reportable <- nrow(gap_ci) == 1L && (gap_ci$ci_high < 0 || gap_ci$ci_low > 0)
    point_num <- point_metrics[cohort == cohort_name & metric == block_name, point_estimate]
    point_den <- point_metrics[cohort == cohort_name & metric == "population_fdfle_gap", point_estimate]
    ratios <- 100 * ratio_wide[cohort == cohort_name][[block_name]] /
      ratio_wide[cohort == cohort_name][["population_fdfle_gap"]]
    ratio_parts[[length(ratio_parts) + 1L]] <- data.table(
      cohort = cohort_name,
      component = block_name,
      denominator_gap_ci_excludes_zero = reportable,
      point_percent = if (reportable) 100 * point_num / point_den else NA_real_,
      ci_low = if (reportable) qfun(ratios, 0.025) else NA_real_,
      ci_high = if (reportable) qfun(ratios, 0.975) else NA_real_
    )
  }
}
percent_contributions <- rbindlist(ratio_parts)

heterogeneity_one <- function(metric_name) {
  s <- bootstrap_ci[metric == metric_name & is.finite(bootstrap_sd) & bootstrap_sd > 0]
  if (nrow(s) < 2L) return(NULL)
  w <- 1 / s$bootstrap_sd^2
  pooled <- sum(w * s$point_estimate) / sum(w)
  q <- sum(w * (s$point_estimate - pooled)^2)
  df <- nrow(s) - 1L
  data.table(
    metric = metric_name,
    cohorts = nrow(s),
    inverse_variance_pooled = pooled,
    cochran_q = q,
    df = df,
    p_heterogeneity = pchisq(q, df = df, lower.tail = FALSE),
    i2_percent = if (q > 0) max(0, 100 * (q - df) / q) else 0
  )
}
heterogeneity_metrics <- c("population_fdfle_gap", "population_recovery_relapse")
heterogeneity <- rbindlist(lapply(heterogeneity_metrics, heterogeneity_one), fill = TRUE)

pairwise_parts <- list()
for (metric_name in heterogeneity_metrics) {
  s <- bootstrap_ci[metric == metric_name]
  pairs <- combn(cohorts[cohorts %chin% s$cohort], 2, simplify = FALSE)
  for (pair in pairs) {
    a <- s[cohort == pair[[1L]]]
    b <- s[cohort == pair[[2L]]]
    contrast <- a$point_estimate - b$point_estimate
    se <- sqrt(a$bootstrap_sd^2 + b$bootstrap_sd^2)
    z <- contrast / se
    pairwise_parts[[length(pairwise_parts) + 1L]] <- data.table(
      metric = metric_name,
      cohort_1 = pair[[1L]],
      cohort_2 = pair[[2L]],
      contrast_1_minus_2 = contrast,
      se_independent_cohorts = se,
      ci_low = contrast - 1.96 * se,
      ci_high = contrast + 1.96 * se,
      p_two_sided = 2 * pnorm(abs(z), lower.tail = FALSE)
    )
  }
}
pairwise <- rbindlist(pairwise_parts)

# Pre-specified Monte Carlo stop rule: compare the first half of replicates with
# the complete run for the key year-scale endpoints. Increase to 1000 when a CI
# endpoint moves by >0.10 years or the CI half-width changes by >20%.
key_year_metrics <- c(
  "population_fdfle_gap", "conditional_fdfle_gap", "population_recovery_relapse",
  "population_contribution_initial_state", "population_contribution_onset",
  "population_contribution_recovery", "population_contribution_relapse",
  "population_contribution_post_difficulty_mortality",
  "population_contribution_pre_difficulty_mortality"
)
stability <- metrics[metric %chin% key_year_metrics, {
  reps <- sort(unique(replicate))
  half_reps <- reps[seq_len(floor(length(reps) / 2))]
  full_low <- qfun(estimate, 0.025)
  full_high <- qfun(estimate, 0.975)
  half_low <- qfun(estimate[replicate %in% half_reps], 0.025)
  half_high <- qfun(estimate[replicate %in% half_reps], 0.975)
  full_hw <- (full_high - full_low) / 2
  half_hw <- (half_high - half_low) / 2
  endpoint_drift <- max(abs(c(half_low - full_low, half_high - full_high)))
  relative_halfwidth_change <- if (full_hw > 0) abs(half_hw - full_hw) / full_hw else NA_real_
  .(
    half_replicates = length(half_reps),
    full_replicates = length(reps),
    half_ci_low = half_low,
    half_ci_high = half_high,
    full_ci_low = full_low,
    full_ci_high = full_high,
    max_endpoint_drift_years = endpoint_drift,
    relative_halfwidth_change = relative_halfwidth_change,
    requires_1000 = endpoint_drift > 0.10 || relative_halfwidth_change > 0.20
  )
}, by = .(cohort, metric)]

setorder(bootstrap_ci, cohort, metric)
setorder(percent_contributions, cohort, component)
setorder(stability, cohort, metric)
fwrite(replicate_qc, file.path(out_dir, "bootstrap_replicate_qc.csv"))
fwrite(bootstrap_ci, file.path(out_dir, "bootstrap_percentile_ci_summary.csv"))
fwrite(percent_contributions, file.path(out_dir, "population_percent_contributions.csv"))
fwrite(heterogeneity, file.path(out_dir, "cross_cohort_heterogeneity.csv"))
fwrite(pairwise, file.path(out_dir, "cross_cohort_pairwise_contrasts.csv"))
fwrite(stability, file.path(out_dir, "bootstrap_monte_carlo_stability.csv"))

key_ci <- bootstrap_ci[metric %chin% c("population_fdfle_gap", "population_recovery_relapse")]
report <- c(
  "# Household-cluster bootstrap inference summary",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Analysis module: %s", analysis_module),
  sprintf("Functional-difficulty definition: %s", sensitivity),
  sprintf("Valid replicates: %s", paste(replicate_qc[, paste0(cohort, "=", valid)], collapse = "; ")),
  "",
  "## Key percentile intervals",
  "",
  paste(capture.output(print(key_ci[, .(cohort, metric, point_estimate, ci_low, ci_high, bootstrap_p_two_sided)])), collapse = "\n"),
  "",
  "## Cross-cohort heterogeneity",
  "",
  paste(capture.output(print(heterogeneity)), collapse = "\n"),
  "",
  "## Monte Carlo stop rule",
  "",
  sprintf("Endpoints triggering extension to 1000 replicates: %d", sum(stability$requires_1000)),
  paste(capture.output(print(stability[requires_1000 == TRUE])), collapse = "\n")
)
writeLines(report, file.path(out_dir, "bootstrap_inference_summary.md"), useBytes = TRUE)
writeLines(
  capture.output(sessionInfo()),
  file.path(
    log_dir,
    paste0("sessionInfo_42_bootstrap_summary_", analysis_module, "_", sensitivity, ".txt")
  )
)

cat(paste(report, collapse = "\n"), "\n")
