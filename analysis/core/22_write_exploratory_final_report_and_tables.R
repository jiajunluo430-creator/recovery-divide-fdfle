#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
v1_dir <- file.path(root, "03_outputs", "06_bootstrap_summary_final")
v2_dir <- file.path(root, "03_outputs", "13_exploratory_bootstrap_summary_final")
report_dir <- file.path(root, "05_reports")
table_dir <- file.path(root, "03_outputs", "17_exploratory_submission_tables_final")
log_dir <- file.path(root, "06_logs")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

required_v2 <- c(
  "bootstrap_validity_qc.csv",
  "bootstrap_percentile_intervals.csv",
  "low_education_wealth_promotion_by_cohort.csv",
  "low_education_wealth_mhas_contrast_screen.csv",
  "recovery_phase_promotion_by_cohort.csv",
  "sex_heterogeneity_promotion_by_cohort.csv",
  "exploratory_upgrade_decision.csv"
)
missing_v2 <- required_v2[!file.exists(file.path(v2_dir, required_v2))]
if (length(missing_v2)) stop("Missing final v2 summary files: ", paste(missing_v2, collapse = "; "))

cohort_order <- c("CHARLS", "ELSA", "HRS", "MHAS")
fmt <- function(x, digits = 2L) ifelse(is.finite(x), sprintf(paste0("%.", digits, "f"), x), "NA")
yn <- function(x) ifelse(isTRUE(x), "PASS", "NO")
set_cohort_order <- function(z, extra = character()) {
  z[, .cohort_order := match(cohort, cohort_order)]
  setorderv(z, c(".cohort_order", extra))
  z[, .cohort_order := NULL]
  invisible(z)
}

validity <- fread(file.path(v2_dir, "bootstrap_validity_qc.csv"))
intervals <- fread(file.path(v2_dir, "bootstrap_percentile_intervals.csv"))
lowedu <- fread(file.path(v2_dir, "low_education_wealth_promotion_by_cohort.csv"))
mhas_heterogeneity <- fread(file.path(v2_dir, "low_education_wealth_mhas_contrast_screen.csv"))
phase <- fread(file.path(v2_dir, "recovery_phase_promotion_by_cohort.csv"))
sex <- fread(file.path(v2_dir, "sex_heterogeneity_promotion_by_cohort.csv"))
decision <- fread(file.path(v2_dir, "exploratory_upgrade_decision.csv"))

stopifnot(nrow(decision) == 1L)
stopifnot(decision$decision_scope == "binding_final")
stopifnot(decision$target_replicates == 500L)
stopifnot(all(validity$target_replicates == 500L))
stopifnot(all(validity$minimum_valid_required == 450L))

set_cohort_order(validity, "module")
set_cohort_order(lowedu)
set_cohort_order(phase, "stratum")
set_cohort_order(sex)

# Machine-readable final v2 manuscript tables.
fwrite(validity, file.path(table_dir, "Table_S_v2_bootstrap_validity.csv"))
fwrite(lowedu, file.path(table_dir, "Table_3_low_education_wealth_decomposition.csv"))
fwrite(mhas_heterogeneity, file.path(table_dir, "Table_S_low_education_wealth_country_heterogeneity.csv"))
fwrite(phase, file.path(table_dir, "Table_S_recovery_durability.csv"))
fwrite(sex, file.path(table_dir, "Table_S_sex_heterogeneity.csv"))
fwrite(decision, file.path(table_dir, "Table_S_v2_binding_decision.csv"))

v1_promotion <- fread(file.path(v1_dir, "mechanism_promotion_by_cohort.csv"))
v1_intervals <- fread(file.path(v1_dir, "bootstrap_percentile_intervals.csv"))
v1_gap <- v1_intervals[metric == "gap_dfle", .(
  cohort, exposure,
  v1_dfle_gap = point_estimate,
  v1_dfle_gap_ci_low = ci_low,
  v1_dfle_gap_ci_high = ci_high
)]
v1_rr <- v1_promotion[, .(
  cohort, exposure,
  v1_recovery_relapse = recovery_relapse_contribution_years,
  v1_recovery_relapse_ci_low = rr_ci_low,
  v1_recovery_relapse_ci_high = rr_ci_high,
  v1_recovery_relapse_percent = recovery_relapse_percent,
  v1_mechanism_trigger = mechanism_trigger_cohort
)]
v1_integrated <- merge(v1_gap, v1_rr, by = c("cohort", "exposure"), all = TRUE)
v1_integrated <- dcast(
  v1_integrated,
  cohort ~ exposure,
  value.var = setdiff(names(v1_integrated), c("cohort", "exposure"))
)
v2_integrated <- lowedu[, .(
  cohort,
  v2_lowedu_wealth_dfle_gap = dfle_gap,
  v2_lowedu_wealth_dfle_gap_ci_low = ci_low_gap_dfle,
  v2_lowedu_wealth_dfle_gap_ci_high = ci_high_gap_dfle,
  v2_lowedu_wealth_recovery_relapse = recovery_relapse_contribution_years,
  v2_lowedu_wealth_recovery_relapse_ci_low = ci_low_recovery_relapse_contribution,
  v2_lowedu_wealth_recovery_relapse_ci_high = ci_high_recovery_relapse_contribution,
  v2_lowedu_wealth_recovery_relapse_percent = recovery_relapse_percent,
  v2_lowedu_wealth_promotion = cohort_promotion_gate
)]
integrated <- merge(v1_integrated, v2_integrated, by = "cohort", all = TRUE)
set_cohort_order(integrated)
fwrite(integrated, file.path(table_dir, "Table_2_integrated_primary_and_loweducation_results.csv"))

report_lines <- c(
  "# v2 final 500-replicate bootstrap results",
  "",
  "Status: binding final inference under exploratory addenda v2.1-v2.4. Preview intervals are superseded by this report.",
  "",
  "## Computational validity",
  "",
  paste0(
    "All 12 cohort-module validity gates: `", yn(decision$all_cohort_module_validity_gates_pass),
    "`. Valid replicates ranged from ", min(validity$valid_replicates), " to ",
    max(validity$valid_replicates), " of 500; the frozen minimum was 450."
  ),
  "",
  "| Cohort | Module | QC replicates | Valid | Failed | Gate |",
  "|---|---|---:|---:|---:|---|"
)
for (i in seq_len(nrow(validity))) {
  z <- validity[i]
  report_lines <- c(report_lines, sprintf(
    "| %s | %s | %d | %d | %d | %s |",
    z$cohort, z$module, z$qc_replicates, z$valid_replicates,
    z$failed_replicates, yn(z$validity_gate_pass)
  ))
}

report_lines <- c(
  report_lines,
  "",
  "## Wealth-related recovery and relapse within low education",
  "",
  "| Cohort | Low-minus-high wealth DFLE, y (95% percentile interval) | Recovery+relapse contribution, y (95% percentile interval) | Share, % | Promotion |",
  "|---|---:|---:|---:|---|"
)
for (i in seq_len(nrow(lowedu))) {
  z <- lowedu[i]
  report_lines <- c(report_lines, sprintf(
    "| %s | %s (%s to %s) | %s (%s to %s) | %s | %s |",
    z$cohort,
    fmt(z$dfle_gap), fmt(z$ci_low_gap_dfle), fmt(z$ci_high_gap_dfle),
    fmt(z$recovery_relapse_contribution_years),
    fmt(z$ci_low_recovery_relapse_contribution),
    fmt(z$ci_high_recovery_relapse_contribution),
    fmt(z$recovery_relapse_percent, 1L), yn(z$cohort_promotion_gate)
  ))
}
report_lines <- c(
  report_lines,
  "",
  paste0(
    "Binding C2 promotion across at least two cohorts: `",
    yn(decision$low_education_wealth_promoted), "`."
  ),
  "",
  "## Formal country contrasts against MHAS",
  "",
  "| Contrast | Difference in recovery+relapse contribution, y | 95% percentile interval | Paired valid replicates | Promotion |",
  "|---|---:|---:|---:|---|"
)
for (i in seq_len(nrow(mhas_heterogeneity))) {
  z <- mhas_heterogeneity[i]
  report_lines <- c(report_lines, sprintf(
    "| %s | %s | %s to %s | %d | %s |",
    z$contrast, fmt(z$point_difference), fmt(z$ci_low), fmt(z$ci_high),
    z$paired_valid_replicates, yn(z$supported_country_contrast)
  ))
}
report_lines <- c(
  report_lines,
  "",
  paste0(
    "At least two same-direction country-versus-MHAS contrasts: `",
    yn(decision$low_education_wealth_cross_country_heterogeneity_promoted), "`."
  ),
  "",
  "## Recovery durability",
  "",
  "| Cohort | Recovery phase | Low/high wealth relapse HR | 95% percentile interval | Cohort-phase gate |",
  "|---|---|---:|---:|---|"
)
for (i in seq_len(nrow(phase[stratum != "sustained_vs_early_modification"]))) {
  z <- phase[stratum != "sustained_vs_early_modification"][i]
  report_lines <- c(report_lines, sprintf(
    "| %s | %s | %s | %s to %s | %s |",
    z$cohort, z$stratum, fmt(z$hazard_ratio),
    fmt(z$bootstrap_hr_ci_low), fmt(z$bootstrap_hr_ci_high),
    yn(z$cohort_phase_gate)
  ))
}
report_lines <- c(
  report_lines,
  "",
  paste0(
    "Three-cohort primary durability promotion: `", yn(decision$recovery_phase_promoted_primary),
    "`; exactly-two-cohort secondary replication: `",
    yn(decision$recovery_phase_same_phase_two_cohorts_secondary), "`."
  ),
  "",
  "## Sex heterogeneity",
  "",
  "| Cohort | Women-minus-men recovery+relapse contribution, y | 95% percentile interval | Interval support | Direction |",
  "|---|---:|---:|---|---|"
)
for (i in seq_len(nrow(sex))) {
  z <- sex[i]
  report_lines <- c(report_lines, sprintf(
    "| %s | %s | %s to %s | %s | %s |",
    z$cohort, fmt(z$female_minus_male_recovery_relapse),
    fmt(z$bootstrap_ci_low), fmt(z$bootstrap_ci_high),
    ifelse(z$interval_supported, "YES", "NO"), z$direction
  ))
}
report_lines <- c(
  report_lines,
  "",
  paste0(
    "Sign reversal across cohorts: `", ifelse(decision$sex_sign_reversal_across_cohorts, "YES", "NO"),
    "`; universal sex-vulnerability claim blocked: `",
    ifelse(decision$universal_sex_vulnerability_claim_blocked, "YES", "NO"), "`."
  ),
  "",
  "## Binding v2 decision",
  "",
  paste0("- Low-education wealth module: `", yn(decision$low_education_wealth_promoted), "`."),
  paste0("- Formal non-MHAS versus MHAS recovery/relapse heterogeneity: `", yn(decision$low_education_wealth_cross_country_heterogeneity_promoted), "`."),
  paste0("- Recovery durability as a four-country primary mechanism: `", yn(decision$recovery_phase_promoted_primary), "`."),
  paste0("- Recovery durability as a two-context secondary result: `", yn(decision$recovery_phase_same_phase_two_cohorts_secondary), "`."),
  paste0("- Sex heterogeneity module: `", yn(decision$sex_heterogeneity_promoted), "`, with a universal sex claim blocked."),
  paste0("- Overall exploratory upgrade gate: `", yn(decision$exploratory_upgrade_gate), "`."),
  "",
  "Negative and stopped modules remain binding: four-cell education-by-wealth discordance, age-varying SES, and recovery exhaustion after relapse were not rescued by recoding, cohort deletion, or model tuning.",
  "",
  "## Interpretation boundary",
  "",
  "These estimates are descriptive transition-process contributions. Shapley replacement is accounting of a modelled contrast, not an intervention effect or causal mediation. The estimates do not identify a causal effect of wealth, education, health care, rehabilitation, or national policy. MHAS is reported wherever the replicated adverse pattern in CHARLS, ELSA, and HRS is described."
)
writeLines(report_lines, file.path(report_dir, "15_exploratory_bootstrap_final.md"))

writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, "sessionInfo_exploratory_final_report_and_tables.txt")
)
cat("Wrote binding v2 final report and manuscript tables.\n")
