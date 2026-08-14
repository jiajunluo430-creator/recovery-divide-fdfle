#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
report_dir <- file.path(root, "05_reports")
log_dir <- file.path(root, "06_logs")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cohort_order <- c("CHARLS", "ELSA", "HRS", "MHAS")
fmt <- function(x, digits = 2L) ifelse(is.finite(x), sprintf(paste0("%.", digits, "f"), x), "NA")
set_cohort_order <- function(z, extra = character()) {
  z[, .cohort_order := match(cohort, cohort_order)]
  setorderv(z, c(".cohort_order", extra))
  z[, .cohort_order := NULL]
  invisible(z)
}

point_dir <- file.path(root, "03_outputs", "11_exploratory_upgrade_point")
support_dir <- file.path(root, "03_outputs", "10_exploratory_upgrade_support")
preview_dir <- file.path(root, "03_outputs", "13_exploratory_bootstrap_summary_preview")
sensitivity_dir <- file.path(root, "03_outputs", "16_exploratory_sensitivity_summary")

low_point <- fread(file.path(point_dir, "low_education_wealth_compensation_point_screen.csv"))
phase_point <- fread(file.path(point_dir, "recovery_phase_model_contrasts.csv"))
sex_point <- fread(file.path(point_dir, "sex_difference_point_screen.csv"))
joint_people <- fread(file.path(support_dir, "joint_ses_person_support.csv"))
exhaustion_gate <- fread(file.path(root, "03_outputs", "14_recovery_exhaustion_support", "recovery_exhaustion_support_gate_by_cohort.csv"))

set_cohort_order(low_point)
set_cohort_order(sex_point)

point_lines <- c(
  "# v2 exploratory support and point-model results",
  "",
  "Status: point estimates and support decisions; uncertainty claims require the final person bootstrap.",
  "",
  "## Support-driven module decisions",
  "",
  paste0(
    "The original four-cell education-by-wealth discordance module was stopped before outcome fitting. ",
    "The CHARLS high-education/low-wealth cell contained ",
    joint_people[cohort == "CHARLS" & joint_ses4 == "highedu_lowwealth", persons],
    " people, so the v2 four-cell support contract could not be met."
  ),
  "",
  "The replacement C2 contrast—low versus high wealth within the frozen low-education stratum—had at least 100 recovery and 50 relapse events in both decisive cells in all four cohorts. It is interpreted as residual wealth stratification within low education, not causal education compensation.",
  "",
  "## Wealth gradient within low education",
  "",
  "| Cohort | Low-minus-high DFLE, y | Recovery+relapse contribution, y | Share of DFLE gap, % | Point gate |",
  "|---|---:|---:|---:|---|"
)
for (i in seq_len(nrow(low_point))) {
  z <- low_point[i]
  point_lines <- c(point_lines, sprintf(
    "| %s | %s | %s | %s | %s |",
    z$cohort, fmt(z$dfle_gap), fmt(z$recovery_relapse_contribution_years),
    fmt(z$recovery_relapse_percent, 1L), ifelse(z$point_threshold_met, "PASS", "NO")
  ))
}

point_lines <- c(
  point_lines,
  "",
  "The point pattern is replicated in CHARLS, ELSA, and HRS. MHAS is null/opposing and is retained as substantive heterogeneity rather than tuned away.",
  "",
  "## Recovery durability",
  "",
  "| Cohort | Early recovery: low/high relapse HR | Sustained recovery: low/high relapse HR | Sustained-vs-early modification |",
  "|---|---:|---:|---:|"
)
phase_wide <- dcast(phase_point, cohort ~ contrast, value.var = "hazard_ratio")
set_cohort_order(phase_wide)
for (i in seq_len(nrow(phase_wide))) {
  z <- phase_wide[i]
  point_lines <- c(point_lines, sprintf(
    "| %s | %s | %s | %s |",
    z$cohort,
    fmt(z$low_vs_high_early_recovery),
    fmt(z$low_vs_high_sustained_recovery),
    fmt(z$sustained_vs_early_modification)
  ))
}

point_lines <- c(
  point_lines,
  "",
  "ELSA and HRS show low-wealth excess relapse both immediately after observed recovery and after recovery persisted for another observed wave. Under v2.1 this is a two-context secondary result, not a four-country primary durability mechanism.",
  "",
  "## Sex heterogeneity",
  "",
  "| Cohort | Female-minus-male recovery+relapse contribution, y | Point direction |",
  "|---|---:|---|"
)
for (i in seq_len(nrow(sex_point))) {
  z <- sex_point[i]
  direction <- if (z$female_minus_male_recovery_relapse < 0) "more adverse in women" else "more adverse in men"
  point_lines <- c(point_lines, sprintf(
    "| %s | %s | %s |", z$cohort, fmt(z$female_minus_male_recovery_relapse), direction
  ))
}

point_lines <- c(
  point_lines,
  "",
  "The sign reverses across countries, which blocks a universal female or male vulnerability claim.",
  "",
  "## Binding stops",
  "",
  "- Age-varying wealth effects: `STOP_SECOND_IDENTICAL_FAILURE`; the unchanged ELSA onset model returned non-finite coefficients twice.",
  paste0(
    "- Recovery exhaustion after relapse: `STOP_TWO_OR_MORE_COHORT_FAILURES`; minimum decisive-cell recoveries were ",
    paste(exhaustion_gate[, paste0(cohort, " ", minimum_recovery_events)], collapse = ", "),
    ". No effect model was fitted."
  ),
  "- Four-cell education-by-wealth discordance: `STOP_STRUCTURAL_SUPPORT`; it is not revived by the C2 low-education contrast.",
  "",
  "## Point-stage decision",
  "",
  "`PILOT_POSITIVE` for C2 and recovery durability; `PILOT_HETEROGENEITY` for sex; age, four-cell discordance, and recurrent-recovery exhaustion remain stopped."
)
writeLines(point_lines, file.path(report_dir, "11_exploratory_upgrade_point_results.md"))

validity <- fread(file.path(preview_dir, "bootstrap_validity_qc.csv"))
low_preview <- fread(file.path(preview_dir, "low_education_wealth_promotion_by_cohort.csv"))
phase_preview <- fread(file.path(preview_dir, "recovery_phase_promotion_by_cohort.csv"))
sex_preview <- fread(file.path(preview_dir, "sex_heterogeneity_promotion_by_cohort.csv"))
decision <- fread(file.path(preview_dir, "exploratory_upgrade_decision.csv"))
set_cohort_order(low_preview)
set_cohort_order(sex_preview)

preview_lines <- c(
  "# Exploratory bootstrap preview",
  "",
  "Status: 100-person-bootstrap computational preview; final claims require 500 replicates and at least 450 valid replicates per cohort-module.",
  "",
  "## Computational validity",
  "",
  paste0(
    "All 12 cohort-module validity gates passed. Valid life-table replicates ranged from ",
    min(validity$valid_replicates), " to ", max(validity$valid_replicates),
    "; all recovery-phase models had 100 valid replicates. Failed resamples remain recorded and were not repaired by model tuning."
  ),
  "",
  "## Wealth gradient within low education",
  "",
  "| Cohort | Recovery+relapse, y | Preview 95% percentile interval | Low-minus-high DFLE, y | Preview 95% percentile interval | Promotion |",
  "|---|---:|---:|---:|---:|---|"
)
for (i in seq_len(nrow(low_preview))) {
  z <- low_preview[i]
  preview_lines <- c(preview_lines, sprintf(
    "| %s | %s | %s to %s | %s | %s to %s | %s |",
    z$cohort,
    fmt(z$recovery_relapse_contribution_years),
    fmt(z$ci_low_recovery_relapse_contribution), fmt(z$ci_high_recovery_relapse_contribution),
    fmt(z$dfle_gap), fmt(z$ci_low_gap_dfle), fmt(z$ci_high_gap_dfle),
    ifelse(z$cohort_promotion_gate, "PASS", "NO")
  ))
}

phase_pass <- phase_preview[cohort_phase_gate == TRUE]
preview_lines <- c(
  preview_lines,
  "",
  paste0(
    "Recovery durability passed in exactly two cohorts and both phases: ",
    paste(unique(phase_pass$cohort), collapse = " and "),
    ". The three-cohort primary durability gate was not met."
  ),
  "",
  "## Sex result",
  "",
  "| Cohort | Female-minus-male contribution, y | Preview 95% percentile interval | Interval support |",
  "|---|---:|---:|---|"
)
for (i in seq_len(nrow(sex_preview))) {
  z <- sex_preview[i]
  preview_lines <- c(preview_lines, sprintf(
    "| %s | %s | %s to %s | %s |",
    z$cohort, fmt(z$female_minus_male_recovery_relapse),
    fmt(z$bootstrap_ci_low), fmt(z$bootstrap_ci_high),
    ifelse(z$interval_supported, "YES", "NO")
  ))
}

preview_lines <- c(
  preview_lines,
  "",
  "Only ELSA had interval support for the sex contrast. HRS supplied same-direction point replication, but the cross-country sign reversal blocks a universal sex claim.",
  "",
  "## Preview decision",
  "",
  paste0("- C2 low-education wealth promotion: `", ifelse(decision$low_education_wealth_promoted, "PASS", "NO"), "`."),
  paste0("- Recovery durability primary promotion: `", ifelse(decision$recovery_phase_promoted_primary, "PASS", "NO"), "`; two-cohort secondary replication: `", ifelse(decision$recovery_phase_same_phase_two_cohorts_secondary, "PASS", "NO"), "`."),
  paste0("- Exploratory upgrade gate: `", ifelse(decision$exploratory_upgrade_gate, "PASS_PREVIEW", "NO"), "`."),
  "",
  "These decisions are computational previews, not binding final inference."
)
writeLines(preview_lines, file.path(report_dir, "13_exploratory_bootstrap_preview.md"))

sensitivity <- fread(file.path(sensitivity_dir, "low_education_wealth_sensitivity_point_estimates.csv"))
robustness <- fread(file.path(sensitivity_dir, "low_education_wealth_sensitivity_robustness.csv"))
sensitivity_decision <- fread(file.path(sensitivity_dir, "exploratory_sensitivity_decision.csv"))
set_cohort_order(sensitivity, "sensitivity")
set_cohort_order(robustness)

sensitivity_lines <- c(
  "# v2 exploratory sensitivity results",
  "",
  "The C2 wealth-within-low-education model was repeated without survey weights, with only 1–3-year intervals, and after excluding explicitly identified proxy interviews. Functional items, state history, SES cuts, and estimands were unchanged.",
  "",
  "| Cohort | Implementation | Low-minus-high DFLE, y | Recovery+relapse, y | Share, % |",
  "|---|---|---:|---:|---:|"
)
for (i in seq_len(nrow(sensitivity))) {
  z <- sensitivity[i]
  sensitivity_lines <- c(sensitivity_lines, sprintf(
    "| %s | %s | %s | %s | %s |",
    z$cohort, z$sensitivity, fmt(z$dfle_gap),
    fmt(z$recovery_relapse_contribution_years), fmt(z$recovery_relapse_percent, 1L)
  ))
}

sensitivity_lines <- c(
  sensitivity_lines,
  "",
  "## Robustness decision",
  ""
)
for (i in seq_len(nrow(robustness))) {
  z <- robustness[i]
  sensitivity_lines <- c(sensitivity_lines, sprintf(
    "- %s: directionally robust=%s; recovery+relapse magnitude range %s–%s years; all four implementations at least 0.50 adverse years=%s.",
    z$cohort,
    ifelse(z$cohort_robustness_gate, "yes", "no"),
    fmt(z$minimum_absolute_recovery_relapse_years), fmt(z$maximum_absolute_recovery_relapse_years),
    ifelse(z$all_implementations_meet_adverse_0_50y, "yes", "no")
  ))
}
sensitivity_lines <- c(
  sensitivity_lines,
  "",
  paste0(
    "All sensitivity models converged: ",
    ifelse(sensitivity_decision$all_sensitivity_models_converged, "yes", "no"),
    ". At least two cohorts retained the adverse C2 mechanism direction: ",
    ifelse(sensitivity_decision$low_education_wealth_directionally_robust_in_at_least_two_cohorts, "yes", "no"),
    "."
  ),
  "",
  "MHAS remains a negative/heterogeneous result: its recovery-plus-relapse direction changes across implementations, while its total low-wealth DFLE gap remains adverse."
)
writeLines(sensitivity_lines, file.path(report_dir, "14_exploratory_sensitivity_results.md"))

writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, "sessionInfo_exploratory_interim_report_writer.txt")
)
cat("Wrote exploratory interim reports 11, 13, and 14.\n")
