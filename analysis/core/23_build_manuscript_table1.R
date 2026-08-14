#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
derived_dir <- file.path(root, "02_derived")
gate_dir <- file.path(root, "03_outputs", "01_gate0")
risk_dir <- file.path(root, "03_outputs", "03_formal_risksets")
out_dir <- file.path(root, "03_outputs", "17_exploratory_submission_tables_final")
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cohort_order <- c("CHARLS", "ELSA", "HRS", "MHAS")
panel <- as.data.table(readRDS(file.path(derived_dir, "individual_wave_panel_gate0.rds")))
person_ses <- as.data.table(readRDS(file.path(derived_dir, "person_fixed_ses_internal.rds")))[, .(
  cohort, person_id, education3_fixed, wealth3
)]

# First valid functional observation at or after age 60 is the descriptive
# baseline. This table is unweighted; model pseudo-weights are reported as
# available rather than used to redefine the sample.
setorder(panel, cohort, person_id, interview_time, wave)
baseline <- panel[
  in_wave == 1 & observed_state == TRUE & is.finite(age) & age >= 60 & age <= 95,
  .SD[1L], by = .(cohort, country, person_id)
]
baseline <- merge(
  baseline,
  person_ses,
  by = c("cohort", "person_id"), all.x = TRUE, sort = FALSE
)

count_percent <- function(condition, denominator) {
  c(n = sum(condition, na.rm = TRUE), percent = 100 * sum(condition, na.rm = TRUE) / denominator)
}

summary <- baseline[, {
  n_total <- .N
  proxy_denom <- sum(!is.na(proxy))
  list(
    baseline_persons = n_total,
    age_mean = mean(age),
    age_sd = sd(age),
    women_n = sum(female == 1, na.rm = TRUE),
    women_percent = 100 * mean(female == 1, na.rm = TRUE),
    functional_difficulty_n = sum(difficulty == 1, na.rm = TRUE),
    functional_difficulty_percent = 100 * mean(difficulty == 1, na.rm = TRUE),
    proxy_status_observed_n = proxy_denom,
    proxy_interview_n = sum(proxy == 1, na.rm = TRUE),
    proxy_interview_percent_among_observed = if (proxy_denom) 100 * sum(proxy == 1, na.rm = TRUE) / proxy_denom else NA_real_,
    positive_origin_weight_n = sum(is.finite(respondent_weight) & respondent_weight > 0),
    positive_origin_weight_percent = 100 * mean(is.finite(respondent_weight) & respondent_weight > 0),
    education_low_n = sum(education3_fixed == "low", na.rm = TRUE),
    education_middle_n = sum(education3_fixed == "middle", na.rm = TRUE),
    education_high_n = sum(education3_fixed == "high", na.rm = TRUE),
    education_missing_n = sum(is.na(education3_fixed)),
    wealth_low_n = sum(wealth3 == "low", na.rm = TRUE),
    wealth_middle_n = sum(wealth3 == "middle", na.rm = TRUE),
    wealth_high_n = sum(wealth3 == "high", na.rm = TRUE),
    wealth_missing_n = sum(is.na(wealth3))
  )
}, by = .(cohort, country)]

for (prefix in c("education_low", "education_middle", "education_high", "education_missing",
                 "wealth_low", "wealth_middle", "wealth_high", "wealth_missing")) {
  summary[, paste0(prefix, "_percent") := 100 * get(paste0(prefix, "_n")) / baseline_persons]
}

event_gate <- fread(file.path(gate_dir, "gate0_event_counts_by_cohort.csv"))
living_intervals <- fread(file.path(gate_dir, "interval_length_audit.csv"))[
  interval_exclusion == "included",
  .(cohort, living_transition_intervals = intervals, living_transition_persons = persons)
]
mortality_gate <- fread(file.path(risk_dir, "formal_mortality_gate.csv"))[, .(
  cohort,
  mortality_supported_intervals = supported_intervals,
  mortality_supported_persons = supported_persons,
  mortality_deaths = deaths,
  mortality_disabled_origin_deaths = disabled_deaths
)]

summary <- Reduce(
  function(x, y) merge(x, y, by = "cohort", all.x = TRUE, sort = FALSE),
  list(summary, living_intervals, mortality_gate, event_gate[, -"country"])
)
summary[, .cohort_order := match(cohort, cohort_order)]
setorder(summary, .cohort_order)
summary[, .cohort_order := NULL]

fwrite(summary, file.path(out_dir, "Table1_cohort_overview_and_transition_support.csv"))

long <- melt(
  summary,
  id.vars = c("cohort", "country"),
  variable.name = "characteristic", value.name = "value"
)
fwrite(long, file.path(out_dir, "Table1_cohort_overview_and_transition_support_long.csv"))

caption <- c(
  "# Table 1 caption",
  "",
  "**Table 1 | Descriptive baseline and transition support by cohort.** Baseline is the first valid nine-item functional observation at age 60–95 and is summarised without survey weighting. Functional difficulty denotes at least one difficulty in the common five-ADL plus four-IADL definition. Women's percentages use observations with non-missing sex; proxy percentages use interviews with observable proxy status. Living transition intervals satisfy the frozen 1–4-year rule. Mortality counts come from the formal adjacent-scheduled-wave risk sets, which permit verified deaths within 1 year and therefore differ by design from the Gate-0 common-transition death endpoints. Each formal death row represents one distinct person. A person can contribute more than one recovery or relapse."
)
writeLines(caption, file.path(out_dir, "Table1_caption.md"))

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_manuscript_table1.txt"))
cat("Wrote manuscript Table 1.\n")
