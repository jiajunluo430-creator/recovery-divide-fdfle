#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
out_dir <- file.path(root, "03_outputs", "14_recovery_exhaustion_support")
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

risk <- as.data.table(readRDS(file.path(root, "02_derived", "formal_function_transition_riskset.rds")))[
  origin_state %in% c("D1", "D2") & wealth3 %in% c("high", "middle", "low")
]
risk[, event_recovery := as.integer(destination == "R1")]

support <- risk[, .(
  intervals = .N,
  persons = uniqueN(person_id),
  recovery_events = sum(event_recovery),
  person_years = sum(interval_years)
), by = .(cohort, origin_state, wealth3)]

grid <- CJ(
  cohort = c("CHARLS", "ELSA", "HRS", "MHAS"),
  origin_state = c("D1", "D2"),
  wealth3 = c("high", "middle", "low"),
  unique = TRUE
)
support <- merge(grid, support, by = c("cohort", "origin_state", "wealth3"), all.x = TRUE)
for (column_name in c("intervals", "persons", "recovery_events", "person_years")) {
  set(support, which(is.na(support[[column_name]])), column_name, 0)
}

decisive <- support[wealth3 %in% c("high", "low")]
gate <- decisive[, .(
  decisive_cells = .N,
  minimum_recovery_events = min(recovery_events),
  all_four_cells_ge_100_recoveries = .N == 4L & min(recovery_events) >= 100L
), by = cohort]
gate[, support_decision := fifelse(
  all_four_cells_ge_100_recoveries,
  "SUPPORTED_FOR_POINT_MODEL",
  "STOP_DECISIVE_CELL_SUPPORT"
)]

overall <- data.table(
  supported_cohorts = gate[all_four_cells_ge_100_recoveries == TRUE, .N],
  required_supported_cohorts = 3L
)
overall[, module_decision := fifelse(
  supported_cohorts >= required_supported_cohorts,
  "SUPPORTED_FOR_POINT_MODEL",
  "STOP_TWO_OR_MORE_COHORT_FAILURES"
)]

support[, .wealth_order := match(wealth3, c("high", "middle", "low"))]
setorder(support, cohort, origin_state, .wealth_order)
support[, .wealth_order := NULL]
setorder(gate, cohort)
fwrite(support, file.path(out_dir, "recovery_origin_wealth_support.csv"))
fwrite(gate, file.path(out_dir, "recovery_exhaustion_support_gate_by_cohort.csv"))
fwrite(overall, file.path(out_dir, "recovery_exhaustion_support_decision.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, "sessionInfo_recovery_exhaustion_support.txt")
)

print(gate)
print(overall)
