#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
out_dir <- file.path(root, "03_outputs", "11_exploratory_upgrade_point")
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

function_risk <- as.data.table(readRDS(file.path(root, "02_derived", "formal_function_transition_riskset.rds")))
z <- copy(function_risk[cohort == "ELSA" & origin_state == "I0"])
z[, `:=`(
  event = as.integer(event_onset),
  ses = factor(wealth3, levels = c("high", "middle", "low")),
  sex_factor = factor(fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown"))),
  proxy_factor = factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self"))),
  period_factor = factor(origin_wave),
  age_model = pmin(pmax(origin_age, 60), 95),
  age_interaction = (pmin(pmax(origin_age, 60), 90) - 70) / 10,
  fit_weight_raw = respondent_weight
)]
z <- z[
  !is.na(ses) & is.finite(interval_years) & interval_years > 0 &
    is.finite(fit_weight_raw) & fit_weight_raw > 0 & !is.na(age_model)
]
z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
for (column_name in c("ses", "sex_factor", "proxy_factor", "period_factor")) {
  set(z, j = column_name, value = droplevels(z[[column_name]]))
}

formula <- event ~ ses +
  ns(age_model, knots = c(70, 80, 90), Boundary.knots = c(60, 95)) +
  ses:age_interaction + sex_factor + period_factor + proxy_factor +
  offset(log(interval_years))

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

finite_coefficients <- all(is.finite(coef(fit)))
status <- if (!isTRUE(fit$converged) || !finite_coefficients) {
  "STOP_SECOND_IDENTICAL_FAILURE"
} else {
  "NONREPRODUCIBLE_FAILURE_REVIEW"
}

qc <- data.table(
  module = "wealth_age_varying",
  cohort = "ELSA",
  process = "onset",
  rows = nrow(z),
  persons = uniqueN(z$person_id),
  events = sum(z$event),
  converged = isTRUE(fit$converged),
  iterations = fit$iter,
  finite_coefficients = finite_coefficients,
  coefficient_abs_max = if (finite_coefficients) max(abs(coef(fit))) else Inf,
  warning_n = length(unique(warnings)),
  warnings = paste(unique(warnings), collapse = " | "),
  status = status,
  formula = paste(deparse(formula), collapse = " ")
)

fwrite(qc, file.path(out_dir, "age_interaction_elsa_onset_recheck.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_age_interaction_recheck.txt"))
print(qc)
