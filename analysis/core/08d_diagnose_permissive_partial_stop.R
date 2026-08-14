#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
derived_dir <- file.path(root, "02_derived", "threshold_sensitivity", "permissive_partial")
out_dir <- file.path(root, "03_outputs", "09_threshold_sensitivity", "permissive_partial")
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mortality <- as.data.table(readRDS(file.path(derived_dir, "formal_mortality_riskset.rds")))
z <- copy(mortality[
  cohort == "MHAS" & primary_mortality_interval == TRUE & origin_state %in% c("D1", "R1", "D2")
])
z[, `:=`(
  event = as.integer(event_death),
  ses = factor(education3_fixed, levels = c("high", "middle", "low")),
  sex_factor = factor(fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown"))),
  proxy_factor = factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self"))),
  period_factor = factor(origin_wave),
  origin_factor = factor(origin_state, levels = c("D1", "R1", "D2")),
  age_model = pmin(pmax(origin_age, 60), 95),
  fit_weight_raw = origin_weight
)]
z <- z[
  !is.na(ses) & is.finite(interval_years) & interval_years > 0 &
    is.finite(fit_weight_raw) & fit_weight_raw > 0 & !is.na(age_model)
]
z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
z[, `:=`(
  ses = droplevels(ses), sex_factor = droplevels(sex_factor),
  proxy_factor = droplevels(proxy_factor), period_factor = droplevels(period_factor),
  origin_factor = droplevels(origin_factor)
)]

support <- z[, .(
  intervals = .N, persons = uniqueN(person_id), deaths = sum(event),
  person_years = sum(interval_years), positive_weight_intervals = sum(fit_weight_raw > 0)
), by = .(ses, origin_state)]

formula <- event ~ ses + ns(age_model, knots = c(70, 80, 90), Boundary.knots = c(60, 95)) +
  origin_factor + sex_factor + period_factor + offset(log(interval_years))
warnings <- character()
fit <- withCallingHandlers(
  glm(
    formula, data = z, family = binomial(link = "cloglog"), weights = fit_weight,
    control = glm.control(maxit = 100, epsilon = 1e-9)
  ),
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

diagnostic <- data.table(
  threshold = "permissive_partial", cohort = "MHAS", exposure = "education",
  process = "death_post", rows = nrow(z), persons = uniqueN(z$person_id),
  deaths = sum(z$event), converged = isTRUE(fit$converged), iterations = fit$iter,
  all_coefficients_finite = all(is.finite(coef(fit))),
  coefficient_abs_max = suppressWarnings(max(abs(coef(fit)), na.rm = TRUE)),
  warning_n = uniqueN(warnings), warnings = paste(unique(warnings), collapse = " | "),
  frozen_decision = "STOP_threshold_model_no_rescue"
)

coefficients <- data.table(
  term = names(coef(fit)), estimate = as.numeric(coef(fit)),
  finite = is.finite(as.numeric(coef(fit)))
)
fwrite(support, file.path(out_dir, "permissive_partial_stop_event_support.csv"))
fwrite(diagnostic, file.path(out_dir, "permissive_partial_model_stop_diagnostic.csv"))
fwrite(coefficients, file.path(out_dir, "permissive_partial_isolated_fit_coefficients.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_permissive_partial_stop_diagnostic.txt"))

print(support)
print(diagnostic)
cat("Frozen decision: STOP this threshold model; no category collapse or alternate fit.\n")
