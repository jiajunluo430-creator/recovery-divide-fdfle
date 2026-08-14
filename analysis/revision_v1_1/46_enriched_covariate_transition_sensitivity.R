#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
  library(splines)
  library(sandwich)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 46_enriched_covariate_transition_sensitivity.R CHARLS|ELSA|HRS|MHAS")
}
cohort_name <- toupper(args[[1L]])
stopifnot(cohort_name %in% c("CHARLS", "ELSA", "HRS", "MHAS"))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
derived_dir <- file.path(revision_root, "02_derived")
out_dir <- file.path(revision_root, "03_outputs", "08_enriched_covariate_sensitivity")
log_dir <- file.path(revision_root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(
  log_dir,
  paste0("46_enriched_covariates_", tolower(cohort_name), "_", stamp, ".log")
)
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

chr_id <- function(z) {
  if (is.character(z)) return(z)
  if (is.numeric(z)) return(format(z, scientific = FALSE, trim = TRUE, digits = 22))
  as.character(z)
}
num <- function(z) suppressWarnings(as.numeric(z))
getv <- function(x, name) {
  if (length(name) == 0L || is.na(name) || !name %in% names(x)) {
    return(rep(NA_real_, nrow(x)))
  }
  num(x[[name]])
}
valid_binary <- function(z) {
  z <- num(z)
  fifelse(z %in% 0:1, z, NA_real_)
}
valid_partnered <- function(z) {
  z <- num(z)
  fifelse(z %in% c(1, 2, 3), 1, fifelse(z %in% 4:8, 0, NA_real_))
}

configs <- list(
  CHARLS = list(
    path = Sys.getenv("RECOVERY_DIVIDE_CHARLS_FILE", unset = ""),
    id = "ID", waves = 1:4,
    rural = function(w) paste0("h", w, "rural"),
    insurance = function(w) c(paste0("r", w, "higov"), paste0("r", w, "hipriv")),
    race = NA_character_, hispanic = NA_character_
  ),
  ELSA = list(
    path = Sys.getenv("RECOVERY_DIVIDE_ELSA_FILE", unset = ""),
    id = "idauniq", waves = 1:10,
    rural = function(w) NA_character_,
    insurance = function(w) paste0("r", w, "hipriv"),
    race = "raracem", hispanic = NA_character_
  ),
  HRS = list(
    path = Sys.getenv("RECOVERY_DIVIDE_HRS_FILE", unset = ""),
    id = "hhidpn", waves = 2:16,
    rural = function(w) paste0("r", w, "urbrur"),
    insurance = function(w) paste0("r", w, "henum"),
    race = "raracem", hispanic = "rahispan"
  ),
  MHAS = list(
    path = Sys.getenv("RECOVERY_DIVIDE_MHAS_FILE", unset = ""),
    id = "rahhidnp", waves = 1:6,
    rural = function(w) paste0("h", w, "rural"),
    insurance = function(w) c(paste0("r", w, "htnum"), paste0("r", w, "hipriv")),
    race = NA_character_, hispanic = NA_character_
  )
)
cfg <- configs[[cohort_name]]

person_ses <- as.data.table(readRDS(file.path(derived_dir, "person_fixed_ses_revision_v1_1.rds")))
entry <- person_ses[
  cohort == cohort_name & education3_fixed == "low" & wealth3 %chin% c("high", "low"),
  .(
    person_id,
    wealth3,
    wealth_entry_wave,
    wealth_entry_time,
    wealth_entry_household_id
  )
]
stopifnot(nrow(entry) > 0L, !anyDuplicated(entry$person_id))

requested <- c(cfg$id, cfg$race, cfg$hispanic)
for (w in cfg$waves) {
  requested <- c(
    requested,
    paste0("r", w, c("mstat", "hibpe", "diabe", "cancre", "stroke", "arthre")),
    cfg$rural(w),
    cfg$insurance(w)
  )
}
requested <- unique(na.omit(requested))
meta <- read_dta(cfg$path, n_max = 0)
available <- intersect(requested, names(meta))
missing_schema <- setdiff(requested, names(meta))
raw <- as.data.table(read_dta(cfg$path, col_select = all_of(available)))
raw[, person_id := chr_id(get(cfg$id))]
setkey(raw, person_id)

parts <- list()
for (w in cfg$waves) {
  ids <- entry[wealth_entry_wave == w]
  if (!nrow(ids)) next
  x <- raw[ids, on = "person_id", nomatch = 0]
  conditions <- c("hibpe", "diabe", "cancre", "stroke", "arthre")
  chronic <- do.call(cbind, lapply(conditions, function(s) {
    valid_binary(getv(x, paste0("r", w, s)))
  }))
  colnames(chronic) <- conditions

  rural_raw <- getv(x, cfg$rural(w))
  insurance_names <- cfg$insurance(w)
  insurance_values <- lapply(insurance_names, function(v) getv(x, v))
  insurance_any <- if (cohort_name %in% c("CHARLS", "MHAS")) {
    a <- insurance_values[[1L]]
    b <- insurance_values[[2L]]
    fifelse(
      a > 0 | b == 1, 1,
      fifelse((a == 0 | is.na(a)) & (b == 0 | is.na(b)) & !(is.na(a) & is.na(b)), 0, NA_real_)
    )
  } else if (cohort_name == "HRS") {
    a <- insurance_values[[1L]]
    fifelse(a >= 1, 1, fifelse(a == 0, 0, NA_real_))
  } else {
    valid_binary(insurance_values[[1L]])
  }

  race_raw <- getv(x, cfg$race)
  hispanic_raw <- getv(x, cfg$hispanic)
  race_group <- if (cohort_name == "HRS") {
    fcase(
      hispanic_raw == 1, "Hispanic",
      race_raw == 1, "White non-Hispanic",
      race_raw == 2, "Black non-Hispanic",
      race_raw == 3, "Other non-Hispanic",
      default = "Unknown"
    )
  } else if (cohort_name == "ELSA") {
    fcase(race_raw == 1, "White", race_raw > 1, "Non-White", default = "Unknown")
  } else {
    rep("Not collected in harmonized source", nrow(x))
  }

  rural_group <- if (cohort_name %in% c("CHARLS", "MHAS")) {
    fcase(rural_raw == 1, "Rural", rural_raw == 0, "Urban", default = "Unknown")
  } else if (cohort_name == "HRS") {
    fcase(
      rural_raw == 1, "Urban",
      rural_raw == 2, "Suburban",
      rural_raw == 3, "Rural",
      default = "Unknown"
    )
  } else {
    rep("Not available in harmonized source", nrow(x))
  }

  parts[[length(parts) + 1L]] <- data.table(
    person_id = x$person_id,
    wealth3 = x$wealth3,
    wealth_entry_wave = w,
    partnered = valid_partnered(getv(x, paste0("r", w, "mstat"))),
    chronic_hibp = chronic[, "hibpe"],
    chronic_diabetes = chronic[, "diabe"],
    chronic_cancer = chronic[, "cancre"],
    chronic_stroke = chronic[, "stroke"],
    chronic_arthritis = chronic[, "arthre"],
    rural_group,
    race_group,
    insurance_any
  )
}
covariates <- rbindlist(parts, use.names = TRUE, fill = TRUE)
stopifnot(!anyDuplicated(covariates$person_id))

condition_cols <- c(
  "chronic_hibp", "chronic_diabetes", "chronic_cancer",
  "chronic_stroke", "chronic_arthritis"
)
for (v in condition_cols) {
  m <- mean(covariates[[v]], na.rm = TRUE)
  if (!is.finite(m)) m <- 0
  covariates[, paste0(v, "_imp") := fifelse(is.na(get(v)), m, get(v))]
}
covariates[, chronic_missing_n := rowSums(is.na(.SD)), .SDcols = condition_cols]
covariates[, chronic_count_imp := rowSums(.SD), .SDcols = paste0(condition_cols, "_imp")]
covariates[, `:=`(
  partnered_factor = factor(fcase(
    partnered == 1, "Partnered",
    partnered == 0, "Not partnered",
    default = "Unknown"
  )),
  rural_factor = factor(rural_group),
  race_factor = factor(race_group),
  insurance_factor = factor(fcase(
    insurance_any == 1, "Insured indicator yes",
    insurance_any == 0, "Insured indicator no",
    default = "Unknown"
  ))
)]

missingness <- covariates[, .(
  people = .N,
  partnered_missing_n = sum(is.na(partnered)),
  partnered_missing_percent = 100 * mean(is.na(partnered)),
  any_chronic_item_missing_n = sum(chronic_missing_n > 0),
  any_chronic_item_missing_percent = 100 * mean(chronic_missing_n > 0),
  rural_unknown_n = sum(grepl("Unknown|Not available|Not collected", rural_group)),
  rural_unknown_percent = 100 * mean(grepl("Unknown|Not available|Not collected", rural_group)),
  race_unknown_n = sum(grepl("Unknown|Not available|Not collected", race_group)),
  race_unknown_percent = 100 * mean(grepl("Unknown|Not available|Not collected", race_group)),
  insurance_unknown_n = sum(is.na(insurance_any)),
  insurance_unknown_percent = 100 * mean(is.na(insurance_any)),
  chronic_count_mean = mean(chronic_count_imp),
  chronic_count_sd = sd(chronic_count_imp)
), by = wealth3]
missingness[, `:=`(
  cohort = cohort_name,
  covariate_measurement = "at_first_valid_wealth_wave",
  insurance_semantics = fifelse(
    cohort_name == "ELSA", "private_health_insurance_indicator",
    fifelse(cohort_name == "HRS", "number_of_health_insurance_plans_gt_zero", "public_or_private_coverage_indicator")
  )
)]
fwrite(
  missingness,
  file.path(out_dir, paste0("covariate_missingness_", tolower(cohort_name), ".csv"))
)

schema_audit <- data.table(
  cohort = cohort_name,
  source_path = cfg$path,
  requested_variables = length(requested),
  available_variables = length(available),
  missing_schema_n = length(missing_schema),
  missing_schema = paste(missing_schema, collapse = ";")
)
fwrite(
  schema_audit,
  file.path(out_dir, paste0("covariate_schema_audit_", tolower(cohort_name), ".csv"))
)

function_risk <- as.data.table(readRDS(file.path(derived_dir, "formal_function_transition_riskset_revision_v1_1.rds")))
mortality_risk <- as.data.table(readRDS(file.path(derived_dir, "formal_mortality_riskset_revision_v1_1.rds")))

process_specs <- list(
  onset = list(source = "function", origins = "I0", event = "event_onset"),
  recovery = list(source = "function", origins = c("D1", "D2"), event = "event_recovery"),
  relapse = list(source = "function", origins = "R1", event = "event_relapse"),
  death_pre = list(source = "mortality", origins = "I0", event = "event_death"),
  death_post = list(source = "mortality", origins = c("D1", "R1", "D2"), event = "event_death")
)

prepare_process <- function(process_name) {
  spec <- process_specs[[process_name]]
  if (spec$source == "function") {
    z <- copy(function_risk[
      cohort == cohort_name & wealth_time_eligible == TRUE & education3_fixed == "low" &
        wealth3 %chin% c("high", "low") & origin_state %chin% spec$origins
    ])
    z[, fit_weight_raw := respondent_weight]
  } else {
    z <- copy(mortality_risk[
      cohort == cohort_name & primary_mortality_interval == TRUE &
        wealth_time_eligible == TRUE & education3_fixed == "low" &
        wealth3 %chin% c("high", "low") & origin_state %chin% spec$origins
    ])
    z[, fit_weight_raw := origin_weight]
  }
  z <- covariates[z, on = "person_id", nomatch = 0]
  z[, `:=`(
    event = as.integer(get(spec$event)),
    ses = factor(wealth3, levels = c("high", "low")),
    age_model = pmin(pmax(origin_age, 60), 95),
    age_centered = pmin(pmax(origin_age, 60), 95) - 70,
    sex_factor = factor(fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown"))),
    proxy_factor = factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self"))),
    period_factor = factor(origin_wave),
    origin_factor = factor(origin_state, levels = spec$origins),
    cluster_id = factor(household_id)
  )]
  z <- z[
    !is.na(ses) & is.finite(interval_years) & interval_years > 0 &
      is.finite(fit_weight_raw) & fit_weight_raw > 0 & !is.na(cluster_id)
  ]
  z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
  factor_vars <- c(
    "ses", "sex_factor", "proxy_factor", "period_factor", "origin_factor",
    "partnered_factor", "rural_factor", "race_factor", "insurance_factor"
  )
  for (v in factor_vars) set(z, j = v, value = droplevels(z[[v]]))
  z
}

build_formula <- function(z, enriched, include_period = TRUE, include_proxy = TRUE, use_spline = TRUE) {
  age_range <- range(z$age_model, na.rm = TRUE)
  knots <- c(70, 80, 90)
  knots <- knots[knots > age_range[[1L]] & knots < age_range[[2L]]]
  age_term <- if (use_spline && length(knots)) {
    paste0(
      "ns(age_model, knots=c(", paste(knots, collapse = ","),
      "), Boundary.knots=c(60,95))"
    )
  } else {
    "age_centered"
  }
  terms <- c("ses", age_term)
  for (v in c("origin_factor", "sex_factor")) {
    if (nlevels(z[[v]]) > 1L) terms <- c(terms, v)
  }
  if (include_period && nlevels(z$period_factor) > 1L) terms <- c(terms, "period_factor")
  if (include_proxy && nlevels(z$proxy_factor) > 1L) terms <- c(terms, "proxy_factor")
  if (enriched) {
    terms <- c(terms, "chronic_count_imp", "chronic_missing_n")
    for (v in c("partnered_factor", "rural_factor", "race_factor", "insurance_factor")) {
      if (nlevels(z[[v]]) > 1L) terms <- c(terms, v)
    }
  }
  as.formula(paste("event ~", paste(terms, collapse = " + "), "+ offset(log(interval_years))"))
}

fit_one <- function(z, enriched) {
  include_period <- TRUE
  include_proxy <- TRUE
  use_spline <- TRUE
  fit_once <- function() suppressWarnings(glm(
      build_formula(z, enriched, include_period, include_proxy, use_spline),
      data = z,
      family = binomial(link = "cloglog"),
      weights = fit_weight,
      control = glm.control(maxit = 500, epsilon = 1e-8)
    ))
  acceptable <- function(fit) {
    isTRUE(fit$converged) && "seslow" %in% names(coef(fit)) && is.finite(coef(fit)[["seslow"]])
  }
  fit <- fit_once()
  if (!acceptable(fit) && nlevels(z$proxy_factor) > 1L) {
    include_proxy <- FALSE
    fit <- fit_once()
  }
  if (!acceptable(fit) && nlevels(z$period_factor) > 1L) {
    include_period <- FALSE
    fit <- fit_once()
  }
  if (!acceptable(fit) && use_spline) {
    use_spline <- FALSE
    fit <- fit_once()
  }
  if (!acceptable(fit)) {
    coef_rank <- sort(abs(coef(fit)), decreasing = TRUE, na.last = TRUE)
    coef_diag <- paste(
      names(head(coef_rank, 8L)),
      format(head(coef_rank, 8L), digits = 5),
      collapse = " | "
    )
    stop(
      "Enriched transition model convergence/term failure: process=", process_name,
      " enriched=", enriched, " converged=", fit$converged,
      " iterations=", fit$iter, " ses_term_present=", "seslow" %in% names(coef(fit)),
      " largest_abs_coefficients=", coef_diag,
      " formula=", paste(deparse(formula(fit)), collapse = " ")
    )
  }
  vc <- vcovCL(fit, cluster = z$cluster_id, type = "HC0")
  beta <- coef(fit)[["seslow"]]
  se <- sqrt(vc["seslow", "seslow"])
  data.table(
    adjustment = if (enriched) "enriched_at_wealth_entry" else "base_age_sex_wave_proxy",
    log_intensity_ratio = beta,
    robust_se = se,
    intensity_ratio = exp(beta),
    ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se),
    intervals = nrow(z),
    people = uniqueN(z$person_id),
    households = uniqueN(z$household_id),
    events = sum(z$event),
    converged = fit$converged,
    period_adjusted = include_period && nlevels(z$period_factor) > 1L,
    proxy_adjusted = include_proxy && nlevels(z$proxy_factor) > 1L,
    age_function = if (use_spline) "restricted_cubic_spline" else "linear_centered_at_70",
    formula_used = paste(deparse(formula(fit)), collapse = " ")
  )
}

results <- list()
support <- list()
for (process_name in names(process_specs)) {
  cat("Process: ", process_name, "\n", sep = "")
  z <- prepare_process(process_name)
  base <- fit_one(z, FALSE)
  enriched <- fit_one(z, TRUE)
  ans <- rbindlist(list(base, enriched))
  ans[, `:=`(
    cohort = cohort_name,
    process = process_name,
    contrast = "low_vs_high_wealth_within_low_education"
  )]
  results[[process_name]] <- ans
  support[[process_name]] <- z[, .(
    intervals = .N,
    people = uniqueN(person_id),
    households = uniqueN(household_id),
    events = sum(event),
    person_years = sum(interval_years)
  ), by = wealth3][, `:=`(cohort = cohort_name, process = process_name)]
}

results <- rbindlist(results, use.names = TRUE, fill = TRUE)
results[, attenuation_percent_log_scale := {
  b <- log_intensity_ratio[adjustment == "base_age_sex_wave_proxy"]
  100 * (log_intensity_ratio - b) / abs(b)
}, by = process]
setcolorder(results, c(
  "cohort", "process", "contrast", "adjustment", "intensity_ratio",
  "ci_low", "ci_high", "robust_se", "intervals", "people", "households",
  "events", "attenuation_percent_log_scale", "converged"
))
fwrite(
  results,
  file.path(out_dir, paste0("enriched_transition_results_", tolower(cohort_name), ".csv"))
)
fwrite(
  rbindlist(support),
  file.path(out_dir, paste0("enriched_transition_support_", tolower(cohort_name), ".csv"))
)

writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, paste0("sessionInfo_46_enriched_", tolower(cohort_name), ".txt"))
)

cat("\nEnriched covariate transition sensitivity:\n")
print(results)
cat("Completed: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
