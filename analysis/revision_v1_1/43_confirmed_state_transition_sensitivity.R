options(stringsAsFactors = FALSE, width = 220, warn = 1)

if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(splines)
  library(sandwich)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
derived_dir <- file.path(revision_root, "02_derived")
out_dir <- file.path(revision_root, "03_outputs", "05_confirmed_state_sensitivity")
log_dir <- file.path(revision_root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

d <- as.data.table(readRDS(file.path(
  derived_dir, "formal_function_transition_riskset_revision_v1_1.rds"
)))
cohorts <- c("CHARLS", "ELSA", "HRS", "MHAS")
ses_levels <- c("high", "middle", "low")

# For each interval ending at wave t, locate the formal adjacent interval that
# begins at t. Its living destination is the next observed functional state.
# This uses only observed states and does not carry a state through death or
# non-response.
confirmation_map <- d[, .(
  cohort,
  person_id,
  destination_wave = origin_wave,
  confirmation_next_state = dest_state_living,
  confirmation_next_time = dest_time_living,
  confirmation_death_before_next = death_before_next
)]
if (confirmation_map[, anyDuplicated(paste(cohort, person_id, destination_wave, sep = "|"))] > 0L) {
  stop("Confirmation map is not unique by cohort/person/destination wave")
}
d <- merge(
  d, confirmation_map,
  by = c("cohort", "person_id", "destination_wave"),
  all.x = TRUE, sort = FALSE
)
d[, confirmation_next_class := fifelse(
  confirmation_next_state %chin% c("D1", "D2"), "difficulty",
  fifelse(confirmation_next_state %chin% c("I0", "R1"), "independent", NA_character_)
)]
d[, `:=`(
  event_onset_confirmed = as.integer(
    event_onset == 1L & !is.na(confirmation_next_class) & confirmation_next_class == "difficulty"
  ),
  event_recovery_confirmed = as.integer(
    event_recovery == 1L & !is.na(confirmation_next_class) & confirmation_next_class == "independent"
  ),
  event_relapse_confirmed = as.integer(
    event_relapse == 1L & !is.na(confirmation_next_class) & confirmation_next_class == "difficulty"
  )
)]

process_specs <- list(
  onset = list(origins = "I0", original_event = "event_onset", confirmed_event = "event_onset_confirmed"),
  recovery = list(origins = c("D1", "D2"), original_event = "event_recovery", confirmed_event = "event_recovery_confirmed"),
  relapse = list(origins = "R1", original_event = "event_relapse", confirmed_event = "event_relapse_confirmed")
)

make_formula <- function(z, include_period = TRUE, include_proxy = TRUE) {
  age_range <- range(z$age_model, na.rm = TRUE)
  knots <- c(70, 80, 90)
  knots <- knots[knots > age_range[[1L]] & knots < age_range[[2L]]]
  age_term <- if (length(knots)) {
    paste0(
      "ns(age_model, knots=c(", paste(knots, collapse = ","),
      "), Boundary.knots=c(60,95))"
    )
  } else {
    "age_model"
  }
  terms <- c("ses", age_term)
  if (nlevels(z$origin_factor) > 1L) terms <- c(terms, "origin_factor")
  if (nlevels(z$sex_factor) > 1L) terms <- c(terms, "sex_factor")
  if (include_period && nlevels(z$period_factor) > 1L) terms <- c(terms, "period_factor")
  if (include_proxy && nlevels(z$proxy_factor) > 1L) terms <- c(terms, "proxy_factor")
  as.formula(paste("event ~", paste(terms, collapse = " + "), "+ offset(log(interval_years))"))
}

support_parts <- list()
estimate_parts <- list()
for (cohort_name in cohorts) {
  for (process_name in names(process_specs)) {
    spec <- process_specs[[process_name]]
    z0 <- copy(d[
      cohort == cohort_name &
        origin_state %chin% spec$origins &
        wealth_time_eligible == TRUE &
        education3_fixed == "low" &
        wealth3 %chin% ses_levels
    ])
    z0[, original_event_value := as.integer(get(spec$original_event))]
    z0[, confirmed_event_value := as.integer(get(spec$confirmed_event))]
    z0[, candidate_unconfirmed := original_event_value == 1L & confirmed_event_value == 0L]
    support_row <- z0[, .(
      intervals_before_confirmation = .N,
      people_before_confirmation = uniqueN(person_id),
      original_events = sum(original_event_value),
      confirmed_events = sum(confirmed_event_value),
      candidate_events_excluded = sum(candidate_unconfirmed),
      candidate_events_without_next_observed_state = sum(
        original_event_value == 1L & is.na(confirmation_next_class)
      ),
      households = uniqueN(household_id)
    )]
    support_row[, `:=`(cohort = cohort_name, process = process_name)]
    setcolorder(support_row, c("cohort", "process", setdiff(names(support_row), c("cohort", "process"))))
    support_parts[[length(support_parts) + 1L]] <- support_row

    # Non-events remain at risk. A candidate event is included only if the
    # destination status is observed again at the next formal adjacent visit.
    z <- z0[candidate_unconfirmed == FALSE]
    z <- z[
      is.finite(interval_years) & interval_years > 0 &
        is.finite(origin_age) &
        is.finite(respondent_weight) & respondent_weight > 0
    ]
    z[, `:=`(
      event = confirmed_event_value,
      ses = factor(wealth3, levels = ses_levels),
      age_model = pmin(pmax(origin_age, 60), 95),
      origin_factor = factor(origin_state, levels = spec$origins),
      sex_factor = factor(fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown"))),
      proxy_factor = factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self"))),
      period_factor = factor(origin_wave),
      fit_weight = respondent_weight / mean(respondent_weight)
    )]
    for (v in c("ses", "origin_factor", "sex_factor", "proxy_factor", "period_factor")) {
      set(z, j = v, value = droplevels(z[[v]]))
    }
    if (!all(c("high", "low") %in% levels(z$ses))) {
      stop("Missing high/low wealth support: ", cohort_name, " / ", process_name)
    }
    cat(
      "Fit ", cohort_name, " / ", process_name,
      ": n=", nrow(z),
      "; events=", sum(z$event),
      "; levels(ses/origin/sex/proxy/period)=",
      paste(vapply(
        c("ses", "origin_factor", "sex_factor", "proxy_factor", "period_factor"),
        function(v) nlevels(z[[v]]), integer(1)
      ), collapse = "/"),
      "\n", sep = ""
    )

    include_period <- TRUE
    include_proxy <- TRUE
    fit_once <- function() {
      model_formula <- make_formula(z, include_period, include_proxy)
      suppressWarnings(glm(
      model_formula,
      data = z,
      family = binomial(link = "cloglog"),
      weights = fit_weight,
      control = glm.control(maxit = 100, epsilon = 1e-9)
      ))
    }
    fit <- fit_once()
    if ((!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) && nlevels(z$proxy_factor) > 1L) {
      include_proxy <- FALSE
      fit <- fit_once()
    }
    if ((!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) && nlevels(z$period_factor) > 1L) {
      include_period <- FALSE
      fit <- fit_once()
    }
    if (!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) {
      stop("Confirmed-state model failure: ", cohort_name, " / ", process_name)
    }

    cluster_type <- "household"
    cluster <- z$household_id
    if (anyNA(cluster) || uniqueN(cluster) < 30L) {
      cluster_type <- "person"
      cluster <- z$person_id
    }
    robust_vcov <- vcovCL(fit, cluster = cluster, type = "HC0", fix = TRUE)
    beta <- coef(fit)[["seslow"]]
    robust_se <- sqrt(diag(robust_vcov))[["seslow"]]
    estimate_parts[[length(estimate_parts) + 1L]] <- data.table(
      cohort = cohort_name,
      process = process_name,
      intervals = nrow(z),
      people = uniqueN(z$person_id),
      confirmed_events = sum(z$event),
      clusters = uniqueN(cluster),
      cluster_type,
      period_adjusted = include_period && nlevels(z$period_factor) > 1L,
      proxy_adjusted = include_proxy && nlevels(z$proxy_factor) > 1L,
      converged = fit$converged,
      low_vs_high_rate_ratio = exp(beta),
      robust_ci_low = exp(beta - 1.96 * robust_se),
      robust_ci_high = exp(beta + 1.96 * robust_se),
      robust_p = 2 * pnorm(abs(beta / robust_se), lower.tail = FALSE)
    )
  }
}

support <- rbindlist(support_parts)
estimates <- rbindlist(estimate_parts)
primary <- fread(file.path(
  revision_root, "03_outputs", "02_revision_point", "transition_model_coefficients.csv"
))[
  module == "wealth_within_low_education" & term == "seslow" & process %chin% names(process_specs),
  .(
    cohort,
    process,
    primary_low_vs_high_rate_ratio = hazard_ratio,
    primary_model_ci_low = model_ci_low,
    primary_model_ci_high = model_ci_high
  )
]
estimates <- primary[estimates, on = .(cohort, process)]

fwrite(support, file.path(out_dir, "confirmed_state_event_support.csv"))
fwrite(estimates, file.path(out_dir, "confirmed_state_cluster_robust_estimates.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_43_confirmed_state.txt"))

cat("Confirmed-state event support:\n")
print(support[order(cohort, process)])
cat("\nConfirmed-state low-vs-high wealth rate ratios:\n")
print(estimates[order(cohort, process)])
