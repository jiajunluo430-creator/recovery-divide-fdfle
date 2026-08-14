#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: 06_bootstrap_cohort.R COHORT REPLICATES PHASE [START_REPLICATE] [END_REPLICATE]")
}
cohort_name <- toupper(args[[1]])
replicates <- as.integer(args[[2]])
phase <- tolower(args[[3]])
replicate_start <- if (length(args) >= 4L) as.integer(args[[4]]) else 1L
replicate_end <- if (length(args) >= 5L) as.integer(args[[5]]) else replicates
stopifnot(cohort_name %in% c("CHARLS", "HRS", "ELSA", "MHAS"))
stopifnot(is.finite(replicates), replicates >= 1L)
stopifnot(phase %in% c("test", "preview", "final"))
stopifnot(is.finite(replicate_start), is.finite(replicate_end), replicate_start >= 1L,
  replicate_end <= replicates, replicate_start <= replicate_end)

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
derived_dir <- file.path(root, "02_derived")
out_dir <- Sys.getenv(
  "D4_BOOTSTRAP_OUT_DIR",
  unset = file.path(root, "03_outputs", paste0("05_bootstrap_", phase))
)
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

result_path <- file.path(out_dir, paste0("bootstrap_", tolower(cohort_name), "_replicates.csv"))
qc_path <- file.path(out_dir, paste0("bootstrap_", tolower(cohort_name), "_qc.csv"))
cut_path <- file.path(out_dir, paste0("bootstrap_", tolower(cohort_name), "_wealth_cutpoint_audit.csv"))
session_path <- file.path(log_dir, paste0("sessionInfo_bootstrap_", phase, "_", tolower(cohort_name), ".txt"))

function_base <- as.data.table(readRDS(file.path(derived_dir, "formal_function_transition_riskset.rds")))[cohort == cohort_name]
mortality_base <- as.data.table(readRDS(file.path(derived_dir, "formal_mortality_riskset.rds")))[
  cohort == cohort_name & primary_mortality_interval == TRUE
]
person_ses <- as.data.table(readRDS(file.path(derived_dir, "person_fixed_ses_internal.rds")))[cohort == cohort_name]

function_base[, event_recovery := as.integer(origin_state %in% c("D1", "D2") & destination == "R1")]
eligible_ids <- sort(unique(c(function_base$person_id, mortality_base$person_id)))
person_frame <- data.table(person_id = eligible_ids, person_index = seq_along(eligible_ids))
person_ses <- merge(person_frame, person_ses, by = "person_id", all.x = TRUE, sort = FALSE)
function_base[, person_index := person_frame$person_index[match(person_id, person_frame$person_id)]]
mortality_base[, person_index := person_frame$person_index[match(person_id, person_frame$person_id)]]

states <- c("I0", "D1", "R1", "D2", "DEAD")
living_states <- states[1:4]
ses_levels <- c("high", "middle", "low")
ages <- 60:99
allowed <- list(
  I0 = c("D1", "DEAD"),
  D1 = c("R1", "DEAD"),
  R1 = c("D2", "DEAD"),
  D2 = c("R1", "DEAD")
)
process_specs <- list(
  onset = list(source = "function", origins = "I0", destination = "D1", event = "event_onset"),
  recovery = list(source = "function", origins = c("D1", "D2"), destination = "R1", event = "event_recovery"),
  relapse = list(source = "function", origins = "R1", destination = "D2", event = "event_relapse"),
  death_pre = list(source = "mortality", origins = "I0", destination = "DEAD", event = "event_death"),
  death_post = list(source = "mortality", origins = c("D1", "R1", "D2"), destination = "DEAD", event = "event_death")
)

weighted_cut <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) / sum(w) >= p)[1L]]
}

bootstrap_wealth <- function(mult, replicate_id) {
  z <- copy(person_ses)
  z[, boot_mult := mult[person_index]]
  entry <- z[boot_mult > 0 & is.finite(wealth_entry_value)]
  cuts <- entry[, {
    weighted_ok <- any(is.finite(wealth_entry_weight) & wealth_entry_weight > 0)
    ww <- if (weighted_ok) wealth_entry_weight * boot_mult else boot_mult
    list(
      q33 = weighted_cut(wealth_entry_value, ww, 1 / 3),
      q67 = weighted_cut(wealth_entry_value, ww, 2 / 3),
      cut_method = if (weighted_ok) "weighted_empirical" else "unweighted_fallback",
      resampled_person_frequency = sum(boot_mult),
      unique_sampled_persons = .N
    )
  }, by = wealth_entry_wave]
  entry <- cuts[entry, on = "wealth_entry_wave"]
  entry[, wealth3_boot := fifelse(
    wealth_entry_value <= q33, "low",
    fifelse(wealth_entry_value <= q67, "middle", "high")
  )]
  cuts[, `:=`(replicate = replicate_id, cohort = cohort_name)]
  list(map = entry[, .(person_index, wealth3_boot)], cuts = cuts)
}

make_formula <- function(z) {
  terms <- c("ses", "ns(age_model, knots = c(70, 80, 90), Boundary.knots = c(60, 95))")
  if (nlevels(z$origin_factor) > 1L) terms <- c(terms, "origin_factor")
  if (nlevels(z$sex_factor) > 1L) terms <- c(terms, "sex_factor")
  if (nlevels(z$period_factor) > 1L) terms <- c(terms, "period_factor")
  if (nlevels(z$proxy_factor) > 1L) terms <- c(terms, "proxy_factor")
  as.formula(paste("event ~", paste(terms, collapse = " + "), "+ offset(log(interval_years))"))
}

prepare_model_data <- function(fboot, mboot, exposure_name, spec) {
  z <- if (spec$source == "function") {
    copy(fboot[origin_state %in% spec$origins])
  } else {
    copy(mboot[origin_state %in% spec$origins])
  }
  z[, event := as.integer(get(spec$event))]
  z[, ses_value := if (exposure_name == "education") education3_fixed else wealth3_boot]
  z[, ses := factor(ses_value, levels = ses_levels)]
  z[, sex_factor := factor(fifelse(female == 1, "female", fifelse(female == 0, "male", "unknown")))]
  z[, proxy_factor := factor(fifelse(is.na(proxy), "unavailable", fifelse(proxy == 1, "proxy", "self")))]
  z[, period_factor := factor(origin_wave)]
  z[, origin_factor := factor(origin_state, levels = spec$origins)]
  z[, age_model := pmin(pmax(origin_age, 60), 95)]
  z <- z[
    boot_mult > 0 & !is.na(ses) & is.finite(interval_years) & interval_years > 0 &
    is.finite(base_weight) & base_weight > 0 & !is.na(age_model)
  ]
  z[, fit_weight_raw := base_weight * boot_mult]
  z[, fit_weight := fit_weight_raw / mean(fit_weight_raw)]
  z[, `:=`(
    ses = droplevels(ses),
    sex_factor = droplevels(sex_factor),
    proxy_factor = droplevels(proxy_factor),
    period_factor = droplevels(period_factor),
    origin_factor = droplevels(origin_factor)
  )]
  z
}

fit_process_rates <- function(fboot, mboot, exposure_name, process_name, spec) {
  z <- prepare_model_data(fboot, mboot, exposure_name, spec)
  if (!all(c("high", "low") %in% levels(z$ses))) stop("low/high SES level absent")
  formula <- make_formula(z)
  fit <- suppressWarnings(glm(
    formula, data = z, family = binomial(link = "cloglog"),
    weights = fit_weight, control = glm.control(maxit = 100, epsilon = 1e-9)
  ))
  if (!isTRUE(fit$converged) || any(!is.finite(coef(fit)))) stop("nonconvergent/nonfinite model")

  rate_parts <- list()
  for (origin_name in spec$origins) {
    ref <- z[origin_state == origin_name, .(
      fit_weight = sum(fit_weight),
      standardization_rows = .N
    ), by = .(sex_factor, period_factor, proxy_factor, origin_factor)]
    if (!nrow(ref)) stop("empty origin standardization set")
    ref[, `:=`(
      sex_factor = factor(sex_factor, levels = levels(z$sex_factor)),
      period_factor = factor(period_factor, levels = levels(z$period_factor)),
      proxy_factor = factor(proxy_factor, levels = levels(z$proxy_factor)),
      origin_factor = factor(origin_factor, levels = levels(z$origin_factor))
    )]
    for (ses_name in c("high", "low")) {
      for (age_value in ages) {
        nd <- copy(ref)
        nd[, `:=`(
          ses = factor(ses_name, levels = levels(z$ses)),
          age_model = min(age_value, 95),
          interval_years = 1
        )]
        hazard <- exp(as.numeric(predict(fit, newdata = nd, type = "link")))
        annual_hazard <- sum(hazard * nd$fit_weight) / sum(nd$fit_weight)
        if (!is.finite(annual_hazard) || annual_hazard < 0) stop("invalid standardized hazard")
        rate_parts[[length(rate_parts) + 1L]] <- data.table(
          exposure = exposure_name,
          ses = ses_name,
          age = age_value,
          process = process_name,
          origin_state = origin_name,
          destination = spec$destination,
          annual_hazard = annual_hazard
        )
      }
    }
  }
  rbindlist(rate_parts)
}

lambda_array <- function(z) {
  a <- array(
    NA_real_, c(length(ages), length(living_states), length(states)),
    dimnames = list(as.character(ages), living_states, states)
  )
  for (i in seq_len(nrow(z))) {
    a[as.character(z$age[[i]]), z$origin_state[[i]], z$destination[[i]]] <- z$annual_hazard[[i]]
  }
  for (origin_name in names(allowed)) {
    for (destination_name in allowed[[origin_name]]) {
      if (any(!is.finite(a[, origin_name, destination_name]))) stop("missing lambda")
    }
  }
  a
}

transition_matrix <- function(lambda_age) {
  m <- matrix(0, length(states), length(states), dimnames = list(states, states))
  m["DEAD", "DEAD"] <- 1
  for (origin_name in names(allowed)) {
    destinations <- allowed[[origin_name]]
    lambda <- lambda_age[origin_name, destinations]
    total <- sum(lambda)
    move <- if (total > 0) 1 - exp(-total) else 0
    m[origin_name, origin_name] <- 1 - move
    if (total > 0) m[origin_name, destinations] <- move * lambda / total
  }
  if (max(abs(rowSums(m) - 1)) > 1e-10) stop("matrix row-sum failure")
  m
}

life_metrics <- function(lambda) {
  v <- setNames(c(1, 0, 0, 0, 0), states)
  occupancy <- setNames(rep(0, length(states)), states)
  for (age_value in ages) {
    m <- transition_matrix(lambda[as.character(age_value), , , drop = TRUE])
    v_next <- as.numeric(v %*% m)
    names(v_next) <- states
    occupancy <- occupancy + (v + v_next) / 2
    v <- v_next
  }
  c(
    tle = sum(occupancy[living_states]),
    dfle = sum(occupancy[c("I0", "R1")]),
    disabled = sum(occupancy[c("D1", "D2")]),
    residual_age100 = sum(v[living_states])
  )
}

block_pairs <- list(
  onset = list(c("I0", "D1")),
  recovery = list(c("D1", "R1"), c("D2", "R1")),
  relapse = list(c("R1", "D2")),
  post_disability_mortality = list(c("D1", "DEAD"), c("R1", "DEAD"), c("D2", "DEAD")),
  pre_disability_mortality = list(c("I0", "DEAD"))
)
blocks <- names(block_pairs)
permutations <- function(x) {
  if (length(x) == 1L) return(list(x))
  out <- list()
  for (i in seq_along(x)) for (tail in permutations(x[-i])) out[[length(out) + 1L]] <- c(x[[i]], tail)
  out
}
block_permutations <- permutations(blocks)

replace_block <- function(base, donor, block_name) {
  out <- base
  for (pair in block_pairs[[block_name]]) out[, pair[[1]], pair[[2]]] <- donor[, pair[[1]], pair[[2]]]
  out
}

decompose_dfle <- function(high_lambda, low_lambda) {
  high <- life_metrics(high_lambda)
  low <- life_metrics(low_lambda)
  contribution <- setNames(rep(0, length(blocks)), blocks)
  for (perm in block_permutations) {
    hybrid <- high_lambda
    previous <- high[["dfle"]]
    for (block_name in perm) {
      hybrid <- replace_block(hybrid, low_lambda, block_name)
      current <- life_metrics(hybrid)[["dfle"]]
      contribution[[block_name]] <- contribution[[block_name]] + current - previous
      previous <- current
    }
  }
  contribution <- contribution / length(block_permutations)
  gap <- low[["dfle"]] - high[["dfle"]]
  list(high = high, low = low, gap = gap, contribution = contribution,
    closure = sum(contribution) - gap)
}

cohort_offsets <- c(CHARLS = 1L, HRS = 2L, ELSA = 3L, MHAS = 4L)
seed_base <- if (phase == "final") 51000L else if (phase == "preview") 41000L else 31000L

existing_results <- if (file.exists(result_path)) fread(result_path) else data.table()
existing_qc <- if (file.exists(qc_path)) fread(qc_path) else data.table()
completed_replicates <- if (nrow(existing_qc)) unique(existing_qc$replicate) else integer()
result_parts <- if (nrow(existing_results)) list(existing_results) else list()
qc_parts <- if (nrow(existing_qc)) list(existing_qc) else list()
cut_parts <- if (file.exists(cut_path)) list(fread(cut_path)) else list()

cat("Bootstrap cohort=", cohort_name, " phase=", phase, " target=", replicates,
  " range=", replicate_start, "-", replicate_end,
  " eligible_people=", nrow(person_frame), " resume_completed=", length(completed_replicates), "\n", sep = "")

for (b in seq.int(replicate_start, replicate_end)) {
  if (b %in% completed_replicates) next
  start_time <- proc.time()[["elapsed"]]
  seed <- seed_base + cohort_offsets[[cohort_name]] * 100000L + b
  set.seed(seed)
  sampled <- sample.int(nrow(person_frame), nrow(person_frame), replace = TRUE)
  mult <- tabulate(sampled, nbins = nrow(person_frame))
  wealth_boot <- bootstrap_wealth(mult, b)
  cut_parts[[length(cut_parts) + 1L]] <- wealth_boot$cuts

  fboot <- function_base[mult[person_index] > 0]
  fboot[, `:=`(
    boot_mult = mult[person_index],
    base_weight = respondent_weight,
    wealth3_boot = wealth_boot$map$wealth3_boot[match(person_index, wealth_boot$map$person_index)]
  )]
  mboot <- mortality_base[mult[person_index] > 0]
  mboot[, `:=`(
    boot_mult = mult[person_index],
    base_weight = origin_weight,
    wealth3_boot = wealth_boot$map$wealth3_boot[match(person_index, wealth_boot$map$person_index)]
  )]

  for (exposure_name in c("education", "wealth")) {
    status <- "valid"
    error_message <- ""
    row <- NULL
    tryCatch({
      rate_list <- list()
      for (process_name in names(process_specs)) {
        rate_list[[length(rate_list) + 1L]] <- fit_process_rates(
          fboot, mboot, exposure_name, process_name, process_specs[[process_name]]
        )
      }
      rr <- rbindlist(rate_list)
      high_lambda <- lambda_array(rr[ses == "high"])
      low_lambda <- lambda_array(rr[ses == "low"])
      dec <- decompose_dfle(high_lambda, low_lambda)
      if (!is.finite(dec$closure) || abs(dec$closure) > 0.01) stop("Shapley closure failure")
      cval <- dec$contribution
      recovery_relapse <- cval[["recovery"]] + cval[["relapse"]]
      row <- data.table(
        cohort = cohort_name, phase = phase, replicate = b, seed = seed,
        exposure = exposure_name,
        high_tle = dec$high[["tle"]], low_tle = dec$low[["tle"]], gap_tle = dec$low[["tle"]] - dec$high[["tle"]],
        high_dfle = dec$high[["dfle"]], low_dfle = dec$low[["dfle"]], gap_dfle = dec$gap,
        high_disabled = dec$high[["disabled"]], low_disabled = dec$low[["disabled"]],
        gap_disabled = dec$low[["disabled"]] - dec$high[["disabled"]],
        high_residual_age100 = dec$high[["residual_age100"]],
        low_residual_age100 = dec$low[["residual_age100"]],
        contribution_onset = cval[["onset"]],
        contribution_recovery = cval[["recovery"]],
        contribution_relapse = cval[["relapse"]],
        contribution_post_disability_mortality = cval[["post_disability_mortality"]],
        contribution_pre_disability_mortality = cval[["pre_disability_mortality"]],
        recovery_relapse_contribution = recovery_relapse,
        recovery_relapse_percent = if (abs(dec$gap) > 1e-10) 100 * recovery_relapse / dec$gap else NA_real_,
        closure_error = dec$closure
      )
    }, error = function(e) {
      status <<- "failed"
      error_message <<- conditionMessage(e)
    })
    if (!is.null(row)) result_parts[[length(result_parts) + 1L]] <- row
    qc_parts[[length(qc_parts) + 1L]] <- data.table(
      cohort = cohort_name, phase = phase, replicate = b, seed = seed,
      exposure = exposure_name, status = status, error_message = error_message,
      unique_sampled_people = sum(mult > 0),
      elapsed_seconds_replicate_exposure = proc.time()[["elapsed"]] - start_time
    )
    rm(row)
    invisible(gc())
  }

  if (b %% 5L == 0L || b == replicate_end) {
    fwrite(rbindlist(result_parts, fill = TRUE), result_path)
    fwrite(rbindlist(qc_parts, fill = TRUE), qc_path)
    fwrite(rbindlist(cut_parts, fill = TRUE), cut_path)
  }
  latest_qc <- rbindlist(qc_parts, fill = TRUE)[replicate == b]
  cat("replicate=", b, " valid=", sum(latest_qc$status == "valid"), "/2 elapsed=",
    round(proc.time()[["elapsed"]] - start_time, 2), "s\n", sep = "")
  flush.console()
  rm(fboot, mboot, wealth_boot)
  invisible(gc())
}

final_qc <- rbindlist(qc_parts, fill = TRUE)
cat("Completed cohort=", cohort_name, " valid counts:\n", sep = "")
print(final_qc[, .(valid = sum(status == "valid"), failed = sum(status != "valid")), by = exposure])
writeLines(capture.output(sessionInfo()), session_path)
