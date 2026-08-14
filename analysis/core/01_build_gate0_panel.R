#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
derived_dir <- file.path(root, "02_derived")
out_dir <- file.path(root, "03_outputs", "01_gate0")
log_dir <- file.path(root, "06_logs")
for (d in c(derived_dir, out_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("01_build_gate0_panel_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

num <- function(z) suppressWarnings(as.numeric(z))
chr_id <- function(z) {
  if (is.character(z)) return(z)
  if (is.numeric(z)) return(format(z, scientific = FALSE, trim = TRUE, digits = 22))
  as.character(z)
}
getv <- function(x, name) {
  if (!length(name) || is.na(name) || !name %in% names(x)) return(rep(NA_real_, nrow(x)))
  num(x[[name]])
}
valid01 <- function(z) fifelse(num(z) %in% 0:1, num(z), NA_real_)
valid_age <- function(z) fifelse(is.finite(num(z)) & num(z) >= 0 & num(z) <= 120, num(z), NA_real_)
valid_weight <- function(z) fifelse(is.finite(num(z)) & num(z) > 0, num(z), NA_real_)
valid_wealth <- function(z) fifelse(is.finite(num(z)), num(z), NA_real_)
valid_proxy <- function(z) {
  z <- num(z)
  fifelse(z == 0, 0, fifelse(z %in% 1:3, 1, NA_real_))
}
valid_year <- function(z) fifelse(is.finite(num(z)) & num(z) >= 1900 & num(z) <= 2030, num(z), NA_real_)
valid_month <- function(z) fifelse(num(z) %in% 1:12, num(z), NA_real_)
last_nonmissing <- function(z) {
  idx <- which(!is.na(z))
  if (length(idx)) z[idx[[length(idx)]]] else z[NA_integer_]
}

education3 <- function(z, cohort) {
  z <- num(z)
  if (cohort == "HRS") {
    return(fifelse(z == 1, "low", fifelse(z %in% 2:4, "middle", fifelse(z == 5, "high", NA_character_))))
  }
  fifelse(z %in% c(0, 1), "low", fifelse(z == 2, "middle", fifelse(z == 3, "high", NA_character_)))
}

configs <- list(
  CHARLS = list(
    cohort = "CHARLS", country = "China",
    path = Sys.getenv("RECOVERY_DIVIDE_CHARLS_FILE", unset = ""),
    id = "ID", waves = 1:4, years = c(`1` = 2011, `2` = 2013, `3` = 2015, `4` = 2018),
    age = function(w) paste0("r", w, "agey"), educ = "raeducl",
    wealth = function(w) paste0("h", w, "atotb"), proxy = function(w) NA_character_,
    weight = function(w) paste0("r", w, "wtresp"),
    iwyear = function(w) NA_character_, iwmonth = function(w) NA_character_,
    death_year = "radyear", death_month = "radmonth"
  ),
  HRS = list(
    cohort = "HRS", country = "United States",
    path = Sys.getenv("RECOVERY_DIVIDE_HRS_FILE", unset = ""),
    id = "hhidpn", waves = 2:16, years = setNames(seq(1994, 2022, by = 2), 2:16),
    age = function(w) paste0("r", w, "agey_b"), educ = "raeduc",
    wealth = function(w) paste0("h", w, "atotw"), proxy = function(w) paste0("r", w, "proxy"),
    weight = function(w) paste0("r", w, "wtresp"),
    iwyear = function(w) NA_character_, iwmonth = function(w) NA_character_,
    death_year = "radyear", death_month = "radmonth"
  ),
  ELSA = list(
    cohort = "ELSA", country = "England",
    path = Sys.getenv("RECOVERY_DIVIDE_ELSA_FILE", unset = ""),
    id = "idauniq", waves = 1:10,
    years = c(`1` = 2002, `2` = 2004, `3` = 2006, `4` = 2008, `5` = 2010, `6` = 2012, `7` = 2014, `8` = 2016, `9` = 2018, `10` = 2021),
    age = function(w) paste0("r", w, "agey"), educ = "raeducl",
    wealth = function(w) paste0("h", w, "atotb"), proxy = function(w) paste0("r", w, "proxy"),
    weight = function(w) paste0("r", w, "wtresp"),
    iwyear = function(w) paste0("r", w, "iwy"), iwmonth = function(w) paste0("r", w, "iwm"),
    death_year = NA_character_, death_month = NA_character_
  ),
  MHAS = list(
    cohort = "MHAS", country = "Mexico",
    path = Sys.getenv("RECOVERY_DIVIDE_MHAS_FILE", unset = ""),
    id = "rahhidnp", waves = 1:6,
    years = c(`1` = 2001, `2` = 2003, `3` = 2012, `4` = 2015, `5` = 2018, `6` = 2021),
    age = function(w) paste0("r", w, "agey"), educ = "raeducl",
    wealth = function(w) paste0("h", w, "atotb"), proxy = function(w) paste0("r", w, "proxy"),
    weight = function(w) paste0("r", w, "wtresp"),
    iwyear = function(w) paste0("r", w, "iwy"), iwmonth = function(w) paste0("r", w, "iwm"),
    death_year = "radyear", death_month = "radmonth"
  )
)

adl_tokens <- c("dressa", "batha", "eata", "beda", "toilta")
iadl_tokens <- c("shopa", "mealsa", "medsa", "moneya")
all_tokens <- c(adl_tokens, iadl_tokens)

panel_parts <- list()
schema_parts <- list()
value_count_parts <- list()
value_label_parts <- list()
source_parts <- list()

cat("Started: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
cat("Contract: FROZEN_ANALYSIS_CONTRACT_v1.0.md\n")

for (cfg in configs) {
  cat("\nReading ", cfg$cohort, " metadata and selected columns...\n", sep = "")
  meta <- read_dta(cfg$path, n_max = 0)
  requested <- c(cfg$id, "ragender", cfg$educ, cfg$death_year, cfg$death_month)
  for (w in cfg$waves) {
    requested <- c(
      requested, paste0("inw", w), paste0("r", w, "iwstat"), cfg$age(w),
      paste0("r", w, all_tokens), cfg$wealth(w), cfg$proxy(w), cfg$weight(w),
      cfg$iwyear(w), cfg$iwmonth(w)
    )
  }
  requested <- unique(na.omit(requested))
  available <- intersect(requested, names(meta))
  schema_parts[[length(schema_parts) + 1L]] <- data.table(
    cohort = cfg$cohort, variable = requested, available = requested %in% available,
    label = vapply(requested, function(v) {
      if (!v %in% names(meta)) return(NA_character_)
      z <- attr(meta[[v]], "label", exact = TRUE)
      if (is.null(z) || !length(z)) NA_character_ else as.character(z[[1]])
    }, character(1))
  )

  label_vars <- available[grepl("iwstat$|proxy$|^raeduc|dressa$|batha$|eata$|beda$|toilta$|shopa$|mealsa$|medsa$|moneya$", available)]
  for (v in label_vars) {
    labs <- attr(meta[[v]], "labels", exact = TRUE)
    if (!is.null(labs) && length(labs)) {
      value_label_parts[[length(value_label_parts) + 1L]] <- data.table(
        cohort = cfg$cohort, variable = v, value = as.numeric(labs), value_label = names(labs)
      )
    }
  }

  x <- as.data.table(read_dta(cfg$path, col_select = all_of(available)))
  person_id <- chr_id(x[[cfg$id]])
  female <- fifelse(getv(x, "ragender") == 2, 1, fifelse(getv(x, "ragender") == 1, 0, NA_real_))
  educ3 <- education3(getv(x, cfg$educ), cfg$cohort)
  death_year <- if (!is.na(cfg$death_year) && cfg$death_year %in% names(x)) valid_year(getv(x, cfg$death_year)) else rep(NA_real_, nrow(x))
  death_month <- if (!is.na(cfg$death_month) && cfg$death_month %in% names(x)) valid_month(getv(x, cfg$death_month)) else rep(NA_real_, nrow(x))
  static_death_time <- fifelse(!is.na(death_year), death_year + (fifelse(!is.na(death_month), death_month, 6.5) - 0.5) / 12, NA_real_)

  source_parts[[length(source_parts) + 1L]] <- data.table(
    cohort = cfg$cohort, country = cfg$country, source_path = cfg$path,
    source_rows = nrow(x), source_columns_total = ncol(meta),
    requested_columns = length(requested), available_columns = length(available),
    exact_schema_pass = all(c(cfg$id, "ragender", cfg$educ) %in% available)
  )

  semantic_vars <- available[grepl("iwstat$|proxy$|^raeduc|dressa$|batha$|eata$|beda$|toilta$|shopa$|mealsa$|medsa$|moneya$", available)]
  for (v in semantic_vars) {
    vv <- num(x[[v]])
    tab <- data.table(value = fifelse(is.na(vv), "<NA>", format(vv, scientific = FALSE, trim = TRUE)))
    tab <- tab[, .(n = .N), by = value]
    tab[, `:=`(cohort = cfg$cohort, variable = v)]
    value_count_parts[[length(value_count_parts) + 1L]] <- tab[, .(cohort, variable, value, n)]
  }

  cohort_panel <- vector("list", length(cfg$waves))
  for (j in seq_along(cfg$waves)) {
    w <- cfg$waves[[j]]
    items <- lapply(all_tokens, function(token) valid01(getv(x, paste0("r", w, token))))
    names(items) <- all_tokens
    item_matrix <- do.call(cbind, items)
    complete_items <- rowSums(!is.na(item_matrix)) == length(all_tokens)
    difficulty_count <- rowSums(item_matrix, na.rm = FALSE)
    difficulty <- fifelse(complete_items, as.integer(difficulty_count > 0), NA_integer_)

    iwy <- valid_year(getv(x, cfg$iwyear(w)))
    iwm <- valid_month(getv(x, cfg$iwmonth(w)))
    nominal_year <- unname(cfg$years[as.character(w)])
    interview_time <- fifelse(!is.na(iwy), iwy + (fifelse(!is.na(iwm), iwm, 6.5) - 0.5) / 12, nominal_year + 0.5)
    time_source <- fifelse(!is.na(iwy), "interview_year_month", "nominal_midyear")

    z <- data.table(
      cohort = cfg$cohort, country = cfg$country, person_id = person_id,
      wave = as.integer(w), nominal_year = nominal_year,
      interview_time = interview_time, time_source = time_source,
      in_wave = valid01(getv(x, paste0("inw", w))),
      iwstat = getv(x, paste0("r", w, "iwstat")),
      age = valid_age(getv(x, cfg$age(w))), female = female,
      education3 = educ3, wealth = valid_wealth(getv(x, cfg$wealth(w))),
      proxy = if (is.na(cfg$proxy(w))) rep(NA_real_, nrow(x)) else valid_proxy(getv(x, cfg$proxy(w))),
      respondent_weight = valid_weight(getv(x, cfg$weight(w))),
      difficulty_count = difficulty_count, function_complete9 = complete_items,
      difficulty = difficulty, static_death_time = static_death_time,
      static_death_source = fifelse(!is.na(static_death_time), "static_death_year_month", NA_character_)
    )
    for (token in all_tokens) z[[token]] <- items[[token]]
    cohort_panel[[j]] <- z
  }
  panel_parts[[length(panel_parts) + 1L]] <- rbindlist(cohort_panel, use.names = TRUE, fill = TRUE)
  rm(x, meta, cohort_panel)
  invisible(gc())
}

panel <- rbindlist(panel_parts, use.names = TRUE, fill = TRUE)

# ELSA wave 11 is available as a raw UKDS core file rather than in the
# Harmonized ELSA source. The frozen 9-item definition is reconstructed from
# the exact raw items. The wave 11 EOL file is a 46%-response subset of known
# deaths and therefore supplements confirmed death timing only; absence from
# EOL is never treated as survival.
elsa_w11_root <- file.path(
  root,
  "00_staging/elsa_wave11/UKDA-5050-stata/stata/stata13_se"
)
elsa_w11_core_path <- file.path(elsa_w11_root, "wave_11_elsa_data_eul.dta")
elsa_w11_eol_path <- file.path(elsa_w11_root, "wave_11_elsa_eol_eul.dta")
stopifnot(file.exists(elsa_w11_core_path), file.exists(elsa_w11_eol_path))

w11_raw_items <- c(
  dressa = "headldr", batha = "headlba", eata = "headlea",
  beda = "headlbe", toilta = "headlwc", shopa = "headlsh",
  mealsa = "headlpr", medsa = "headlme", moneya = "headlmo"
)
w11_core_vars <- c(
  "idauniq", "iintdatm", "iintdaty", "indager", "diint", "w11xwgt",
  unname(w11_raw_items)
)
w11_core_meta <- read_dta(elsa_w11_core_path, n_max = 0)
w11_core_available <- intersect(w11_core_vars, names(w11_core_meta))
if (!all(w11_core_vars %in% w11_core_available)) {
  stop("ELSA wave 11 core is missing frozen required fields: ",
    paste(setdiff(w11_core_vars, w11_core_available), collapse = ", "))
}
w11_core <- as.data.table(read_dta(elsa_w11_core_path, col_select = all_of(w11_core_vars)))
w11_core[, person_id := chr_id(idauniq)]
w11_core[, `:=`(
  core_present = TRUE,
  interview_year = valid_year(iintdaty),
  interview_month = valid_month(iintdatm),
  age_raw = valid_age(indager),
  proxy_raw = fifelse(num(diint) %in% 1:5, 1,
    fifelse(num(diint) == -1, 0, NA_real_)),
  weight_raw = valid_weight(w11xwgt)
)]
for (token in names(w11_raw_items)) {
  w11_core[[token]] <- valid01(w11_core[[w11_raw_items[[token]]]])
}
w11_core <- w11_core[, c(
  list(
    person_id = person_id, core_present = core_present,
    interview_year = interview_year, interview_month = interview_month,
    age_raw = age_raw, proxy_raw = proxy_raw, weight_raw = weight_raw
  ),
  mget(names(w11_raw_items))
)]

w11_eol_vars <- c("idauniq", "wave", "eidatey", "dveidates")
w11_eol_meta <- read_dta(elsa_w11_eol_path, n_max = 0)
w11_eol_available <- intersect(w11_eol_vars, names(w11_eol_meta))
if (!all(c("idauniq", "eidatey") %in% w11_eol_available)) {
  stop("ELSA wave 11 EOL is missing ID or death-year fields")
}
w11_eol <- as.data.table(read_dta(elsa_w11_eol_path, col_select = all_of(w11_eol_available)))
w11_eol[, person_id := chr_id(idauniq)]
w11_eol[, eol_death_year := valid_year(eidatey)]
w11_eol_map <- w11_eol[, .(
  eol_death_year = if (any(!is.na(eol_death_year))) min(eol_death_year, na.rm = TRUE) else NA_real_,
  eol_record_n = .N
), by = person_id]

elsa_prior <- panel[cohort == "ELSA"]
setorder(elsa_prior, person_id, wave)
elsa_lookup <- elsa_prior[, .(
  female_prior = last_nonmissing(female),
  education3_prior = last_nonmissing(education3),
  last_age_prior = last_nonmissing(age),
  last_time_prior = last_nonmissing(interview_time)
), by = person_id]

w11_ids <- data.table(person_id = unique(c(
  elsa_prior$person_id, w11_core$person_id, w11_eol_map$person_id
)))
w11 <- merge(w11_ids, w11_core, by = "person_id", all.x = TRUE, sort = FALSE)
w11 <- merge(w11, w11_eol_map, by = "person_id", all.x = TRUE, sort = FALSE)
w11 <- merge(w11, elsa_lookup, by = "person_id", all.x = TRUE, sort = FALSE)
w11[, interview_time := fifelse(
  !is.na(interview_year),
  interview_year + (fifelse(!is.na(interview_month), interview_month, 6.5) - 0.5) / 12,
  2024.5
)]
w11[, age_derived := fifelse(
  !is.na(age_raw), age_raw,
  fifelse(!is.na(last_age_prior) & !is.na(last_time_prior),
    last_age_prior + interview_time - last_time_prior, NA_real_)
)]
w11_item_matrix <- as.matrix(w11[, ..all_tokens])
w11[, function_complete9 := rowSums(!is.na(w11_item_matrix)) == length(all_tokens)]
w11[, difficulty_count := rowSums(w11_item_matrix, na.rm = FALSE)]
w11[, difficulty := fifelse(function_complete9, as.integer(difficulty_count > 0), NA_integer_)]

w11_panel <- w11[, c(
  list(
    cohort = "ELSA", country = "England", person_id = person_id,
    wave = 11L, nominal_year = 2024,
    interview_time = interview_time,
    time_source = fifelse(!is.na(interview_year), "interview_year_month", "nominal_midyear"),
    in_wave = fifelse(!is.na(core_present) & core_present == TRUE, 1, 0),
    iwstat = fifelse(!is.na(core_present) & core_present == TRUE, 1,
      fifelse(!is.na(eol_death_year), 5, NA_real_)),
    age = valid_age(age_derived), female = female_prior,
    education3 = education3_prior, wealth = NA_real_,
    proxy = proxy_raw, respondent_weight = weight_raw,
    difficulty_count = difficulty_count,
    function_complete9 = function_complete9,
    difficulty = difficulty,
    static_death_time = fifelse(!is.na(eol_death_year), eol_death_year + 0.5, NA_real_),
    static_death_source = fifelse(!is.na(eol_death_year), "elsa_wave11_eol_year", NA_character_)
  ),
  mget(all_tokens)
)]
panel <- rbindlist(list(panel, w11_panel), use.names = TRUE, fill = TRUE)

schema_parts[[length(schema_parts) + 1L]] <- data.table(
  cohort = "ELSA", variable = c(w11_core_vars, w11_eol_vars),
  available = c(w11_core_vars %in% w11_core_available, w11_eol_vars %in% w11_eol_available),
  label = c(
    vapply(w11_core_vars, function(v) {
      z <- attr(w11_core_meta[[v]], "label", exact = TRUE)
      if (is.null(z) || !length(z)) NA_character_ else as.character(z[[1]])
    }, character(1)),
    vapply(w11_eol_vars, function(v) {
      if (!v %in% names(w11_eol_meta)) return(NA_character_)
      z <- attr(w11_eol_meta[[v]], "label", exact = TRUE)
      if (is.null(z) || !length(z)) NA_character_ else as.character(z[[1]])
    }, character(1))
  )
)
source_parts[[length(source_parts) + 1L]] <- data.table(
  cohort = "ELSA", country = "England", source_path = c(elsa_w11_core_path, elsa_w11_eol_path),
  source_rows = c(nrow(w11_core), nrow(w11_eol)),
  source_columns_total = c(ncol(w11_core_meta), ncol(w11_eol_meta)),
  requested_columns = c(length(w11_core_vars), length(w11_eol_vars)),
  available_columns = c(length(w11_core_available), length(w11_eol_available)),
  exact_schema_pass = c(all(w11_core_vars %in% w11_core_available), all(c("idauniq", "eidatey") %in% w11_eol_available))
)

# CHARLS wave 5 (2020) is reconstructed from canonical raw modules under the
# versioned v1.0.2 adapter. Locally derived wave-5 files were used only for the
# pre-transition semantic validation and do not enter this construction.
charls_w5_root <- Sys.getenv("RECOVERY_DIVIDE_CHARLS2020_RAW_ROOT", unset = "")
charls_sample_path <- file.path(charls_w5_root, "Sample_Infor.dta")
charls_health_path <- file.path(charls_w5_root, "Health_Status_and_Functioning.dta")
charls_weight_path <- file.path(charls_w5_root, "Weights.dta")
charls_exit_path <- file.path(charls_w5_root, "Exit_Module.dta")
stopifnot(all(file.exists(c(charls_sample_path, charls_health_path, charls_weight_path, charls_exit_path))))

charls_w5_raw_items <- c(
  dressa = "db001", batha = "db003", eata = "db005", beda = "db007", toilta = "db009",
  shopa = "db016", mealsa = "db014", medsa = "db020", moneya = "db022"
)
charls_sample_vars <- c("ID", "died", "iyear", "imonth")
charls_health_vars <- c("ID", "proxy", unname(charls_w5_raw_items))
charls_weight_vars <- c("ID", "INDV_weight")
charls_exit_vars <- c("ID", "exb001_1", "exb001_2", "exb002")

charls_sample_meta <- read_dta(charls_sample_path, n_max = 0)
charls_health_meta <- read_dta(charls_health_path, n_max = 0)
charls_weight_meta <- read_dta(charls_weight_path, n_max = 0)
charls_exit_meta <- read_dta(charls_exit_path, n_max = 0)
if (!all(charls_sample_vars %in% names(charls_sample_meta)) ||
    !all(charls_health_vars %in% names(charls_health_meta)) ||
    !all(charls_weight_vars %in% names(charls_weight_meta)) ||
    !all(charls_exit_vars %in% names(charls_exit_meta))) {
  stop("CHARLS 2020 canonical raw modules are missing one or more frozen adapter fields")
}

charls_sample <- as.data.table(read_dta(charls_sample_path, col_select = all_of(charls_sample_vars)))
charls_sample[, person_id := chr_id(ID)]
charls_sample <- charls_sample[, .(
  person_id,
  sample_present = TRUE,
  died_w5 = valid01(died),
  interview_year_w5 = valid_year(iyear),
  interview_month_w5 = valid_month(imonth)
)]

charls_health <- as.data.table(read_dta(charls_health_path, col_select = all_of(charls_health_vars)))
charls_health[, person_id := chr_id(ID)]
charls_health[, health_present := TRUE]
charls_health[, proxy_raw_w5 := as.integer(!is.na(num(proxy)) & num(proxy) == 1)]
charls_w5_value_counts <- lapply(unname(charls_w5_raw_items), function(v) {
  vv <- num(charls_health[[v]])
  tab <- data.table(value = fifelse(is.na(vv), "<NA>", format(vv, scientific = FALSE, trim = TRUE)))[, .(n = .N), by = value]
  tab[, `:=`(cohort = "CHARLS", variable = paste0("w5_raw_", v))]
  tab[, .(cohort, variable, value, n)]
})
for (token in names(charls_w5_raw_items)) {
  raw_var <- charls_w5_raw_items[[token]]
  charls_health[[token]] <- fifelse(num(charls_health[[raw_var]]) == 1, 0,
    fifelse(num(charls_health[[raw_var]]) %in% 2:4, 1, NA_real_))
}
charls_health <- charls_health[, c(
  list(person_id = person_id, health_present = health_present, proxy_raw_w5 = proxy_raw_w5),
  mget(all_tokens)
)]

charls_weight <- as.data.table(read_dta(charls_weight_path, col_select = all_of(charls_weight_vars)))
charls_weight[, person_id := chr_id(ID)]
charls_weight <- charls_weight[, .(person_id, weight_raw_w5 = valid_weight(INDV_weight))]

charls_exit <- as.data.table(read_dta(charls_exit_path, col_select = all_of(charls_exit_vars)))
charls_exit[, person_id := chr_id(ID)]
charls_exit[, `:=`(
  exit_year_w5 = valid_year(exb001_1),
  exit_month_w5 = valid_month(exb001_2),
  exit_calendar_w5 = num(exb002)
)]
charls_exit[, static_death_time_w5 := fifelse(
  !is.na(exit_year_w5) & exit_calendar_w5 == 1 & !is.na(exit_month_w5),
  exit_year_w5 + (exit_month_w5 - 0.5) / 12,
  fifelse(!is.na(exit_year_w5), exit_year_w5 + 0.5, NA_real_)
)]
charls_exit[, static_death_source_w5 := fifelse(
  !is.na(exit_year_w5) & exit_calendar_w5 == 1 & !is.na(exit_month_w5),
  "charls2020_exit_solar_year_month",
  fifelse(!is.na(exit_year_w5), "charls2020_exit_year_midpoint_lunar_or_year_only", NA_character_)
)]
charls_exit <- charls_exit[, .(
  person_id, exit_record_w5 = TRUE, static_death_time_w5, static_death_source_w5
)]

charls_prior <- panel[cohort == "CHARLS"]
setorder(charls_prior, person_id, wave)
charls_lookup <- charls_prior[, .(
  female_prior = last_nonmissing(female),
  education3_prior = last_nonmissing(education3),
  last_age_prior = last_nonmissing(age),
  last_time_prior = last_nonmissing(interview_time)
), by = person_id]

charls_w5_ids <- data.table(person_id = unique(c(
  charls_prior$person_id, charls_sample$person_id, charls_exit$person_id
)))
charls_w5 <- merge(charls_w5_ids, charls_sample, by = "person_id", all.x = TRUE, sort = FALSE)
charls_w5 <- merge(charls_w5, charls_health, by = "person_id", all.x = TRUE, sort = FALSE)
charls_w5 <- merge(charls_w5, charls_weight, by = "person_id", all.x = TRUE, sort = FALSE)
charls_w5 <- merge(charls_w5, charls_exit, by = "person_id", all.x = TRUE, sort = FALSE)
charls_w5 <- merge(charls_w5, charls_lookup, by = "person_id", all.x = TRUE, sort = FALSE)
charls_w5[, interview_time := fifelse(
  !is.na(interview_year_w5),
  interview_year_w5 + (fifelse(!is.na(interview_month_w5), interview_month_w5, 6.5) - 0.5) / 12,
  2020.5
)]
charls_w5[, age_derived := fifelse(
  !is.na(last_age_prior) & !is.na(last_time_prior),
  last_age_prior + interview_time - last_time_prior,
  NA_real_
)]
charls_w5_item_matrix <- as.matrix(charls_w5[, ..all_tokens])
charls_w5[, function_complete9 := rowSums(!is.na(charls_w5_item_matrix)) == length(all_tokens)]
charls_w5[, difficulty_count := rowSums(charls_w5_item_matrix, na.rm = FALSE)]
charls_w5[, difficulty := fifelse(function_complete9, as.integer(difficulty_count > 0), NA_integer_)]

charls_w5_panel <- charls_w5[, c(
  list(
    cohort = "CHARLS", country = "China", person_id = person_id,
    wave = 5L, nominal_year = 2020,
    interview_time = interview_time,
    time_source = fifelse(!is.na(interview_year_w5), "interview_year_month", "nominal_midyear"),
    in_wave = fifelse(!is.na(died_w5) & died_w5 == 0, 1, 0),
    iwstat = fifelse(!is.na(died_w5) & died_w5 == 1, 5,
      fifelse(!is.na(died_w5) & died_w5 == 0, 1, NA_real_)),
    age = valid_age(age_derived), female = female_prior,
    education3 = education3_prior, wealth = NA_real_,
    proxy = fifelse(!is.na(health_present) & health_present == TRUE, proxy_raw_w5, NA_real_),
    respondent_weight = weight_raw_w5,
    difficulty_count = difficulty_count,
    function_complete9 = function_complete9,
    difficulty = difficulty,
    static_death_time = static_death_time_w5,
    static_death_source = static_death_source_w5
  ),
  mget(all_tokens)
)]
panel <- rbindlist(list(panel, charls_w5_panel), use.names = TRUE, fill = TRUE)

charls_meta_roles <- list(
  sample = list(meta = charls_sample_meta, vars = charls_sample_vars, path = charls_sample_path, rows = nrow(charls_sample)),
  health = list(meta = charls_health_meta, vars = charls_health_vars, path = charls_health_path, rows = nrow(charls_health)),
  weight = list(meta = charls_weight_meta, vars = charls_weight_vars, path = charls_weight_path, rows = nrow(charls_weight)),
  exit = list(meta = charls_exit_meta, vars = charls_exit_vars, path = charls_exit_path, rows = nrow(charls_exit))
)
for (role in names(charls_meta_roles)) {
  spec <- charls_meta_roles[[role]]
  schema_parts[[length(schema_parts) + 1L]] <- data.table(
    cohort = "CHARLS", variable = paste0("w5_", role, ":", spec$vars), available = TRUE,
    label = vapply(spec$vars, function(v) {
      z <- attr(spec$meta[[v]], "label", exact = TRUE)
      if (is.null(z) || !length(z)) NA_character_ else as.character(z[[1]])
    }, character(1))
  )
  source_parts[[length(source_parts) + 1L]] <- data.table(
    cohort = "CHARLS", country = "China", source_path = spec$path,
    source_rows = spec$rows, source_columns_total = ncol(spec$meta),
    requested_columns = length(spec$vars), available_columns = length(spec$vars),
    exact_schema_pass = TRUE
  )
}
for (v in unname(charls_w5_raw_items)) {
  value_count_parts[[length(value_count_parts) + 1L]] <- charls_w5_value_counts[[match(v, unname(charls_w5_raw_items))]]
  labs <- attr(charls_health_meta[[v]], "labels", exact = TRUE)
  if (!is.null(labs) && length(labs)) {
    value_label_parts[[length(value_label_parts) + 1L]] <- data.table(
      cohort = "CHARLS", variable = paste0("w5_raw_", v),
      value = as.numeric(labs), value_label = names(labs)
    )
  }
}
rm(
  charls_sample_meta, charls_health_meta, charls_weight_meta, charls_exit_meta,
  charls_sample, charls_health, charls_weight, charls_exit,
  charls_w5_ids, charls_w5, charls_w5_item_matrix, charls_w5_value_counts,
  charls_prior, charls_lookup, charls_meta_roles
)
invisible(gc())

panel <- panel[!is.na(person_id) & person_id != ""]
setorder(panel, cohort, person_id, wave)

# Broadcast the earliest documented static/EOL death time across all records
# for each person so that the last living interval can terminate at death.
panel[, c("static_death_time", "static_death_source") := {
  idx <- which(!is.na(static_death_time))
  if (!length(idx)) list(NA_real_, NA_character_) else {
    k <- idx[[which.min(static_death_time[idx])]]
    list(static_death_time[[k]], static_death_source[[k]])
  }
}, by = .(cohort, person_id)]

# Attach a person-level raw death time. Static death year/month is preferred;
# otherwise use the first wave carrying the verified Gateway/RAND death code 5.
panel[, iwstat_death_time := {
  z <- interview_time[iwstat == 5 & !is.na(iwstat)]
  if (length(z)) min(z, na.rm = TRUE) else NA_real_
}, by = .(cohort, person_id)]
panel[, death_time_raw := fifelse(!is.na(static_death_time), static_death_time, iwstat_death_time)]
panel[, death_source_raw := fifelse(!is.na(static_death_time), static_death_source,
  fifelse(!is.na(iwstat_death_time), "iwstat_5", NA_character_))]

# Observed-history state construction; no lifetime pre-cohort history is inferred.
panel[, observed_state := !is.na(difficulty) & in_wave == 1]
panel[, ever_difficulty_prior := shift(cummax(fifelse(observed_state & difficulty == 1, 1L, 0L)), fill = 0L), by = .(cohort, person_id)]
panel[, recovery_observation := as.integer(observed_state & difficulty == 0 & ever_difficulty_prior == 1L)]
panel[, ever_recovery_prior := shift(cummax(recovery_observation), fill = 0L), by = .(cohort, person_id)]
panel[, history_state := fifelse(!observed_state, NA_character_,
  fifelse(difficulty == 0,
    fifelse(ever_difficulty_prior == 1L, "R1", "I0"),
    fifelse(ever_recovery_prior == 1L, "D2", "D1")
  ))]

# Death dates are bounded after the last living interview as required by the
# contract. For sources without interview month (CHARLS/RAND HRS), a raw death
# in the same/adjacent nominal survey year can precede the artificial midyear
# timestamp even though the true interview preceded death; those cases are
# moved one day after the last observed interview and flagged. Exact-date
# contradictions are not repaired: the raw death time is set unresolved unless
# a later iwstat-5 time supplies a valid alternative.
panel[, `:=`(
  last_living_time = if (any(observed_state)) max(interview_time[observed_state], na.rm = TRUE) else NA_real_,
  last_living_time_source = if (any(observed_state)) time_source[which.max(fifelse(observed_state, interview_time, -Inf))] else NA_character_
), by = .(cohort, person_id)]
panel[, death_conflict_raw := !is.na(death_time_raw) & !is.na(last_living_time) & death_time_raw <= last_living_time]
panel[, nominal_time_boundary_adjustment := death_conflict_raw & last_living_time_source == "nominal_midyear" &
  death_time_raw >= last_living_time - 1.05]
panel[, coarse_death_year_boundary_adjustment := death_conflict_raw &
  death_source_raw == "elsa_wave11_eol_year" &
  floor(death_time_raw) == floor(last_living_time)]
panel[, valid_later_iwstat_death := !is.na(iwstat_death_time) & !is.na(last_living_time) & iwstat_death_time > last_living_time]
panel[, death_time := fifelse(!death_conflict_raw, death_time_raw,
  fifelse(nominal_time_boundary_adjustment, last_living_time + (1 / 365.25),
    fifelse(coarse_death_year_boundary_adjustment, last_living_time + (1 / 365.25),
      fifelse(valid_later_iwstat_death, iwstat_death_time, NA_real_))))]
panel[, death_source := fifelse(!death_conflict_raw, death_source_raw,
  fifelse(nominal_time_boundary_adjustment, "bounded_nominal_death_after_last_interview",
    fifelse(coarse_death_year_boundary_adjustment, "bounded_coarse_death_year_after_last_interview",
      fifelse(valid_later_iwstat_death, "later_iwstat_5_after_raw_conflict", NA_character_))))]
panel[, death_time_qc := fifelse(!death_conflict_raw, "raw_order_valid",
  fifelse(nominal_time_boundary_adjustment, "bounded_nominal_time_artifact",
    fifelse(coarse_death_year_boundary_adjustment, "bounded_coarse_death_year",
      fifelse(valid_later_iwstat_death, "replaced_by_later_iwstat", "unresolved_exact_date_conflict"))))]

# Flag any residual living interview after the final accepted death time.
panel[, living_after_death := observed_state & !is.na(death_time) & interview_time > death_time + (1 / 365.25)]

# The next scheduled wave exists in the wide file even when function is missing;
# it supplies the administrative observation boundary.
panel[, `:=`(
  next_scheduled_wave = shift(wave, type = "lead"),
  next_scheduled_time = shift(interview_time, type = "lead")
), by = .(cohort, person_id)]

# Build consecutive observed-state intervals, then replace the next living state by
# death if a verified death lies first. Death after the final living observation is
# added when it lies within the cohort observation window.
obs <- panel[observed_state == TRUE]
setorder(obs, cohort, person_id, interview_time, wave)
obs[, `:=`(
  dest_wave_living = shift(wave, type = "lead"),
  dest_time_living = shift(interview_time, type = "lead"),
  dest_age_living = shift(age, type = "lead"),
  dest_state_living = shift(history_state, type = "lead")
), by = .(cohort, person_id)]

study_end <- panel[, .(study_end_time = max(interview_time, na.rm = TRUE)), by = cohort]
intervals <- study_end[obs, on = "cohort"]
intervals[, death_before_next := !is.na(death_time) & death_time > interview_time &
  death_time <= study_end_time & (is.na(dest_time_living) | death_time < dest_time_living)]
intervals[, destination := fifelse(death_before_next, "DEAD", dest_state_living)]
intervals[, destination_time := fifelse(death_before_next, death_time, dest_time_living)]
intervals[, destination_wave := fifelse(death_before_next, NA_integer_, dest_wave_living)]
intervals[, interval_years := destination_time - interview_time]
intervals[, origin_state := history_state]
intervals[, origin_wave := wave]
intervals[, origin_age := age]
intervals[, origin_year := nominal_year]
intervals[, long_mhas_gap := !is.na(destination_wave) & cohort == "MHAS" & origin_wave == 2L & destination_wave == 3L]
intervals[, primary_interval := !is.na(destination) & !is.na(interval_years) &
  interval_years >= 1 & interval_years <= 4 & !long_mhas_gap &
  !is.na(origin_age) & origin_age >= 60 & origin_age <= 95 & !living_after_death]
intervals[, interval_exclusion := fifelse(is.na(destination), "no_observed_destination",
  fifelse(is.na(interval_years), "missing_interval",
    fifelse(interval_years <= 0, "nonpositive_interval",
      fifelse(interval_years < 1, "interval_lt_1y",
        fifelse(interval_years > 4, "interval_gt_4y",
          fifelse(long_mhas_gap, "prespecified_mhas_2003_2012_gap",
            fifelse(is.na(origin_age), "missing_origin_age",
            fifelse(origin_age < 60, "origin_age_lt_60",
              fifelse(origin_age > 95, "origin_age_gt_95",
                fifelse(living_after_death, "living_after_recorded_death", "included"))))))))))]

allowed_pairs <- c(
  "I0->I0", "I0->D1", "I0->DEAD",
  "D1->D1", "D1->R1", "D1->DEAD",
  "R1->R1", "R1->D2", "R1->DEAD",
  "D2->D2", "D2->R1", "D2->DEAD"
)
intervals[, transition := paste0(origin_state, "->", destination)]
intervals[, impossible_transition := !is.na(destination) & !transition %in% allowed_pairs]

# Scheduled-wave attrition audit keeps death separate from non-death missingness.
scheduled <- panel[, .(
  cohort, country, person_id, wave, interview_time, age, history_state, observed_state,
  iwstat, death_time
)]
scheduled[, origin_wave := wave]
scheduled[, next_wave := wave + 1L]
next_panel <- panel[, .(
  cohort, person_id, wave, next_observed_state = observed_state,
  next_history_state = history_state, next_iwstat = iwstat,
  next_interview_time = interview_time, next_in_wave = in_wave
)]
scheduled <- next_panel[scheduled, on = .(cohort, person_id, wave = next_wave)]
scheduled <- scheduled[!is.na(next_interview_time) & observed_state == TRUE & age >= 60 & age <= 95]
scheduled[, destination_wave := wave]
scheduled[, scheduled_outcome := fifelse(next_observed_state == TRUE, "observed_function",
  fifelse((!is.na(next_iwstat) & next_iwstat == 5) |
      (!is.na(death_time) & death_time > interview_time & death_time <= next_interview_time + 0.75), "recorded_death",
    fifelse(!is.na(next_iwstat) & next_iwstat != 5, "known_non_death_nonresponse_or_missing_function", "unknown_vital_or_observation_status")))]

# Aggregate audits.
wave_state <- panel[, .(
  source_rows = .N, in_wave_n = sum(in_wave == 1, na.rm = TRUE),
  age60plus_n = sum(in_wave == 1 & age >= 60, na.rm = TRUE),
  complete9_n = sum(in_wave == 1 & function_complete9, na.rm = TRUE),
  independent_n = sum(in_wave == 1 & difficulty == 0, na.rm = TRUE),
  difficulty_n = sum(in_wave == 1 & difficulty == 1, na.rm = TRUE),
  I0_n = sum(history_state == "I0", na.rm = TRUE), D1_n = sum(history_state == "D1", na.rm = TRUE),
  R1_n = sum(history_state == "R1", na.rm = TRUE), D2_n = sum(history_state == "D2", na.rm = TRUE),
  iwstat_death_n = sum(iwstat == 5, na.rm = TRUE),
  weight_positive_n = sum(respondent_weight > 0, na.rm = TRUE),
  proxy_observed_n = sum(!is.na(proxy)), proxy_n = sum(proxy == 1, na.rm = TRUE),
  wealth_observed_n = sum(in_wave == 1 & !is.na(wealth), na.rm = TRUE),
  education_observed_n = sum(in_wave == 1 & !is.na(education3), na.rm = TRUE),
  exact_interview_time_n = sum(in_wave == 1 & time_source == "interview_year_month", na.rm = TRUE)
), by = .(cohort, country, wave, nominal_year)]

missing_vars <- c("age", "education3", "wealth", "proxy", "respondent_weight", all_tokens)
missingness <- rbindlist(lapply(missing_vars, function(v) {
  panel[in_wave == 1, .(n = .N, missing_n = sum(is.na(get(v))), missing_percent = 100 * mean(is.na(get(v)))),
    by = .(cohort, wave)][, variable := v]
}), use.names = TRUE, fill = TRUE)
setcolorder(missingness, c("cohort", "wave", "variable", "n", "missing_n", "missing_percent"))

interval_audit <- intervals[, .(
  intervals = .N, persons = uniqueN(person_id),
  min_years = suppressWarnings(min(interval_years, na.rm = TRUE)),
  median_years = suppressWarnings(median(interval_years, na.rm = TRUE)),
  max_years = suppressWarnings(max(interval_years, na.rm = TRUE))
), by = .(cohort, interval_exclusion)]

transition_matrix <- intervals[primary_interval == TRUE & impossible_transition == FALSE,
  .(events = .N, persons = uniqueN(person_id)),
  by = .(cohort, origin_state, destination, transition)]

event_gate <- intervals[primary_interval == TRUE & impossible_transition == FALSE,
  .(
    incident_events = sum(origin_state == "I0" & destination == "D1"),
    incident_persons = uniqueN(person_id[origin_state == "I0" & destination == "D1"]),
    recovery_events = sum(origin_state %in% c("D1", "D2") & destination == "R1"),
    recovery_persons = uniqueN(person_id[origin_state %in% c("D1", "D2") & destination == "R1"]),
    relapse_events = sum(origin_state == "R1" & destination == "D2"),
    relapse_persons = uniqueN(person_id[origin_state == "R1" & destination == "D2"]),
    disabled_death_events = sum(origin_state %in% c("D1", "D2") & destination == "DEAD"),
    disabled_death_persons = uniqueN(person_id[origin_state %in% c("D1", "D2") & destination == "DEAD"]),
    total_death_events = sum(destination == "DEAD"),
    total_death_persons = uniqueN(person_id[destination == "DEAD"])
  ), by = .(cohort, country)]

event_by_interval <- intervals[primary_interval == TRUE & impossible_transition == FALSE,
  .(events = .N, persons = uniqueN(person_id)),
  by = .(cohort, origin_wave, destination_wave, origin_year, transition)]

attrition_audit <- scheduled[, .(n = .N, persons = uniqueN(person_id)),
  by = .(cohort, origin_wave, destination_wave, scheduled_outcome)]

iwstat_values <- panel[, .(n = .N), by = .(cohort, wave, iwstat)][order(cohort, wave, iwstat)]
death_source_audit <- unique(panel[, .(
  cohort, person_id, death_time_raw, death_source_raw, death_time, death_source,
  death_conflict_raw, nominal_time_boundary_adjustment, death_time_qc
)])[, 
  .(persons = .N, death_time_observed_n = sum(!is.na(death_time)),
    raw_death_time_n = sum(!is.na(death_time_raw)),
    raw_conflict_persons = sum(death_conflict_raw, na.rm = TRUE),
    bounded_nominal_persons = sum(nominal_time_boundary_adjustment, na.rm = TRUE),
    bounded_coarse_death_year_persons = sum(death_time_qc == "bounded_coarse_death_year", na.rm = TRUE),
    unresolved_conflict_persons = sum(death_time_qc == "unresolved_exact_date_conflict", na.rm = TRUE),
    static_death_n = sum(death_source_raw == "static_death_year_month" & !is.na(death_time), na.rm = TRUE),
    elsa_wave11_eol_death_n = sum(death_source_raw == "elsa_wave11_eol_year" & !is.na(death_time), na.rm = TRUE),
    iwstat_only_death_n = sum(death_source == "iwstat_5", na.rm = TRUE),
    min_death_time = suppressWarnings(min(death_time, na.rm = TRUE)),
    max_death_time = suppressWarnings(max(death_time, na.rm = TRUE))),
  by = cohort]

impossible <- intervals[, .(
  all_intervals = .N,
  impossible_transition_n = sum(impossible_transition, na.rm = TRUE),
  nonpositive_interval_n = sum(!is.na(interval_years) & interval_years <= 0),
  living_after_recorded_death_n = sum(living_after_death, na.rm = TRUE),
  raw_death_time_conflict_n = sum(death_conflict_raw, na.rm = TRUE)
), by = cohort]

schema <- rbindlist(schema_parts, fill = TRUE)
source_audit <- rbindlist(source_parts, fill = TRUE)
value_counts <- rbindlist(value_count_parts, fill = TRUE)
value_labels <- unique(rbindlist(value_label_parts, fill = TRUE))

saveRDS(panel, file.path(derived_dir, "individual_wave_panel_gate0.rds"), compress = "xz")
saveRDS(intervals, file.path(derived_dir, "transition_intervals_gate0.rds"), compress = "xz")
fwrite(source_audit, file.path(out_dir, "source_read_audit.csv"))
fwrite(schema, file.path(out_dir, "exact_schema_availability.csv"))
fwrite(value_counts, file.path(out_dir, "semantic_value_counts.csv"))
fwrite(value_labels, file.path(out_dir, "semantic_value_labels.csv"))
fwrite(wave_state, file.path(out_dir, "wave_state_counts.csv"))
fwrite(missingness, file.path(out_dir, "wave_variable_missingness.csv"))
fwrite(interval_audit, file.path(out_dir, "interval_length_audit.csv"))
fwrite(transition_matrix, file.path(out_dir, "transition_matrix_counts.csv"))
fwrite(event_gate, file.path(out_dir, "gate0_event_counts_by_cohort.csv"))
fwrite(event_by_interval, file.path(out_dir, "gate0_event_counts_by_interval.csv"))
fwrite(attrition_audit, file.path(out_dir, "death_non_death_attrition_audit.csv"))
fwrite(iwstat_values, file.path(out_dir, "iwstat_value_counts.csv"))
fwrite(death_source_audit, file.path(out_dir, "death_source_audit.csv"))
fwrite(impossible, file.path(out_dir, "impossible_transition_audit.csv"))

cat("\nGate-0 event counts:\n")
print(event_gate)
cat("\nImpossible-transition audit:\n")
print(impossible)
cat("\nDeath-source audit:\n")
print(death_source_audit)
cat("Completed: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_gate0.txt"))
