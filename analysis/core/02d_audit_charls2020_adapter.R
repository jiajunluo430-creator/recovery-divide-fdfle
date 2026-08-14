#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
source_root <- Sys.getenv("RECOVERY_DIVIDE_CHARLS_RAW_ROOT", unset = "")
out_dir <- file.path(root, "03_outputs", "02d_charls2020_adapter_audit")
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  harmonized_prior = Sys.getenv("RECOVERY_DIVIDE_CHARLS_FILE", unset = ""),
  working_wave5 = file.path(source_root, "Working_data", "charls_wave5.dta"),
  temp_charls20 = file.path(source_root, "Temp_data", "charls20.dta"),
  raw_health = file.path(source_root, "Raw_data", "2020charls", "Health_Status_and_Functioning.dta"),
  raw_sample = file.path(source_root, "Raw_data", "2020charls", "Sample_Infor.dta"),
  raw_weights = file.path(source_root, "Raw_data", "2020charls", "Weights.dta"),
  raw_exit = file.path(source_root, "Raw_data", "2020charls", "Exit_Module.dta"),
  translated_combined_w5 = file.path(source_root, "Raw_data", "汉化版CHARLS", "w5_charls_CN.dta")
)
stopifnot(all(file.exists(paths)))

chr_id <- function(z) {
  if (is.character(z)) return(trimws(z))
  if (is.numeric(z)) return(format(z, scientific = FALSE, trim = TRUE, digits = 22))
  trimws(as.character(z))
}
num <- function(z) suppressWarnings(as.numeric(z))
valid_year <- function(z) fifelse(num(z) >= 1900 & num(z) <= 2030, num(z), NA_real_)
valid_month <- function(z) fifelse(num(z) %in% 1:12, num(z), NA_real_)
raw_difficulty <- function(z) fifelse(num(z) == 1, 0, fifelse(num(z) %in% 2:4, 1, NA_real_))
valid01 <- function(z) fifelse(num(z) %in% 0:1, num(z), NA_real_)

tokens <- c("dressa", "batha", "eata", "beda", "toilta", "shopa", "mealsa", "medsa", "moneya")
raw_names <- c(
  dressa = "db001", batha = "db003", eata = "db005", beda = "db007", toilta = "db009",
  shopa = "db016", mealsa = "db014", medsa = "db020", moneya = "db022"
)
temp_names <- setNames(paste0("r5", tokens), tokens)

read_selected <- function(path, vars) {
  meta <- read_dta(path, n_max = 0)
  keep <- intersect(vars, names(meta))
  as.data.table(read_dta(path, col_select = all_of(keep)))
}

prior <- read_selected(paths[["harmonized_prior"]], c("ID", "inw4", "r4iwstat", "r4agey", "raeducl"))
working <- read_selected(paths[["working_wave5"]], c(
  "ID", "wave", "iwy", "iwm", "age", "iwstat", "raeducl", tokens
))
temp <- read_selected(paths[["temp_charls20"]], c(
  "ID", "inw5", "died", "r5iwy", "r5iwm", "exb001_1", "exb001_2", "exb001_3",
  "r5wtresp", "r5wtrespb", "r5gender", "zrbirthyear", "r5birthmonth", "r5educ_c", unname(temp_names)
))
health <- read_selected(paths[["raw_health"]], c("ID", "proxy", "proxy_5", unname(raw_names)))
sample <- read_selected(paths[["raw_sample"]], c("ID", "died", "iyear", "imonth"))
weights <- read_selected(paths[["raw_weights"]], c("ID", "INDV_weight", "INDV_weight_ad2"))
exit <- read_selected(paths[["raw_exit"]], c("ID", "exb001_1", "exb001_2", "exb001_3", "exb002"))
combined <- read_selected(paths[["translated_combined_w5"]], c(
  "ID", "iyear", "imonth", "xiwyear", "xiwmonth", "proxy", "proxy_5", "INDV_weight", "INDV_weight_ad2"
))

objects <- list(prior = prior, working = working, temp = temp, health = health, sample = sample, weights = weights, exit = exit, combined = combined)
file_audit <- rbindlist(lapply(names(objects), function(role) {
  z <- objects[[role]]
  z[, person_id_audit := chr_id(ID)]
  out <- data.table(
    source_role = role, rows = nrow(z), unique_ids = uniqueN(z$person_id_audit),
    missing_ids = sum(is.na(z$person_id_audit) | !nzchar(z$person_id_audit)),
    duplicated_id_rows = sum(duplicated(z$person_id_audit) | duplicated(z$person_id_audit, fromLast = TRUE))
  )
  z[, person_id_audit := NULL]
  out
}))

for (z in objects) z[, person_id := chr_id(ID)]
prior_ids <- unique(prior$person_id)
role_ids <- rbindlist(lapply(names(objects)[-1], function(role) {
  ids <- unique(objects[[role]]$person_id)
  data.table(
    source_role = role,
    ids = length(ids),
    overlap_prior_ids = sum(ids %in% prior_ids),
    percent_source_ids_in_prior = 100 * mean(ids %in% prior_ids),
    prior_ids_captured = sum(prior_ids %in% ids),
    percent_prior_ids_captured = 100 * mean(prior_ids %in% ids)
  )
}))

function_value_parts <- list()
agreement_parts <- list()
for (token in tokens) {
  raw_var <- raw_names[[token]]
  temp_var <- temp_names[[token]]
  z_values <- health[, .N, by = .(raw_value = as.character(get(raw_var)))]
  z_values[, `:=`(source = "raw_health", item = token, source_variable = raw_var)]
  function_value_parts[[length(function_value_parts) + 1L]] <- z_values
  z_values <- temp[, .N, by = .(raw_value = as.character(get(temp_var)))]
  z_values[, `:=`(source = "temp_charls20", item = token, source_variable = temp_var)]
  function_value_parts[[length(function_value_parts) + 1L]] <- z_values
  z_values <- working[, .N, by = .(raw_value = as.character(get(token)))]
  z_values[, `:=`(source = "working_wave5", item = token, source_variable = token)]
  function_value_parts[[length(function_value_parts) + 1L]] <- z_values

  compare <- merge(
    health[, .(person_id, raw_recoded = raw_difficulty(get(raw_var)))],
    temp[, .(person_id, temp_recoded = valid01(get(temp_var)))],
    by = "person_id", all = FALSE
  )
  compare <- merge(compare, working[, .(person_id, working_recoded = valid01(get(token)))], by = "person_id", all = TRUE)
  agreement_parts[[length(agreement_parts) + 1L]] <- data.table(
    item = token,
    raw_temp_overlap = sum(!is.na(compare$raw_recoded) & !is.na(compare$temp_recoded)),
    raw_temp_agree = sum(compare$raw_recoded == compare$temp_recoded, na.rm = TRUE),
    raw_temp_agreement_percent = 100 * mean(compare$raw_recoded == compare$temp_recoded,
      na.rm = TRUE),
    raw_working_overlap = sum(!is.na(compare$raw_recoded) & !is.na(compare$working_recoded)),
    raw_working_agree = sum(compare$raw_recoded == compare$working_recoded, na.rm = TRUE),
    raw_working_agreement_percent = 100 * mean(compare$raw_recoded == compare$working_recoded,
      na.rm = TRUE)
  )
}
function_values <- rbindlist(function_value_parts)
function_agreement <- rbindlist(agreement_parts)

health_state <- health[, .(person_id)]
for (token in tokens) health_state[[token]] <- raw_difficulty(health[[raw_names[[token]]]])
health_state[, complete9 := rowSums(!is.na(.SD)) == length(tokens), .SDcols = tokens]
health_state[, any_difficulty9 := rowSums(.SD == 1, na.rm = TRUE) > 0, .SDcols = tokens]
canonical_status <- merge(
  sample[, .(
    person_id, died = valid01(died), interview_year = valid_year(iyear), interview_month = valid_month(imonth)
  )],
  health_state[, .(person_id, complete9, any_difficulty9)],
  by = "person_id", all.x = TRUE
)
canonical_status <- merge(
  canonical_status,
  weights[, .(person_id, base_weight = num(INDV_weight), adjusted_weight = num(INDV_weight_ad2))],
  by = "person_id", all.x = TRUE
)
canonical_status <- canonical_status[, .(
  rows = .N,
  unique_ids = uniqueN(person_id),
  complete9 = sum(complete9, na.rm = TRUE),
  difficulty9 = sum(complete9 & any_difficulty9, na.rm = TRUE),
  independent9 = sum(complete9 & !any_difficulty9, na.rm = TRUE),
  valid_interview_date = sum(!is.na(interview_year) & !is.na(interview_month)),
  positive_base_weight = sum(is.finite(base_weight) & base_weight > 0),
  positive_adjusted_weight = sum(is.finite(adjusted_weight) & adjusted_weight > 0)
), by = .(
  inw5 = as.character(fifelse(died == 0, 1, 0)),
  died = as.character(died)
)]

proxy_audit <- health[, .(
  rows = .N,
  generic_proxy_yes = sum(num(proxy) == 1, na.rm = TRUE),
  generic_proxy_skip_self = sum(is.na(num(proxy))),
  proxy5_yes = sum(num(proxy_5) == 1, na.rm = TRUE),
  proxy5_no = sum(num(proxy_5) == 2, na.rm = TRUE),
  canonical_proxy_observed = .N,
  canonical_proxy = sum(num(proxy) == 1, na.rm = TRUE),
  canonical_self = sum(is.na(num(proxy)))
)]

weight_compare <- merge(
  temp[, .(person_id, temp_base = num(r5wtresp), temp_adjusted = num(r5wtrespb))],
  weights[, .(person_id, raw_base = num(INDV_weight), raw_adjusted = num(INDV_weight_ad2))],
  by = "person_id", all = FALSE
)
weight_audit <- data.table(
  overlap_ids = nrow(weight_compare),
  base_both_nonmissing = sum(is.finite(weight_compare$temp_base) & is.finite(weight_compare$raw_base)),
  base_exact_agree = sum(weight_compare$temp_base == weight_compare$raw_base, na.rm = TRUE),
  base_max_abs_difference = max(abs(weight_compare$temp_base - weight_compare$raw_base), na.rm = TRUE),
  adjusted_both_nonmissing = sum(is.finite(weight_compare$temp_adjusted) & is.finite(weight_compare$raw_adjusted)),
  adjusted_exact_agree = sum(weight_compare$temp_adjusted == weight_compare$raw_adjusted, na.rm = TRUE),
  adjusted_max_abs_difference = max(abs(weight_compare$temp_adjusted - weight_compare$raw_adjusted), na.rm = TRUE)
)

date_compare <- merge(
  temp[, .(person_id, temp_year = valid_year(r5iwy), temp_month = valid_month(r5iwm))],
  sample[, .(person_id, raw_year = valid_year(iyear), raw_month = valid_month(imonth))],
  by = "person_id", all = FALSE
)
date_audit <- data.table(
  overlap_ids = nrow(date_compare),
  year_both_nonmissing = sum(!is.na(date_compare$temp_year) & !is.na(date_compare$raw_year)),
  year_agree = sum(date_compare$temp_year == date_compare$raw_year, na.rm = TRUE),
  month_both_nonmissing = sum(!is.na(date_compare$temp_month) & !is.na(date_compare$raw_month)),
  month_agree = sum(date_compare$temp_month == date_compare$raw_month, na.rm = TRUE)
)

exit_compare <- merge(
  exit[, .(
    person_id, exit_year = valid_year(exb001_1), exit_month = valid_month(exb001_2), exit_day = num(exb001_3)
  )],
  temp[, .(
    person_id, died = num(died), temp_death_year = valid_year(exb001_1),
    temp_death_month = valid_month(exb001_2), temp_death_day = num(exb001_3)
  )],
  by = "person_id", all = TRUE
)
death_audit <- data.table(
  raw_exit_ids = uniqueN(exit$person_id),
  sample_confirmed_died_ids = uniqueN(sample[died == 1, person_id]),
  temp_died_ids = uniqueN(temp[died == 1, person_id]),
  exit_ids_marked_died_in_sample = sum(unique(exit$person_id) %in% sample[died == 1, person_id]),
  sample_died_ids_present_in_exit = sum(unique(sample[died == 1, person_id]) %in% exit$person_id),
  sample_temp_died_exact_agreement = sum(unique(sample[died == 1, person_id]) %in% temp[died == 1, person_id]),
  raw_exit_valid_year = sum(!is.na(valid_year(exit$exb001_1))),
  raw_exit_valid_month = sum(!is.na(valid_month(exit$exb001_2))),
  death_year_agreement = sum(exit_compare$exit_year == exit_compare$temp_death_year, na.rm = TRUE),
  death_month_agreement = sum(exit_compare$exit_month == exit_compare$temp_death_month, na.rm = TRUE)
)

fwrite(file_audit, file.path(out_dir, "charls2020_record_id_audit.csv"))
fwrite(role_ids, file.path(out_dir, "charls2020_prior_id_linkage_audit.csv"))
fwrite(function_values, file.path(out_dir, "charls2020_function_value_counts.csv"))
fwrite(function_agreement, file.path(out_dir, "charls2020_function_recode_agreement.csv"))
fwrite(canonical_status, file.path(out_dir, "charls2020_status_function_date_weight_audit.csv"))
fwrite(proxy_audit, file.path(out_dir, "charls2020_proxy_audit.csv"))
fwrite(weight_audit, file.path(out_dir, "charls2020_weight_source_agreement.csv"))
fwrite(date_audit, file.path(out_dir, "charls2020_interview_date_source_agreement.csv"))
fwrite(death_audit, file.path(out_dir, "charls2020_exit_death_agreement.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_charls2020_adapter_audit.txt"))

cat("CHARLS 2020 adapter audit complete\n")
print(file_audit)
print(role_ids)
print(function_agreement)
print(canonical_status)
print(proxy_audit)
print(weight_audit)
print(date_audit)
print(death_audit)
