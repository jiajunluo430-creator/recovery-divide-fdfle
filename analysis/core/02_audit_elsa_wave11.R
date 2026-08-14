#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(haven)
})

project_root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
stage_root <- file.path(
  project_root,
  "00_staging/elsa_wave11/UKDA-5050-stata/stata/stata13_se"
)
out_dir <- file.path(project_root, "03_outputs/02_elsa_wave11_audit")
log_dir <- file.path(project_root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

core_path <- file.path(stage_root, "wave_11_elsa_data_eul.dta")
eol_path <- file.path(stage_root, "wave_11_elsa_eol_eul.dta")
stopifnot(file.exists(core_path), file.exists(eol_path))

label_string <- function(x) {
  lab <- attr(x, "label", exact = TRUE)
  if (is.null(lab)) "" else as.character(lab)
}

value_label_string <- function(x) {
  labs <- attr(x, "labels", exact = TRUE)
  if (is.null(labs) || length(labs) == 0L) return("")
  paste0(names(labs), "=", unname(labs), collapse = " | ")
}

metadata_frame <- function(dat, file_role) {
  data.frame(
    file_role = file_role,
    variable = names(dat),
    storage_class = vapply(dat, function(x) paste(class(x), collapse = "/"), character(1)),
    variable_label = vapply(dat, label_string, character(1)),
    value_labels = vapply(dat, value_label_string, character(1)),
    nonmissing_n = vapply(dat, function(x) sum(!is.na(x)), integer(1)),
    unique_nonmissing_n = vapply(dat, function(x) length(unique(x[!is.na(x)])), integer(1)),
    stringsAsFactors = FALSE
  )
}

candidate_pattern <- paste(
  c(
    "idauniq", "death", "died", "deceased", "end.?of.?life", "date", "year", "month",
    "interview", "proxy", "respondent", "weight", "nurse", "age", "sex", "wave",
    "adl", "iadl", "difficulty", "status", "outcome", "mortality"
  ),
  collapse = "|"
)

message("Reading wave 11 core metadata and selected data...")
core_meta_only <- read_dta(core_path, n_max = 1L)
core_names <- names(core_meta_only)
core_labels <- vapply(core_meta_only, label_string, character(1))
core_keep <- unique(c(
  "idauniq",
  grep(candidate_pattern, core_names, value = TRUE, ignore.case = TRUE),
  core_names[grepl(candidate_pattern, paste(core_names, core_labels), ignore.case = TRUE, perl = TRUE)],
  grep("^headl(dr|wa|ba|ea|be|pr|sh|ph|me|mo)$", core_names, value = TRUE, ignore.case = TRUE)
))
core_keep <- intersect(core_keep, core_names)
core <- read_dta(core_path, col_select = all_of(core_keep))

message("Reading wave 11 EOL data...")
eol <- read_dta(eol_path)

core_meta <- metadata_frame(core, "wave11_core_selected")
eol_meta <- metadata_frame(eol, "wave11_eol_all")
write.csv(core_meta, file.path(out_dir, "core_selected_variable_metadata.csv"), row.names = FALSE, na = "")
write.csv(eol_meta, file.path(out_dir, "eol_all_variable_metadata.csv"), row.names = FALSE, na = "")

combined <- rbind(core_meta, eol_meta)
hit <- grepl(
  candidate_pattern,
  paste(combined$variable, combined$variable_label),
  ignore.case = TRUE,
  perl = TRUE
)
write.csv(
  combined[hit, , drop = FALSE],
  file.path(out_dir, "death_date_proxy_weight_candidate_metadata.csv"),
  row.names = FALSE,
  na = ""
)

strict_items <- c(
  "headldr", "headlwc", "headlwa", "headlba", "headlea", "headlbe",
  "headlpr", "headlsh", "headlme", "headlmo"
)
missing_items <- setdiff(strict_items, names(core))
if (length(missing_items)) {
  stop("Missing frozen strict functional items in ELSA wave 11: ", paste(missing_items, collapse = ", "))
}

item_valid <- lapply(core[strict_items], function(x) ifelse(as.numeric(x) %in% c(0, 1), as.numeric(x), NA_real_))
item_mat <- do.call(cbind, item_valid)
complete_9 <- rowSums(!is.na(item_mat)) == length(strict_items)
strict_difficulty <- ifelse(complete_9, as.integer(rowSums(item_mat) >= 1), NA_integer_)

function_summary <- data.frame(
  file_role = "wave11_core",
  row_n = nrow(core),
  unique_id_n = if ("idauniq" %in% names(core)) length(unique(core$idauniq[!is.na(core$idauniq)])) else NA_integer_,
  complete_strict_9_n = sum(complete_9),
  strict_independent_n = sum(strict_difficulty == 0, na.rm = TRUE),
  strict_difficulty_n = sum(strict_difficulty == 1, na.rm = TRUE),
  strict_missing_n = sum(is.na(strict_difficulty)),
  phone_difficulty_positive_n = if ("headlph" %in% names(core)) sum(as.numeric(core$headlph) == 1, na.rm = TRUE) else NA_integer_,
  stringsAsFactors = FALSE
)
write.csv(function_summary, file.path(out_dir, "wave11_strict_function_summary.csv"), row.names = FALSE, na = "")

id_audit <- data.frame(
  file_role = c("wave11_core", "wave11_eol"),
  rows = c(nrow(core), nrow(eol)),
  idauniq_present = c("idauniq" %in% names(core), "idauniq" %in% names(eol)),
  unique_id_n = c(
    if ("idauniq" %in% names(core)) length(unique(core$idauniq[!is.na(core$idauniq)])) else NA_integer_,
    if ("idauniq" %in% names(eol)) length(unique(eol$idauniq[!is.na(eol$idauniq)])) else NA_integer_
  ),
  duplicate_id_rows_n = c(
    if ("idauniq" %in% names(core)) sum(duplicated(core$idauniq) & !is.na(core$idauniq)) else NA_integer_,
    if ("idauniq" %in% names(eol)) sum(duplicated(eol$idauniq) & !is.na(eol$idauniq)) else NA_integer_
  ),
  stringsAsFactors = FALSE
)
write.csv(id_audit, file.path(out_dir, "wave11_id_audit.csv"), row.names = FALSE, na = "")

if ("idauniq" %in% names(core) && "idauniq" %in% names(eol)) {
  overlap <- merge(
    unique(data.frame(idauniq = core$idauniq[!is.na(core$idauniq)])),
    unique(data.frame(idauniq = eol$idauniq[!is.na(eol$idauniq)])),
    by = "idauniq"
  )
  write.csv(
    data.frame(core_eol_id_overlap_n = nrow(overlap)),
    file.path(out_dir, "wave11_core_eol_overlap.csv"),
    row.names = FALSE
  )
}

sink(file.path(log_dir, "02_audit_elsa_wave11_sessionInfo.txt"))
print(sessionInfo())
sink()

message("ELSA wave 11 audit complete: ", out_dir)
