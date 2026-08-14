#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
source_root <- Sys.getenv("RECOVERY_DIVIDE_CHARLS_RAW_ROOT", unset = "")
out_dir <- file.path(root, "03_outputs", "02c_charls2020_schema_audit")
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

sources <- c(
  working_wave5 = file.path(source_root, "Working_data", "charls_wave5.dta"),
  temp_charls2020 = file.path(source_root, "Temp_data", "charls2020.dta"),
  temp_charls20 = file.path(source_root, "Temp_data", "charls20.dta"),
  translated_combined_w5 = file.path(source_root, "Raw_data", "汉化版CHARLS", "w5_charls_CN.dta"),
  raw_health = file.path(source_root, "Raw_data", "2020charls", "Health_Status_and_Functioning.dta"),
  raw_sample = file.path(source_root, "Raw_data", "2020charls", "Sample_Infor.dta"),
  raw_weights = file.path(source_root, "Raw_data", "2020charls", "Weights.dta"),
  raw_exit = file.path(source_root, "Raw_data", "2020charls", "Exit_Module.dta"),
  raw_demographic = file.path(source_root, "Raw_data", "2020charls", "Demographic_Background.dta")
)
stopifnot(all(file.exists(sources)))

value_label_string <- function(z) {
  v <- attr(z, "labels")
  if (is.null(v) || !length(v)) return("")
  paste0(names(v), "=", unname(v), collapse = " | ")
}

all_parts <- list()
file_parts <- list()
for (role in names(sources)) {
  path <- sources[[role]]
  cat("Reading metadata: ", role, "\n", sep = "")
  x <- read_dta(path, n_max = 0)
  all_parts[[length(all_parts) + 1L]] <- data.table(
    source_role = role,
    source_path = path,
    variable = names(x),
    storage_class = vapply(x, function(z) paste(class(z), collapse = "/"), character(1)),
    variable_label = vapply(x, function(z) {
      v <- attr(z, "label")
      if (is.null(v)) "" else as.character(v)
    }, character(1)),
    value_labels = vapply(x, value_label_string, character(1))
  )
  info <- file.info(path)
  file_parts[[length(file_parts) + 1L]] <- data.table(
    source_role = role, source_path = path,
    file_size_bytes = info$size,
    modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    variables = length(x)
  )
  rm(x)
  invisible(gc())
}

dictionary <- rbindlist(all_parts)
files <- rbindlist(file_parts)
keyword <- paste(c(
  "(^|_)id($|_)", "hhid", "household", "person", "individual", "sample",
  "adl", "iadl", "dress", "bath", "shower", "eat", "feed", "bed", "transfer",
  "toilet", "lavator", "shop", "grocery", "meal", "cook", "medic", "money", "finance",
  "difficulty", "difficult", "help", "interview", "date", "year", "month", "wave",
  "weight", "proxy", "death", "deceas", "exit", "mortality", "status"
), collapse = "|")
dictionary[, search_text := paste(variable, variable_label)]
candidates <- dictionary[grepl(keyword, search_text, ignore.case = TRUE)]
candidates[, search_text := NULL]

fwrite(files, file.path(out_dir, "charls2020_candidate_file_audit.csv"))
fwrite(dictionary, file.path(out_dir, "charls2020_all_variable_metadata.csv"))
fwrite(candidates, file.path(out_dir, "charls2020_keyword_candidate_metadata.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_charls2020_schema_audit.txt"))

cat("CHARLS 2020 schema audit complete\n")
print(files)
cat("Candidate variables by source:\n")
print(candidates[, .N, by = source_role])
