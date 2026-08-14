#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
out_dir <- file.path(root, "03_outputs", "00_source_audit")
log_dir <- file.path(root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "America/Chicago")
log_path <- file.path(log_dir, paste0("00_audit_source_schema_", stamp, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

sources <- data.table(
  cohort = c("CHARLS", "HRS", "HRS_GATEWAY", "ELSA", "MHAS"),
  country = c("China", "United States", "United States", "England", "Mexico"),
  path = c(
    Sys.getenv("RECOVERY_DIVIDE_CHARLS_FILE", unset = ""),
    Sys.getenv("RECOVERY_DIVIDE_HRS_FILE", unset = ""),
    Sys.getenv("RECOVERY_DIVIDE_HRS_HARMONIZED_FILE", unset = ""),
    Sys.getenv("RECOVERY_DIVIDE_ELSA_FILE", unset = ""),
    Sys.getenv("RECOVERY_DIVIDE_MHAS_FILE", unset = "")
  ),
  role = c("core", "core", "auxiliary", "core", "core")
)

concept_patterns <- list(
  person_id = "^(id|idauniq|hhidpn|rahhidnp)$",
  in_wave = "^inw[0-9]+$",
  interview_status = "^r[0-9]+iwstat$|iwstat|interview status",
  interview_time = "^r[0-9]+iw(year|month|day)$|interview (year|month|date)",
  age = "^r[0-9]+agey(_b)?$|age in years",
  sex = "^ragender$|gender of respondent|sex of respondent",
  education = "^raeduc(l)?$|education level|educational attainment|years of education",
  adl_count = "^r[0-9]+(adl5a|adlfivea?|adla)$|activities of daily living.*(count|summary)|number of adl",
  iadl_count = "^r[0-9]+(iadl5a|iadlfoura|iadlzaa?|iadla)$|instrumental activities of daily living.*(count|summary)|number of iadl",
  adl_item = "^r[0-9]+(dressa?|bathe?a?|eata?|beda?|toilt?a?|walkr?a?)$|difficulty.*(dress|bath|eat|bed|toilet|walk.*room)",
  iadl_item = "^r[0-9]+(phonea?|moneya?|medsa?|shopa?|mealsa?)$|difficulty.*(phone|money|medication|shopping|prepar.*meal)",
  proxy = "^r[0-9]+proxy$|proxy interview|proxy respondent",
  respondent_weight = "^r[0-9]+(l?wtresp|wtrespb)$|respondent weight|individual weight",
  survey_design = "^r[0-9]+(seclust|seclstr|psu|strat).*$|sampling (unit|stratum|cluster)",
  household_wealth = "^h[0-9]+atot[wb]$|total household wealth|net wealth|total wealth",
  death = "^ra(dyear|dmonth|dage|died).*$|death|died|deceased|vital status",
  household_size = "^h[0-9]+hhres$|household size|household residents"
)

safe_label <- function(x) {
  z <- attr(x, "label", exact = TRUE)
  if (is.null(z) || !length(z)) NA_character_ else as.character(z[[1]])
}

dictionary_parts <- list()
source_rows <- list()

cat("Started: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")
cat("R: ", R.version.string, "\n", sep = "")

for (i in seq_len(nrow(sources))) {
  s <- sources[i]
  cat("Reading metadata: ", s$cohort, "\n", sep = "")
  if (!file.exists(s$path)) {
    source_rows[[length(source_rows) + 1L]] <- data.table(
      cohort = s$cohort, path = s$path, exists = FALSE, bytes = NA_real_,
      rows = NA_integer_, columns = NA_integer_, read_status = "MISSING"
    )
    next
  }

  meta <- read_dta(s$path, n_max = 0)
  vars <- names(meta)
  labels <- vapply(meta, safe_label, character(1))
  haystack <- tolower(paste(vars, labels))

  src_hits <- list()
  for (concept in names(concept_patterns)) {
    pattern <- concept_patterns[[concept]]
    hit <- grepl(pattern, vars, ignore.case = TRUE, perl = TRUE) |
      grepl(pattern, labels, ignore.case = TRUE, perl = TRUE)
    if (any(hit)) {
      src_hits[[length(src_hits) + 1L]] <- data.table(
        cohort = s$cohort,
        role = s$role,
        concept = concept,
        variable = vars[hit],
        label = labels[hit]
      )
    }
  }
  if (length(src_hits)) dictionary_parts[[length(dictionary_parts) + 1L]] <- rbindlist(src_hits, fill = TRUE)

  # The row count is read from the ID only so the audit never imports the full wide file.
  id_candidates <- vars[grepl(concept_patterns$person_id, vars, ignore.case = TRUE, perl = TRUE)]
  n_rows <- NA_integer_
  read_status <- "METADATA_ONLY_NO_ID"
  if (length(id_candidates)) {
    id_only <- read_dta(s$path, col_select = all_of(id_candidates[[1]]))
    n_rows <- nrow(id_only)
    read_status <- "PASS"
    rm(id_only)
  }
  source_rows[[length(source_rows) + 1L]] <- data.table(
    cohort = s$cohort, path = s$path, exists = TRUE,
    bytes = as.numeric(file.info(s$path)$size), rows = n_rows,
    columns = length(vars), read_status = read_status
  )
  rm(meta)
  invisible(gc())
}

dictionary <- unique(rbindlist(dictionary_parts, fill = TRUE))
setorder(dictionary, cohort, concept, variable)
source_audit <- rbindlist(source_rows, fill = TRUE)
setorder(source_audit, cohort)

fwrite(dictionary, file.path(out_dir, "source_candidate_dictionary.csv"))
fwrite(source_audit, file.path(out_dir, "source_file_audit.csv"))

core_concepts <- c(
  "person_id", "in_wave", "interview_status", "interview_time", "age", "sex",
  "education", "adl_count", "iadl_count", "adl_item", "iadl_item", "proxy",
  "respondent_weight", "household_wealth", "death"
)
core_cohorts <- c("CHARLS", "HRS", "ELSA", "MHAS")
coverage <- CJ(cohort = core_cohorts, concept = core_concepts)
coverage <- dictionary[cohort %in% core_cohorts,
  .(n_variables = uniqueN(variable), variables = paste(sort(unique(variable)), collapse = ";")),
  by = .(cohort, concept)
][coverage, on = .(cohort, concept)]
coverage[is.na(n_variables), `:=`(n_variables = 0L, variables = "")]
coverage[, available := n_variables > 0]
coverage[, concept_order := match(concept, core_concepts)]
setorder(coverage, cohort, concept_order)
coverage[, concept_order := NULL]
fwrite(coverage, file.path(out_dir, "source_concept_coverage.csv"))

cat("\nSource audit:\n")
print(source_audit)
cat("\nCore concept coverage:\n")
print(coverage[, .(cohort, concept, available, n_variables)])
cat("Completed: ", format(Sys.time(), tz = "America/Chicago"), "\n", sep = "")

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_schema_audit.txt"))
