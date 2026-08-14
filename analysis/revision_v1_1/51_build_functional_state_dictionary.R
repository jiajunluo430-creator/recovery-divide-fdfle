#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
out_dir <- file.path(revision_root, "03_outputs", "13_data_dictionary")
log_dir <- file.path(revision_root, "06_logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

dictionary_path <- file.path(root, "03_outputs", "00_source_audit", "source_candidate_dictionary.csv")
source_dictionary <- fread(dictionary_path)

item_lookup <- c(
  dressa = "Dressing",
  batha = "Bathing or showering",
  eata = "Eating",
  beda = "Getting in or out of bed",
  toilta = "Using the toilet",
  shopa = "Shopping for groceries",
  mealsa = "Preparing a hot meal",
  medsa = "Taking medications",
  moneya = "Managing money"
)
item_domain <- c(
  dressa = "ADL", batha = "ADL", eata = "ADL", beda = "ADL", toilta = "ADL",
  shopa = "IADL", mealsa = "IADL", medsa = "IADL", moneya = "IADL"
)

harmonized <- source_dictionary[concept %chin% c("adl_item", "iadl_item")]
harmonized[, canonical_token := sub("^r[0-9]+", "", variable)]
harmonized <- harmonized[canonical_token %chin% names(item_lookup)]
harmonized[, wave := as.integer(sub("^r([0-9]+).*$", "\\1", variable))]
harmonized[, `:=`(
  canonical_item = unname(item_lookup[canonical_token]),
  domain = unname(item_domain[canonical_token]),
  source_layer = "harmonized cohort file",
  independent_code = "0 (No difficulty)",
  difficulty_code = "1 (Any difficulty)",
  missing_rule = "Special missing/refusal/don't know/skip values are missing; state missing if any required item is unresolved unless a validated complete count yields the same classification"
)]
harmonized <- harmonized[, .(
  cohort, wave, domain, canonical_item, canonical_token,
  source_variable = variable, source_label = label, source_layer,
  independent_code, difficulty_code, missing_rule
)]

raw_adapters <- data.table(
  cohort = c(rep("CHARLS", 9L), rep("ELSA", 9L)),
  wave = c(rep(5L, 9L), rep(11L, 9L)),
  canonical_token = rep(names(item_lookup), 2L),
  source_variable = c(
    "db001", "db003", "db005", "db007", "db009", "db016", "db014", "db020", "db022",
    "headldr", "headlba", "headlea", "headlbe", "headlwc", "headlsh", "headlpr", "headlme", "headlmo"
  )
)
raw_adapters[, `:=`(
  domain = unname(item_domain[canonical_token]),
  canonical_item = unname(item_lookup[canonical_token]),
  source_label = fifelse(
    cohort == "CHARLS",
    "CHARLS 2020 raw Health Status and Functioning item",
    "ELSA wave 11 raw ADL/IADL item"
  ),
  source_layer = fifelse(cohort == "CHARLS", "audited raw 2020 adapter", "audited raw wave 11 adapter"),
  independent_code = fifelse(cohort == "CHARLS", "raw 1 recoded to 0", "raw 0 retained as 0"),
  difficulty_code = fifelse(cohort == "CHARLS", "raw 2-4 recoded to 1", "raw 1 retained as 1"),
  missing_rule = "All other or missing raw values recoded to missing; the same nine-item complete-state rule is applied"
)]
raw_adapters <- raw_adapters[, names(harmonized), with = FALSE]

functional_mapping <- rbindlist(list(harmonized, raw_adapters), use.names = TRUE, fill = TRUE)
functional_mapping[, domain_order := match(domain, c("ADL", "IADL"))]
setorder(functional_mapping, cohort, wave, domain_order, canonical_item)
functional_mapping[, domain_order := NULL]

state_dictionary <- data.table(
  state = c("I0", "D1", "R1", "D2", "DEAD", "CENSORED_OR_UNKNOWN"),
  functional_status = c("Independent", "Difficulty", "Independent", "Difficulty", "Death", "Not a functional state"),
  observed_history_definition = c(
    "No previously observed functional difficulty within the cohort window",
    "First or continuing observed difficulty before any observed recovery",
    "Observed independence after D1 or D2",
    "Observed difficulty after an R1 state",
    "Verified cohort-specific death; absorbing",
    "Non-death loss, known-alive nonresponse, or unknown vital/observation status"
  ),
  contributes_to = c("FDFLE", "Years with difficulty", "FDFLE", "Years with difficulty", "Neither", "Censoring only"),
  caveat = c(
    "I0 means no prior observed difficulty, not lifetime absence of difficulty",
    "Baseline-prevalent D1 has unknown pre-entry duration",
    "Interview-observed recovery; exact transition time and intervening paths are unknown",
    "Interview-observed relapse; exact transition time and intervening paths are unknown",
    "Disappearance is never coded as death",
    "Never treated as death; addressed through retention/IPCW sensitivity"
  )
)

transition_dictionary <- data.table(
  process = c("Onset", "Recovery", "Recovery after relapse", "Relapse", "Pre-difficulty death", "Post-difficulty death"),
  origin = c("I0", "D1", "D2", "R1", "I0", "D1/R1/D2"),
  destination = c("D1", "R1", "R1", "D2", "DEAD", "DEAD"),
  primary_model_block = c("onset", "recovery", "recovery", "relapse", "pre_difficulty_mortality", "post_difficulty_mortality"),
  model_note = c(
    "Origin-specific onset process",
    "D1 and D2 recovery share the SES coefficient; origin state remains adjusted",
    "D1 and D2 recovery share the SES coefficient; origin state remains adjusted",
    "Origin-specific relapse process",
    "Separate pre-difficulty mortality SES coefficient",
    "D1, R1, and D2 share the SES coefficient; origin state remains adjusted"
  )
)

fwrite(functional_mapping, file.path(out_dir, "functional_item_mapping_by_cohort_wave.csv"))
fwrite(state_dictionary, file.path(out_dir, "history_state_dictionary.csv"))
fwrite(transition_dictionary, file.path(out_dir, "allowed_transition_and_model_blocks.csv"))
fwrite(data.table(
  construct = c("Primary broad", "ADL-only sensitivity", "At-least-two sensitivity", "Confirmed-state sensitivity"),
  rule = c(
    "At least one difficulty among all five ADL and four IADL items",
    "At least one difficulty among the five ADL items",
    "At least two difficulties among the nine common ADL/IADL items",
    "A destination living state must persist at the next formal adjacent living interview"
  ),
  interpretation = c(
    "Reported functional difficulty, including mild and potentially reversible difficulty",
    "Difficulty in personal activities of daily living",
    "Broader multi-item functional burden",
    "Misclassification-resistant persistence check; not a redefinition of the primary outcome"
  )
), file.path(out_dir, "functional_construct_definitions.csv"))

cat("Rows in functional item mapping: ", nrow(functional_mapping), "\n", sep = "")
cat("Cohort-wave coverage:\n")
print(functional_mapping[, .(items = uniqueN(canonical_item)), by = .(cohort, wave)][order(cohort, wave)])
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_51_functional_state_dictionary.txt"))
