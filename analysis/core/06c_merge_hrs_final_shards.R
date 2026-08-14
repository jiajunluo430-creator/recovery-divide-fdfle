#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
final_dir <- file.path(root, "03_outputs", "05_bootstrap_final")
shard_dirs <- file.path(root, "03_outputs", c("05_bootstrap_final_hrs_shard_a", "05_bootstrap_final_hrs_shard_b"))
types <- c("replicates", "qc", "wealth_cutpoint_audit")

for (type in types) {
  paths <- file.path(shard_dirs, paste0("bootstrap_hrs_", type, ".csv"))
  if (any(!file.exists(paths))) stop("Missing shard files for ", type)
  z <- rbindlist(lapply(paths, fread), fill = TRUE)
  keys <- if (type == "replicates") {
    c("cohort", "phase", "replicate", "seed", "exposure")
  } else if (type == "qc") {
    c("cohort", "phase", "replicate", "seed", "exposure")
  } else {
    intersect(c("cohort", "replicate", "wealth_entry_wave", "q33", "q67", "cut_method"), names(z))
  }
  setorderv(z, keys)
  z <- unique(z, by = keys)
  fwrite(z, file.path(final_dir, paste0("bootstrap_hrs_", type, ".csv")))
}

qc <- fread(file.path(final_dir, "bootstrap_hrs_qc.csv"))
replicates <- fread(file.path(final_dir, "bootstrap_hrs_replicates.csv"))
expected_seed <- 51000L + 2L * 100000L + qc$replicate
stopifnot(
  all(qc$cohort == "HRS"), all(qc$phase == "final"),
  all(qc$seed == expected_seed),
  identical(sort(unique(qc$replicate)), 1:500),
  all(qc[, .N, by = .(replicate, exposure)]$N == 1L),
  all(qc[, uniqueN(exposure), by = replicate]$V1 == 2L),
  all(replicates[, .N, by = .(replicate, exposure)]$N == 1L),
  all(replicates$seed == 51000L + 2L * 100000L + replicates$replicate)
)

audit <- data.table(
  cohort = "HRS", target_replicates = 500L,
  qc_rows = nrow(qc), result_rows = nrow(replicates),
  valid_education = qc[exposure == "education", sum(status == "valid")],
  valid_wealth = qc[exposure == "wealth", sum(status == "valid")],
  failed_education = qc[exposure == "education", sum(status != "valid")],
  failed_wealth = qc[exposure == "wealth", sum(status != "valid")],
  seed_and_uniqueness_checks = TRUE
)
fwrite(audit, file.path(final_dir, "bootstrap_hrs_shard_merge_audit.csv"))
print(audit)
cat("HRS final shards merged and validated.\n")
