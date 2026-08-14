# Bootstrap execution note v1.2.1: deterministic HRS sharding

Recorded on 2026-08-10 before final summarisation.

The HRS final bootstrap was stopped after replicate 165 had been checkpointed. Remaining deterministic replicates were executed in two non-overlapping shards: 166–333 and 334–500. The bootstrap seed remained exactly `51000 + 2*100000 + replicate`; resampling unit, wealth-cutpoint recomputation, models, transition matrices, life tables, and Shapley algorithm were unchanged.

Each shard began from identical validated replicates 1–165 in a separate output directory. `06c_merge_hrs_final_shards.R` merges by cohort, phase, replicate, seed, and exposure, removes only identical duplicated checkpoint rows, and requires exactly one QC row per replicate/exposure, both exposures for every replicate 1–500, and the frozen seed identity before writing the formal HRS final files.

This is a wall-time execution optimisation, not an analytic amendment.
