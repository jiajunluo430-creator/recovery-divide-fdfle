# Exploratory bootstrap execution note v2.3: append-only checkpoints

Date: 2026-08-10 (America/Chicago)

Status: implementation-only execution note. It changes no cohort, state, exposure, estimand, model, seed, promotion threshold, or uncertainty rule in v2.0–v2.2.

## Trigger

The project directory is synchronised by Baidu Cloud. During the first final v2 bootstrap launch, the sync client created `*.baiduyun.uploading.cfg` sidecars and held persistent write locks on cohort-level CSV files. Rewriting the same checkpoint file then failed after repeated safe-write attempts. This was an operating-system file-lock condition, not a model or data failure.

Before append-only resumption, complete QC existed for 105 CHARLS replicates and five replicates in each of ELSA, HRS, and MHAS. Some module-metric rows had been written without the corresponding QC checkpoint during an interrupted ELSA write. Such rows are not valid results.

## Frozen execution correction

- Each process still uses the deterministic final seed assigned to cohort and replicate.
- Each replicate still reruns wealth cut points, all three exploratory modules, life tables where applicable, and the specified contrasts.
- A completed replicate is discovered only when QC rows exist for all three modules.
- New results are flushed every five newly evaluated replicates to a filename containing the replicate range, launch timestamp, and process ID. A shard is created once and is never rewritten.
- The final summariser reads the legacy cohort file and every append-only shard, audits overlapping keys, and stops if duplicate QC status, seed, or numeric estimates disagree.
- Module metrics are retained only when a matching selected QC row is explicitly `valid`. Orphan results and metrics from failed module-replicates are removed before intervals are calculated and counted in a machine-readable audit.
- Final validity still requires at least 450 valid replicates of the 500 prespecified draws for every cohort-module.

## Process hygiene

The original CHARLS launcher PID 42380 and its surviving x64 child PID 37648 were stopped after their exact task command line was verified. No source data or completed checkpoints were deleted. Any deterministic overlap between the legacy CHARLS file and append-only shards is retained for the consistency audit described above.
