# Implementation clarification v1.0.5: retain both pandemic sensitivities and audit support

Frozen on 2026-08-10 before final bootstrap reporting. This clarification supersedes the decision in v1.0.4 to treat the interval-only exclusion as invalid.

## Reason for clarification

The contract prespecified exclusion of pandemic-crossing intervals. A later origin state can remain a valid observed-history state even when a prior interval is excluded from a sensitivity risk set; therefore the exact interval-only sensitivity must be retained rather than discarded. A stricter pre-pandemic endpoint restriction is useful as a supplementary path-restricted check but does not replace the prespecified analysis.

## Analyses retained

1. `exclude_pandemic_crossing`: remove only intervals with origin before 2020 and endpoint in 2020 or later, as originally frozen.
2. `pre_pandemic_end_before_2020`: supplementary restriction to intervals ending before 2020.

Both use the unchanged five-state histories, SES mappings, model form, standardization, and life-table code. Neither is eligible to rescue or overturn the primary result by model retuning.

## Required numerical-support audit

For each sensitivity, cohort, exposure, and process, persist the maximum and 99th-percentile standardized annual hazard and counts above 10 and 100. The first runs showed that removing key period bridges can create unsupported age-by-period extrapolation in disabled-state mortality despite directionally concordant SES coefficients. Such a run is retained and reported as numerically non-informative rather than interpreted as a null mechanism result.
