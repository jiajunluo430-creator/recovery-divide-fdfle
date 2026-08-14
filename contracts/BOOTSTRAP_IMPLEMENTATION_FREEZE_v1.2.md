# Bootstrap implementation freeze v1.2

Frozen: 2026-08-10 (America/Chicago), after point-estimate QC and before any bootstrap replicate was run.

## Resampling unit and frame

- Resample people with replacement within cohort from the union of people contributing at least one frozen function or supported mortality risk interval.
- A multinomial person-frequency representation is exactly equivalent to duplicating sampled people. Each person's survey pseudo-weight is multiplied by their bootstrap frequency in every contributed interval.
- One cohort-level draw is shared by education and wealth models within a replicate.
- Wealth entry-wave weighted tertile cut points are recalculated in every replicate using bootstrap-frequency-multiplied entry weights. Education mappings remain fixed.

## Refit requirements

Every replicate reruns all five v1.1 mechanism-level models for education and wealth, standardized annual hazards, transition matrices, age-60 life tables, and the five-block Shapley decomposition. Point-model coefficients or cut points are not reused.

## Replicate validity and failure rule

A cohort-exposure replicate is valid only if all five models converge with finite coefficients, all required low/high hazards are finite and nonnegative, transition rows sum to 1, and Shapley closure is within 0.01 years. Failed replicates are retained in the audit and are not selectively replaced.

- Preview requires at least 90 valid replicates out of 100 for each cohort-exposure family.
- Final inference requires at least 450 valid replicates out of 500 for each cohort-exposure family.
- Falling below either threshold is a bootstrap stability failure for that cohort-exposure result; it cannot be promoted as a precise decomposition.

Seeds are deterministic by cohort and replicate. Preview seeds use base 41000; final seeds use base 51000, with a fixed cohort offset.

## Intervals and promotion screen

Percentile 2.5th and 97.5th bootstrap quantiles are the primary confidence intervals. Bootstrap medians, point-minus-median bias, and valid/failed counts are reported.

The recovery-plus-relapse top-journal screen requires the frozen point threshold (absolute contribution at least 0.5 years and at least 20% of the DFLE gap) plus a 95% bootstrap interval for the combined signed contribution that excludes zero. At least two cohorts must meet this in the same direction. Percentage intervals are reported but are not used alone when the total DFLE gap approaches zero.
