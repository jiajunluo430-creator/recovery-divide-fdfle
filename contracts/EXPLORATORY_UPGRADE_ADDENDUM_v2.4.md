# Direction 4 exploratory upgrade addendum v2.4: formal cross-country contrast

Status: frozen after cohort-specific point estimates, the 100-replicate computational preview, and sensitivity directions were known, but before any pairwise bootstrap-difference interval or final 500-replicate interval was calculated.

Date: 2026-08-10 (America/Chicago)

## Scientific question

The C2 wealth-within-low-education point pattern is adverse in CHARLS, ELSA, and HRS and null or opposing in MHAS. This addendum tests whether that apparent difference is quantitatively larger than bootstrap uncertainty rather than describing heterogeneity from signs alone.

## Exact contrast

- Retain the v2.1 module, low-versus-high wealth contrast, functional states, models, age-60 life table, and person bootstrap without modification.
- For each valid replicate shared by a pair of cohorts, subtract the MHAS recovery-plus-relapse contribution from the contribution in CHARLS, ELSA, or HRS.
- The cohort bootstraps are independent; matching deterministic replicate labels is only a reproducible Monte Carlo pairing for the distribution of the difference and does not couple individuals across countries.
- Report analogous pairwise differences for the total DFLE gap and individual Shapley blocks as secondary diagnostics, but promotion uses only the combined recovery-plus-relapse contribution.

## Validity and promotion

- Preview requires at least 90 paired valid replicates; final inference requires at least 450.
- A country-versus-MHAS contrast is supported when the absolute point difference is at least 0.75 years and its 95% percentile interval excludes zero.
- Formal cross-country heterogeneity is promoted when at least two of CHARLS–MHAS, ELSA–MHAS, and HRS–MHAS satisfy that rule in the same direction.
- Model p values and a sign difference alone cannot promote the module.

## Interpretation boundary

A promoted result supports contextual heterogeneity in observed transition-process contributions. It does not identify an effect of a national health system, social programme, insurance regime, family structure, or policy. Country labels are not causal exposures, and the MHAS result is not treated as an outlier to be removed.
