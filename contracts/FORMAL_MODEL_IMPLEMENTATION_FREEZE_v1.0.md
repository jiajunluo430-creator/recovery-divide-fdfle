# Formal-model implementation freeze v1.0

Frozen: 2026-08-09 (America/Chicago), after Gate-0 counts but before any SES coefficient, life-expectancy estimate, or Shapley contribution was inspected.

This file resolves implementation details left open by `FROZEN_ANALYSIS_CONTRACT_v1.0.md`. It does not change the scientific estimands, state definition, exposure cut points, or GO/PIVOT/STOP rules.

## 1. Mortality-supported risk intervals

Mortality models use scheduled adjacent-wave risk intervals and never infer death from disappearance.

- A death event requires the accepted person-level death time after the origin interview and on/before the next scheduled-wave boundary.
- The 1-year minimum applies to observed living-state changes. A verified death may occur within the first year; mortality risk intervals therefore require positive duration and an upper bound of 4 years. This avoids immortal-time selection against early deaths.
- A non-death risk interval is retained only when the next wave contains a valid functional state, `inw==1`, or `iwstat` explicitly indicates respondent/known-alive nonresponse (`1` or `4`). Unknown vital/observation status is excluded from the mortality denominator, not coded alive.
- CHARLS: origin waves 1-3.
- HRS: origin waves 2-15.
- ELSA: origin waves 1-5 (death newly recorded at waves 2-6). Harmonized ELSA has no new `iwstat==5` after wave 6. Wave 11 EOL is a 46%-response subset of known deaths, so its 344 valid death years supplement Gate-0 timing evidence but are not a complete late-period mortality denominator.
- MHAS: origin waves 1-5, excluding the prespecified 2003-2012 interval.

Living-to-living function models may use every primary 1-4-year interval, including ELSA wave 10 to 11, because those transitions condition on an observed destination. Unsupported late ELSA intervals do not enter mortality hazards.

## 2. Fixed SES construction

- Education follows the frozen three-level Gateway/RAND mapping.
- Wealth is the first valid total household wealth observed at/after age 60.
- Wealth tertiles are calculated separately by cohort and wealth-entry wave using the positive origin respondent weight. The weighted empirical distribution is ordered by wealth, and cumulative weight cut points at 1/3 and 2/3 assign low, middle, and high. If no positive weight exists within a cohort-entry-wave stratum, equal-weight ranks are used and the stratum is flagged.
- Ties remain together at the same threshold; no post-result tie breaking or category collapsing is permitted.

## 3. Model specification

For each cohort and SES family, eight origin-specific complementary-log-log models are fitted:

- `I0 -> D1`;
- `D1 -> R1`;
- `D2 -> R1`;
- `R1 -> D2`;
- `I0`, `D1`, `R1`, and `D2` to death.

Each model includes SES, sex, attained-age natural spline (boundary 60 and 95; internal knots 70, 80, 90), origin wave as a categorical period term, and proxy status when it varies. `log(interval_years)` is the offset. Positive origin respondent weights are rescaled to mean 1 within the fitted risk set. Missing/nonpositive-weight rows are excluded from the weighted primary fit and retained for unweighted sensitivity.

Factor levels are retained as frozen categories. A nonconvergent model or a required SES level with zero risk/events is a model failure; it is not repaired by collapsing SES or redefining function.

## 4. Standardization and transition matrices

At each integer attained age, transition hazards are standardized over the empirical weighted distribution of sex, origin wave, and proxy status in the corresponding fitted origin-state risk set. SES is set to the requested low/middle/high level. This origin-specific standardization is descriptive and does not represent a causal intervention.

Competing hazards from each origin are converted to one-year probabilities using the total-hazard formula, so every row of each transition matrix sums to 1 and death is absorbing.

## 5. Life-table integration and tail rule

- Start at exact age 60 in `I0`.
- Multiply one-year transition matrices through age 100.
- Use half-cycle state occupancy for each year.
- Observed-origin age support ends at 95. For ages 95-99, all hazards are held constant at their age-95 value rather than allowing spline extrapolation. These tail years are flagged, and residual survival at age 100 is reported.
- Total life expectancy, DFLE, and disabled years are therefore age-100-truncated estimates. Their sum must agree within numerical tolerance.

## 6. Shapley decomposition

The low-minus-high DFLE gap is decomposed over five fixed blocks: onset; recovery (`D1/D2 -> R1`); relapse; post-disability mortality (`D1/R1/D2 -> death`); and pre-disability mortality (`I0 -> death`). All 120 replacement orders are averaged. Signed contributions must close to the modeled low-minus-high gap within 0.01 years.

## 7. Uncertainty sequence

- Point estimates are run first for computational and convergence QC.
- A 100-replicate person-cluster bootstrap is the locked preview.
- If point models and preview pass, the final run uses 500 person-cluster replicates. Each replicate redraws people within cohort and reruns wealth cut points, models, matrices, life tables, and decomposition.
