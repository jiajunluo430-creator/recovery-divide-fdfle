# Direction 4 exploratory upgrade addendum v2.2

Status: candidate module frozen before its support counts or effect estimates were inspected.

Date: 2026-08-10 (America/Chicago)

## Module E: recovery exhaustion after relapse

Scientific question: is the low-wealth recovery deficit larger after recurrent disability than after the first observed disability episode, consistent with cumulative material constraints on repeated rehabilitation, care, and environmental adaptation?

Rationale: the v1 state history already distinguishes first disability (`D1`) from disability after recovery (`D2`). Estimating a single common recovery gradient may conceal whether socioeconomic inequality becomes more consequential after relapse. This module uses that frozen distinction and does not add or recode any functional item.

Exact risk set and contrast:

- Recovery origins are `D1` and `D2` in the frozen function-transition risk set.
- The exposure is frozen entry-wave wealth (`high`, `middle`, `low`); reported contrasts are low versus high.
- Fit one interval-censored discrete-time complementary log-log recovery model per cohort with `wealth * origin_state`, the v1 age spline, sex, wave/period, proxy, respondent weight, and log interval offset.
- Report the low-versus-high recovery hazard ratio separately for D1 and D2, plus the D2-versus-D1 modification ratio. A ratio below 1 indicates a larger low-wealth recovery disadvantage after relapse.

Support gate: each of the four decisive wealth-by-origin cells (`high/D1`, `low/D1`, `high/D2`, `low/D2`) must contain at least 100 observed recovery events in at least three cohorts. Failure in two or more cohorts is a binding module stop. Middle wealth is retained in model fitting but cannot rescue a decisive-cell failure.

Point promotion gate: either (a) the low-versus-high D2 recovery hazard ratio is at most 0.85 in the same adverse direction in at least two cohorts, or (b) the D2-versus-D1 modification ratio is at most 0.80 in at least two cohorts. Model p-values alone do not promote the module.

Uncertainty gate: primary promotion requires person-bootstrap percentile intervals below 1 for the same qualifying contrast in at least three cohorts. Exactly two interval-supported cohorts may be retained only as a secondary system-specific pattern. Opposing directions block a universal cumulative-exhaustion claim.

Interpretive boundary: even a promoted result is a descriptive transition-history pattern, not proof that wealth causally restores function or that an observed D2 episode is biologically more severe than D1.

## Support result and binding decision

The support audit was run after the rules above were frozen. ELSA and HRS passed, with minimum decisive-cell recovery counts of 248 and 517, respectively. CHARLS failed because the high-wealth D2 cell had 76 recoveries; MHAS failed because the high-wealth D2 cell had 50. With only two of four cohorts passing the required three-cohort gate, Module E is `STOP_TWO_OR_MORE_COHORT_FAILURES`. No effect model or outcome-directed threshold relaxation will be run.
