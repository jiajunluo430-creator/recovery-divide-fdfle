# Implementation correction v1.0.4: path-consistent pandemic sensitivity

Frozen on 2026-08-10 before interpretation or reporting of sensitivity results.

## Affected component

Only the non-primary pandemic sensitivity in `08_sensitivity_point_models.R` was affected. Gate-0, the five-state histories, primary models, point estimates, and bootstrap analyses were not affected.

## Invalid provisional implementation

The first test removed individual intervals whose origin preceded 2020 and whose endpoint was in 2020 or later, but retained subsequent intervals. Subsequent `R1` and `D2` origin states were constructed from the full observed history and could therefore depend on the omitted bridge observation. This row-level filter did not define a coherent path-level estimand. Its near-zero recovery/relapse contributions in several cohorts are implementation artifacts and must not be interpreted.

## Frozen correction

The sensitivity is renamed `pre_pandemic_end_before_2020`. Both functional-transition and mortality risk sets are restricted to intervals with a finite endpoint strictly before 2020. This preserves internally coherent pre-pandemic histories and changes no state definition, SES mapping, model form, or primary estimand.

The invalid provisional output is retained under `90_invalidated/20260810T014500_pandemic_row_filter_path_inconsistent/`.
