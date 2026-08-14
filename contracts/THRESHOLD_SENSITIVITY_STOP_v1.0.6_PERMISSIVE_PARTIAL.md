# Threshold sensitivity STOP v1.0.6: permissive partial-item state

Frozen on 2026-08-10 before final interpretation.

## Gate result

The permissive partial-item state passed event-count Gate-0 in all four cohorts. It was reconstructed from individual-wave items and produced zero impossible transitions.

## Binding model failure

The unchanged 40-model formal pipeline failed at `MHAS / education / death_post` because the complementary-log-log fit did not meet the frozen converged-and-finite requirement. A complete exact rerun failed at the same model. No SES category, origin state, covariate, knot, or weight was changed.

An isolated diagnostic fit on the same 11,507 intervals and 1,996 deaths converged, showing that event count alone was not the issue; however, the required full reproducible pipeline failed twice at the identical checkpoint. The sensitivity therefore receives `STOP_threshold_model_no_rescue`. Partial estimates from the preceding 35 models are not promoted or combined into life tables.

The ADL-only and at-least-two-difficulties sensitivities completed all 40 models and life-table decompositions. The primary any-difficulty analysis is unaffected.
