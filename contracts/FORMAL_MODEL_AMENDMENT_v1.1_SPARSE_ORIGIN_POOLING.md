# Formal-model amendment v1.1: mechanism-level origin pooling

Frozen: 2026-08-10 (America/Chicago), after SES-stratified event-support counts and before any SES model coefficient or life-table result was inspected.

Reason: the v1.0 implementation proposed separate `D1 -> R1`, `D2 -> R1`, and living-state-specific death SES slopes. The pre-model support audit found deterministic separation cells: CHARLS high education had only one `D2` interval, and MHAS high education had zero `D2` deaths. Separate SES slopes are therefore not identifiable. This is an estimability correction triggered by risk-set counts, not by coefficient direction, significance, or journal attractiveness.

Binding v1.1 rule, applied uniformly to every cohort and both SES families:

1. Onset remains an `I0` risk-set model for `I0 -> D1`.
2. Recovery uses the combined `D1` and `D2` risk set with event `destination == R1`, an origin-state indicator (`D1` versus `D2`), and one common SES contrast. Origin-specific baseline recovery hazards remain distinct through the state indicator.
3. Relapse remains an `R1` risk-set model for `R1 -> D2`.
4. Pre-disability mortality remains an `I0 -> DEAD` model.
5. Post-disability-history mortality uses the combined `D1`, `R1`, and `D2` mortality risk set, includes an origin-state indicator, and estimates one common SES contrast. Origin-specific baseline mortality hazards remain distinct.

The state space, allowed transitions, age spline, offset, weights, covariates, SES categories, life-table estimands, and five Shapley blocks are unchanged. No cohort or SES category is deleted. Separate-origin models are retained only as a sensitivity where every low/high cell has at least 20 events and converges without separation; they cannot replace the uniform mechanism-level primary model selectively.
