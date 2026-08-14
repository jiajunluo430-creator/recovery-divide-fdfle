# Direction 4 exploratory upgrade charter v2.0

Status: versioned exploratory layer. This document does not replace or invalidate the frozen v1.0 primary analysis.

Date: 2026-08-10 (America/Chicago)

## Purpose

Identify a reproducible positive increment that raises the manuscript above a generic socioeconomic disability-free life-expectancy analysis. Exploration may iterate, but every tested module, negative result, and stopping decision remains visible. The v1 primary state, mortality semantics, and main four-cohort estimates remain the audit anchor.

## Immutable safety boundaries

- Read-only source data; project-local writes only.
- Cohorts remain CHARLS, HRS, ELSA, and MHAS. SHARE individual-level data remain out of scope.
- Non-death loss to follow-up is never coded as death.
- The v1 common functional definition remains five ADLs plus four IADLs with any reported difficulty. Exploratory modules may expand state history or interactions but may not silently redefine the underlying items.
- Education and entry-wave wealth remain the v1 harmonised SES measures. Joint-SES analyses may combine their frozen low/high categories but may not recut categories to improve results.
- All exploratory claims are descriptive unless a separate causal identification contract is created.

## Candidate modules and scientific questions

### A. Recovery durability

Question: is the socioeconomic recovery divide driven by failure to attain recovery, failure to sustain the first recovered interval, or relapse after recovery has already persisted?

Operational pilot: classify an R1 origin as `early_recovery` when the immediately preceding scheduled observed state is D1/D2, and `sustained_recovery` when it is R1. Gaps or unavailable preceding states are `unclassified` and excluded from phase contrasts, not reassigned.

Promotion evidence: phase-specific relapse or recovered-state death differences must be supported in at least three cohorts. For the point pilot, low-versus-high wealth relapse must have a hazard ratio of at least 1.15 in the same phase in at least two cohorts, or the sustained-versus-early modification of that hazard ratio must differ by at least 20% in at least two cohorts. Full promotion requires person-bootstrap uncertainty.

### B. Sex-specific recovery divide

Question: do women and men accumulate SES gaps through different onset, recovery, relapse, and death processes?

Operational pilot: retain the v1 risk sets and fit SES-by-sex interactions, standardising life tables separately for women and men. No sex-specific state or SES recoding is allowed.

Promotion evidence: a sex difference in recovery-plus-relapse contribution of at least 0.50 disability-free years in two cohorts, or at least 0.75 years with interval support in one cohort plus directional replication in another. Full promotion requires person-bootstrap uncertainty for the sex contrast.

### C. Education-wealth discordance

Question: can later-life wealth compensate for low education, or can education buffer low wealth, and through which transition process?

Operational pilot: use only the frozen low/high education and wealth categories to form four cells: high-high, low-education/high-wealth, high-education/low-wealth, and low-low. Middle categories remain in descriptive support tables but are not recoded into extremes.

Promotion evidence: all four core cells must have usable transition support in at least three cohorts. A compensation or double-disadvantage contrast must contribute at least 0.50 DFLE years through recovery/relapse in two cohorts with the same direction. Full promotion requires person-bootstrap uncertainty.

### D. Age-varying recovery divide

Question: is the SES recovery/relapse gradient constant after age 60, or does it widen at advanced ages?

Operational pilot: retain the v1 nonlinear main age spline and add a parsimonious SES-by-age interaction using age centred at 70 and capped at 60--90 years. Report age-specific standardised hazards at 60, 70, 80, and 90 years. A nonlinear interaction is considered only if this lower-degree pilot is supported and receives a dated addendum.

Promotion evidence: from age 60 to age 80, the low-versus-high SES recovery or relapse hazard ratio must change by at least 20% in the same direction in at least two cohorts and produce a clinically interpretable age profile. Full promotion requires person-bootstrap uncertainty.

## Support and stopping rules

- A module is `SUPPORTED_FOR_POINT_MODEL` when the decisive low/high or joint-SES cells have at least 100 recovery events and 50 relapse events in at least three cohorts; smaller cells may remain descriptive.
- A module is `PILOT_POSITIVE` only when effect magnitude and directional replication meet its module-specific promotion rule. Model p-values alone cannot promote a module.
- Run full bootstrap only for `PILOT_POSITIVE` modules.
- Stop a module if its defining cells fail support in two or more cohorts, if effects are small and directionally inconsistent, if the full model repeatedly fails without a scientifically defensible unchanged specification, or if a live literature collision removes the claimed increment.
- Do not rescue a module by changing SES cuts, item definitions, cohort inclusion, or contrast direction after results are seen.

## Versioning

Every new module or material change requires a dated addendum stating the scientific reason, inputs, exact contrast, and whether the change was made before or after seeing its result. Negative and stopped modules remain in the audit trail.
