# Direction 4 exploratory upgrade addendum v2.1

Status: prospectively binding for bootstrap interpretation; supplements, and does not replace, v1.0 or the v2.0 exploratory charter.

Date: 2026-08-10 (America/Chicago)

## Audit chronology

- The exposure-specific support audit was completed at 10:17 on 2026-08-10 and is preserved in `03_outputs/10_exploratory_upgrade_support`.
- The support audit showed that the original four-cell education-by-wealth discordance module was structurally infeasible before its outcome model was fitted. In CHARLS, the high-education/low-wealth cell contained 9 people, 1 recovery event, and 0 relapse events. The original module C is therefore `STOP_STRUCTURAL_SUPPORT` and will not be rescued by recutting SES or dropping CHARLS.
- Before inspecting a fitted effect estimate, the feasible contrast was narrowed to the wealth gradient within the frozen low-education stratum. The point-model outputs were subsequently written at 10:33 on 2026-08-10. This addendum was recorded after that run to make the earlier support-driven decision explicit; the contrast and thresholds below are frozen before bootstrap uncertainty is inspected.

## Module C2: wealth stratification within low education

Scientific question: among adults already classified in the frozen low-education group, does low rather than high entry-wave wealth create an additional disability-free life-expectancy deficit through impaired recovery and relapse?

Exact contrast: low wealth minus high wealth within `education3_fixed == "low"`. Wealth remains the cohort- and entry-wave-specific v1 tercile classification. Middle wealth contributes to model fitting but is not part of the reported low-versus-high contrast. No education or wealth cut may be changed after seeing results.

Support gate: both decisive cells must contain at least 100 recovery events and 50 relapse events in at least three cohorts. The pre-model audit passed this gate in all four cohorts: the smaller of the two cells contained 646/283 recovery/relapse events in CHARLS, 235/131 in ELSA, 249/159 in HRS, and 570/229 in MHAS.

Point promotion gate: the recovery-plus-relapse Shapley contribution must be at most -0.50 disability-free years in at least two cohorts, indicating a larger low-wealth deficit in the same adverse direction.

Bootstrap promotion gate: in addition to the point gate, the person-bootstrap 95% percentile interval for that contribution must remain below zero in at least two cohorts, with at least 90 valid replicates per cohort in preview and 450 in final analysis. The percentage of the total DFLE gap is reported but is not used to select cohorts.

Interpretive boundary: this module tests residual wealth stratification within low education. It is not evidence that wealth causally compensates for education, and it does not restore the stopped four-cell discordance claim.

## Module A clarification: recovery durability

The phase analysis is wealth-specific because the exposure-specific support audit showed adequate early and sustained recovered-state relapse cells for wealth in all four cohorts. Recovered-state death-by-phase cells were too sparse in some cohorts and remain descriptive/STOP for phase-specific modelling.

Primary durability promotion still requires interval-supported, same-phase wealth differences in at least three cohorts, as specified in v2.0. A replicated result in exactly two cohorts may be reported only as a secondary system-specific finding, not as a four-country durability mechanism.

## Module B clarification: sex heterogeneity

Any promoted sex result must report the signed female-minus-male recovery-plus-relapse contribution in every cohort. Opposing signs across countries block a universal female or male vulnerability claim even if a prespecified magnitude route is met. Such a pattern may support cross-system heterogeneity only when the relevant person-bootstrap interval excludes zero.

## Module D stop record

The frozen linear wealth-by-age interaction specification failed twice for the ELSA onset model: the model reported convergence but at least one coefficient was non-finite on both the main run and the unchanged recheck. The diagnostic is preserved in `03_outputs/11_exploratory_upgrade_point/age_interaction_elsa_onset_recheck.csv`. Module D is `STOP_SECOND_IDENTICAL_FAILURE`; no knot, cohort, or age-range tuning is permitted to rescue it.
