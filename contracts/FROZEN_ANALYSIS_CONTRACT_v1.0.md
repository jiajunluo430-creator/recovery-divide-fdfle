# Frozen analysis contract v1.0

Frozen: 2026-08-09 (America/Chicago), before Direction4 transition counts or SES model estimates were inspected.

Status: **binding for Gate-0 and formal analysis**. Measurement corrections require a versioned amendment with the error, evidence, and affected outputs. No amendment may be triggered by a favourable or unfavourable estimate.

## 1. Scientific question and claim boundary

The primary question is not whether low socioeconomic status is associated with disability. It is whether education and wealth gaps in remaining disability-free life expectancy (DFLE) at age 60 arise from different rates of:

1. first observed disability/difficulty onset;
2. recovery to independence;
3. relapse after observed recovery;
4. mortality after disability.

The four-country comparison evaluates whether the balance of these processes differs across China (CHARLS), the United States (HRS), England (ELSA), and Mexico (MHAS). Estimates are population transition contrasts under interval observation, not causal effects of changing education, wealth, or a national system.

## 2. Cohorts, waves, population, and time scale

- CHARLS: authorised 2011, 2013, 2015, 2018, and 2020 sources when item-level semantic and linkage gates pass. Harmonized waves are used where available; raw 2020 may supplement only through an audited adapter.
- HRS: RAND HRS 1992–2022 v1, waves 2–16 (1994–2022).
- ELSA: Harmonized ELSA H waves 1–10 (2002–2021); raw wave 11 may supplement only if the same nine-item function state, interview timing, and mortality/observation status can be audited.
- MHAS: Harmonized MHAS D waves 1–6 (2001, 2003, 2012, 2015, 2018, 2021). The 2003–2012 nine-year interval is excluded from primary transitions and retained only as a labelled long-interval sensitivity.
- Eligibility begins at attained age 60. The primary life-table origin is disability-free at exact age 60. Observed intervals start at ages 60–95; extrapolation beyond observed support is prohibited. Life-table integration stops at age 100, with estimates beyond the last supported age flagged.
- Original cohort members remain eligible after entry regardless of later residence when the cohort records a valid interview. Nursing-home residence is not selectively excluded in HRS.
- A person contributes every eligible interval and may contribute more than one episode. Standard errors and bootstrap resampling are clustered by person.

## 3. Socioeconomic exposures

### Education

Education is fixed and mapped to the pre-existing Gateway/RAND three-level construct:

- low: CHARLS/ELSA/MHAS `raeducl` 0–1; HRS `raeduc` 1;
- middle: CHARLS/ELSA/MHAS `raeducl` 2; HRS `raeduc` 2–4;
- high: CHARLS/ELSA/MHAS `raeducl` 3; HRS `raeduc` 5.

The primary inequality contrast is low versus high. Middle is retained for gradient plots and tests. Missing education is excluded from education-specific models and reported.

### Wealth

- Wealth is total household wealth (`h#atotb` for Gateway cohorts; `h#atotw` for RAND HRS, or an audited equivalent).
- Primary wealth is the first valid value observed at or after age 60 and is then fixed, avoiding disability-driven reclassification.
- Within each cohort and entry wave, valid wealth is divided into survey-weighted tertiles. If the relevant positive weight is structurally unavailable, unweighted cut points are used and flagged.
- Primary contrast: lowest versus highest tertile. The middle tertile is retained.
- Missing wealth is excluded from wealth-specific models and reported as a separate flow category, not imputed into a tertile.
- Sensitivity only: time-updated within-wave wealth rank and continuous within-cohort wealth rank.

Education and wealth are separate primary exposure families. They are not mutually adjusted in the headline estimands; mutual adjustment is an interpretation sensitivity, not a rescue model.

## 4. Functional-state definition

The primary state uses the minimum common item set and self/proxy reports accepted by each parent cohort:

- ADL (5): dressing, bathing, eating, transferring in/out of bed, toileting;
- IADL (4): shopping, preparing meals, taking medications, managing money.

Telephone use is excluded from the primary IADL set because it is not in the four-cohort minimum-common construct. Walking across a room is not added as a sixth ADL.

- `independent`: no difficulty in all nine primary items;
- `difficulty`: difficulty in at least one primary item;
- `missing state`: any required item is unresolved and no validated harmonized complete-count variable can establish the state.

Direction2 aggregate variables may be used only after item semantics confirm equivalence. The MHAS 4-IADL construct is not broadened; the other cohorts are narrowed to the same four items. Primary threshold is any difficulty. Prespecified sensitivities are ADL-only, at least two difficulties, and a permissive partial-item rule (any observed difficulty = difficulty; otherwise missing unless all observed items are difficulty-free).

## 5. History-augmented states and transitions

Primary full framework:

- `I0`: independent with no prior observed difficulty;
- `D1`: first or continuing observed difficulty before any recovery;
- `R1`: independent after a prior observed `D1` or relapse episode;
- `D2`: difficulty after an observed `R1` (relapse);
- `DEAD`: absorbing death.

Allowed transitions are stays plus `I0→D1`, `D1→R1`, `R1→D2`, `D2→R1`, and every living state to `DEAD`. Direct `I0→R1`, `I0→D2`, `D1→D2`, and post-death transitions are impossible by construction and must equal zero. For state-occupation summaries, `I0` and `R1` are disability-free; `D1` and `D2` are disabled.

An observed relapse requires the sequence difficulty → independent → difficulty; it is never inferred inside a two-wave interval. A recovery requires difficulty → independent. State history begins at the first valid state, so earlier unobserved episodes are not claimed absent; analyses distinguish observed-history states from lifetime history.

## 6. Interval, interview, death, and attrition rules

- Preferred interval length uses interview year/month (and day if available). Nominal survey year is used only when exact timing is absent, with source flagged.
- Primary living-to-living intervals are 1.0–4.0 years. Intervals below 1 year, above 4 years, nonpositive intervals, duplicate person-waves, and out-of-order observations are excluded from primary models and audited.
- Death is assigned only from a verified cohort death status or death year/date. `iwstat==5` is accepted when its cohort codebook/labels confirm “died this wave.” Missing next interview, disappearance from the panel, or unknown vital status is never death.
- If a verified death date lies after the last living interview and before/at the next scheduled wave, destination is `DEAD`; interval ends at death date. With death year only, midpoint-of-year is used, bounded after the origin interview, and a year-only sensitivity shifts the event to January 1 and December 31.
- Known-alive nonresponse, interviewed-but-function-missing, and unknown vital/observation status are separate categories. Non-death loss is censored, not assigned a functional destination.
- A cohort-wave is death-supported only when death and non-death nonresponse can be distinguished for that interval. Unsupported intervals may inform living-to-living transition Gate-0 counts but not formal mortality/DFLE decomposition.

## 7. Proxy and survey weights

- Proxy interviews are included in the primary analysis because excluding them selectively removes disability and terminal decline. `proxy` is retained as an adjustment/audit field where available.
- Self-respondent-only is a sensitivity. CHARLS harmonized proxy absence is reported as structural unavailability, never coded “self.”
- Positive respondent weight at interval origin is used for weighted descriptive rates, wealth cut points, and primary pseudo-likelihood fits. Invalid/nonpositive weights are excluded only from weighted estimates and retained for unweighted sensitivity/flow.
- Longitudinal weights are used only in an interval-specific sensitivity when their target population and wave span are documented.
- Survey strata/cluster variables are retained when available. Person-cluster bootstrap is the common four-cohort uncertainty method; design-based cohort summaries are secondary when compatible design variables exist.

## 8. Primary estimands and model

For education and wealth separately, within each cohort:

1. origin-specific transition rates per person-year for onset, recovery, relapse, and death;
2. remaining total life expectancy, DFLE, and years with disability at age 60, conditional on `I0` at age 60;
3. low-minus-high (or poorest-minus-richest) absolute gaps in DFLE and years with disability;
4. the contribution in years and percentage of the gap from onset, recovery, relapse, and post-disability mortality.

Models are origin-specific piecewise-exponential competing-transition models. For each allowed destination, a complementary-log-log binomial model uses `log(interval_years)` as an offset and includes SES group, sex, attained age (restricted cubic/natural spline with fixed knots at 60, 70, 80, 90), calendar period/wave, and proxy availability/status. Competing fitted hazards are converted to a coherent one-year transition matrix. Unsupported or nonconvergent sparse terms are pooled only according to the Gate-defined state PIVOT; SES categories and functional thresholds are not collapsed after seeing effects.

Age-specific one-year transition matrices are multiplied from age 60 to 100. Uncertainty uses 500 person-cluster bootstrap replicates for final estimates (100 for computational preview). A replicate must rerun wealth cut points, models, transition matrices, life tables, and decomposition.

### Decomposition

Transition-family contributions use a Shapley replacement decomposition between low and high SES matrices. Blocks are onset, recovery, relapse, post-disability mortality, and pre-disability mortality (accounting closure). The headline four requested mechanisms are reported separately; pre-disability mortality is retained rather than hidden. Contributions sum numerically to the modeled DFLE gap within tolerance 0.01 years. In the three-process PIVOT, relapse is not reported as an estimand.

Cross-country results are coordinated cohort-specific estimates. A pooled mean is secondary; country heterogeneity is not erased by a single fixed-effect estimate.

## 9. Missing data and sensitivities

- Functional states and death are not multiply imputed.
- Primary SES analyses require observed exposure; missingness is reported.
- Minimal structural adjustment uses age, sex, period, and proxy. Broader disease/behavior adjustment is not part of the descriptive inequality estimand.
- Prespecified sensitivities: unweighted; self-response only; ADL-only; at least two difficulties; permissive partial-item state; 1–3-year intervals; exclusion of pandemic-crossing intervals; exact versus nominal time; complete death-supported intervals only; time-updated wealth rank; mutual education/wealth adjustment; non-death observation IPCW; sex stratification only under the upgrade rule.

## 10. Gate-0, PIVOT, and STOP rules

All counts are unweighted persons/events and must be reported by cohort and wave/interval.

### Data integrity gates

- Unique ID-wave after deterministic deduplication; duplicate conflicts are unresolved, not arbitrarily selected.
- At least 95% of living intervals have length 1–4 years after the prespecified MHAS long-gap exclusion, or every excluded interval is source-explainable.
- Impossible transitions and post-death records must be 0 after construction.
- Item/state missingness, proxy availability, weight validity, death support, and non-death attrition must be tabulated.

### Full five-state GO

At least three cohorts must each have:

- at least 200 people with incident `I0→D1`;
- at least 100 people with `D1/D2→R1` recovery;
- at least 50 people with `R1→D2` relapse;
- at least 100 deaths originating from a disabled state; and
- reliable interval-specific distinction of death from non-death loss.

The fourth cohort may remain a living-transition replication cohort if mortality support fails, but it cannot contribute to DFLE/death decomposition.

### Prespecified three-process PIVOT

If relapse fails in two or more cohorts but at least three cohorts each meet incident ≥200, recovery ≥100, disabled-death ≥100, and reliable death/non-death distinction, PIVOT to the frozen three-state `independent ↔ difficulty → death` model. Recovery remains primary; relapse is described only as an unstable Gate-0 count and no relapse contribution is claimed.

### STOP

Stop four-country DFLE mechanism decomposition if fewer than three cohorts reliably distinguish death from non-death loss, or fewer than three cohorts meet the three-process incident/recovery/disabled-death gates. A living-transition report may remain as a transparent feasibility result, but it is not relabelled as the target manuscript.

## 11. Top-journal promotion and positive-result upgrades

Top-journal GO requires more than a significant SES coefficient. At least one must hold with bootstrap uncertainty and prespecified sensitivities:

- recovery plus relapse explains at least 20% and at least 0.5 years of an education or wealth DFLE gap in at least two cohorts, with concordant direction; or
- a mechanism contribution differs across countries by at least 0.75 years in a reproducible, interpretable pattern not driven by one interval or threshold.

If neither holds and the gap is overwhelmingly onset-driven, the story is PIVOT/downgrade rather than model expansion.

Allowed versioned upgrades after a positive gate: sex-specific decomposition; education–wealth discordance; sustained recovery (two consecutive independent waves); policy-relevant country contrasts; population-level counterfactual replacement. Each requires a dated amendment before its estimate is inspected. No data-driven cut point, cohort deletion, or state-definition change is allowed.

## 12. Reporting and release

Required outputs are source/data dictionaries, sample flow, individual-wave construction QC, transition matrices, absolute years and confidence intervals, Shapley closure checks, negative findings, code, logs, and session information. Release packages contain aggregate results and code only; no restricted raw or participant-level derived data.

