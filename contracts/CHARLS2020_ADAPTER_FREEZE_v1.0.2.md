# CHARLS 2020 raw adapter freeze v1.0.2

Frozen: 2026-08-10 (America/Chicago), after read-only source/semantic linkage audit and before CHARLS wave-4 to wave-5 transition counts or updated model estimates were inspected.

Status: **binding data implementation of the pre-authorized 2020 supplement in section 2 of contract v1.0**. It adds an audited observation boundary; it does not alter the nine-item state, SES contrast, estimand, or model.

## Canonical raw files

All sources remain read-only under the authorized local CHARLS 2020 directory:

- `Sample_Infor.dta`: ID, confirmed current-wave death indicator, interview year/month;
- `Health_Status_and_Functioning.dta`: nine function items and health-module proxy indicator;
- `Weights.dta`: individual base weight;
- `Exit_Module.dta`: available death year/month and calendar type.

`Temp_data/charls20.dta` and `Working_data/charls_wave5.dta` are validation comparators only. They are not canonical analytic inputs.

## Frozen mapping

| Construct | Raw field | Rule |
|---|---|---|
| ID | `ID` | exact character linkage to Harmonized CHARLS D |
| interview time | `iyear`, `imonth` | year plus month midpoint; nominal 2020.5 only for IDs absent from the sample file |
| death status | `died` | 1 = verified death boundary; 0 = known alive/interviewed boundary; absence is unknown, never death |
| dressing | `db001` | 1 = no difficulty; 2–4 = difficulty |
| bathing | `db003` | same |
| eating | `db005` | same |
| bed transfer | `db007` | same |
| toileting | `db009` | same |
| shopping | `db016` | same |
| meal preparation | `db014` | same |
| medications | `db020` | same |
| money management | `db022` | same |
| proxy | `proxy` | labelled 1 = proxy; among health-module records, unlabeled/missing is the survey skip-pattern for self response and is coded self |
| weight | `INDV_weight` | positive individual base weight, matching the Harmonized `r#wtresp` target |
| death time | `exb001_1`, `exb001_2`, `exb002` | solar year/month when valid; year midpoint for lunar/year-only records; otherwise the verified `died` wave boundary |

## Audit evidence before transition inspection

- All canonical files have unique nonmissing IDs.
- The raw health file has 19,367 people; the sample file has 19,395 known-alive interviews and 785 confirmed deaths.
- 19,086 raw health IDs and 762 exit IDs link to the 25,586 Harmonized CHARLS D people.
- All nine raw recodes agree 100% with both independent local derived comparators among 19,347–19,349 jointly observed records.
- Raw and comparator base/adjusted weights agree exactly for all 17,364 linked records.
- Interview year/month agrees exactly for all 19,395 interviewed records.
- All 770 Exit Module IDs are included among confirmed deaths; 107 have valid reported death year/month. The remaining verified deaths use the current-wave boundary and are flagged as coarse timing.

## Analysis inclusion

- Add wave 5 (2020) as the scheduled boundary after Harmonized wave 4 (2018).
- A living functional state requires all nine items. Known alive with missing function can close a mortality interval but cannot define a functional transition.
- Persons absent from both the 2020 sample and verified death files have unknown status and are censored, not classified dead.
- Fixed sex and education are carried from the last valid harmonized record. Wealth remains the first valid value at/after age 60 from harmonized waves; no post-disability wave-5 wealth reclassification is introduced.
- Wave-4 to wave-5 is eligible for primary function and mortality models if the 1–4-year interval and all existing gates pass. Pandemic-crossing exclusion remains a prespecified sensitivity.

## Required rerun

The corrected ELSA outputs remain valid as inputs, but every combined panel, Gate-0 table, risk set, point estimate, and bootstrap output must be regenerated after adding CHARLS wave 5. Any bootstrap started before this freeze is computationally provisional and is archived without interpretation.
