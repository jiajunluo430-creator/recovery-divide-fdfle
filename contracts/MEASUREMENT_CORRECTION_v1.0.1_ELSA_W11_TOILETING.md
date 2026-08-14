# Measurement correction v1.0.1: ELSA wave 11 toileting item

Frozen: 2026-08-10 (America/Chicago), immediately after semantic discovery and before any bootstrap uncertainty result was summarized or used for promotion.

Status: **binding measurement correction** under section 12 of the frozen v1.0 contract. This correction restores the prespecified construct; it does not change the construct.

## Error

The initial raw ELSA wave 11 adapter mapped internal token `toilta` to `headlwa`. The source label for `headlwa` is walking across a room. The frozen contract explicitly requires toileting as the fifth ADL and explicitly prohibits adding walking across a room.

## Source evidence

Read-only Stata variable-metadata inspection of `wave_11_elsa_data_eul.dta` identified:

- `headlwc`: "ADL: difficulty using the toilet, including getting up or down";
- `headlwa`: walking-across-a-room ADL, which is outside the frozen nine-item common set.

The raw wave 11 ELSA adapter is therefore corrected to map `toilta = headlwc`. `headlwa` remains available only for audit and is excluded from the primary functional state.

## Affected outputs and required invalidation

All outputs downstream of `01_build_gate0_panel.R` created before this correction are provisional and must be regenerated:

1. individual-wave panel and transition intervals;
2. Gate-0 event counts and transition audits;
3. formal function and mortality risk sets;
4. point models, life tables, and Shapley decomposition;
5. bootstrap preview files produced from the provisional risk sets.

The running provisional bootstrap is terminated and its files are retained in a clearly labelled audit archive; they are not eligible for interpretation or reporting.

## Non-result-driven rationale

The correction was triggered by a direct item-label mismatch while constructing the variable dictionary. No bootstrap confidence interval or promotion decision had been inspected. The prespecified state threshold, cohort set, SES definition, event gates, and model are unchanged.
