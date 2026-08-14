# Implementation correction v1.0.3: CHARLS 2020 proxy skip pattern

Frozen: 2026-08-10 (America/Chicago), during missingness QC and before any bootstrap preview interval was summarized.

Status: **binding implementation correction** to make code conform to the already frozen v1.0.2 raw adapter. No measurement definition or model term changes.

## Detected mismatch

The raw CHARLS 2020 health-module `proxy` field has value label `1 = Yes`. Among 19,367 health-module records, 1,784 carry value 1 and 17,583 are structurally blank under the survey skip pattern. The v1.0.2 adapter explicitly froze value 1 as proxy and blank among a present health interview as self response.

The first implementation used `fifelse(proxy == 1, 1, 0)`. In R/data.table, an `NA` test returns `NA`, so structurally blank self responses were incorrectly retained as unavailable rather than coded 0. Wave-5 missingness QC exposed this as 90.8% proxy missingness.

## Correction

For a present raw health-module record:

- nonmissing `proxy == 1` → `proxy = 1`;
- otherwise → `proxy = 0` (self response).

For a 2020 sample record with no health-module record, proxy remains unavailable. All pre-2020 CHARLS proxy fields remain structurally unavailable as specified in v1.0.

## Consequences

Functional states, death status, Gate-0 transition counts, and weights are unchanged. Proxy-adjusted risk sets and every point/bootstrap model are regenerated because the adjustment design matrix changes. Provisional outputs from the misimplemented version are archived and are ineligible for reporting.
