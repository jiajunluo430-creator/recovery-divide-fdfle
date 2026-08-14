# Public release quality control

Release version: **1.0.0**  
Validation date: **2026-08-14**

## Privacy and disclosure boundary

- Explicit whitelist assembly: PASS
- Participant-level cohort source files: 0
- Participant-level panels or transition risk sets: 0
- Internal person-level survival audit files: 0
- Raw bootstrap shards or model checkpoints: 0
- Restricted data extensions (`.rds`, `.dta`, `.sav`, `.sas7bdat`, `.xpt`, `.parquet`, and related formats): 0
- CSV headers containing person, household, respondent, or case identifiers: 0
- Host-specific `F:`, `D:`, or user-profile absolute paths: 0
- Credential-pattern scan: PASS

## Reproducibility checks

- Public R scripts parsed with R 4.5.2: 45/45
- R parse errors: 0
- Main and supplementary vector-figure QC: PASS
- SHA-256 file manifest: present

The public archive contains only code, contracts, variable/state metadata, aggregate results, tables, vector figures, figure source data, and session information. Independent end-to-end execution requires authorised access to the original cohort data.
