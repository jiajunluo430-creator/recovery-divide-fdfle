# Recovery divide: reproducibility materials

This repository contains the non-disclosive reproducibility materials for:

> **Where socioeconomic gaps in functional-difficulty-free life expectancy arise: onset, recovery, relapse, and mortality across four ageing cohorts**

The coordinated analysis uses the China Health and Retirement Longitudinal Study (CHARLS), English Longitudinal Study of Ageing (ELSA), US Health and Retirement Study (HRS), and Mexican Health and Aging Study (MHAS). It decomposes socioeconomic differences in functional-difficulty-free life expectancy into initial functional composition, onset, recovery, relapse, mortality after difficulty, and mortality before difficulty.

## What is included

- `contracts/`: frozen analysis contracts, source manifests, state definitions, and implementation decisions.
- `analysis/core/`: original Gate-0, multistate, life-table, bootstrap, sensitivity, and reporting code.
- `analysis/revision_v1_1/`: population-initialisation, household/PSU bootstrap, confirmed-state, continuous-time, IPCW, enriched-covariate, calibration, and final-figure code.
- `results/original_analysis/`: aggregate Gate-0, model, bootstrap, life-table, decomposition, and sensitivity outputs.
- `results/revision_v1_1/`: aggregate final and sensitivity outputs used in the manuscript.
- `tables/`: source tables for the main manuscript.
- `figures/`: editable vector figures, figure source data, and structural QC records.
- `environment/`: R session information and the project library manifest.
- `SHA256SUMS.csv`: file-level integrity manifest.

## Data-access boundary

No participant-level data are distributed here. The repository excludes cohort source files, individual-wave panels, transition risk sets, internal person-level survival audits, model checkpoints, and raw bootstrap shards. Access to individual-level CHARLS, ELSA, HRS, and MHAS data is governed by the original study repositories and their data-use terms.

The published CSV files are aggregate results or variable/state metadata. They are sufficient to reproduce the reported tables and figures and to audit numerical consistency. Re-running the full participant-level pipeline requires independent authorised access to the corresponding cohort releases.

## Re-running the analysis

The analysis was run with R 4.5.2. Set the following environment variables before running scripts against authorised local holdings:

- `RECOVERY_DIVIDE_ROOT`: local project/output root.
- `RECOVERY_DIVIDE_CHARLS_FILE`: harmonised CHARLS file.
- `RECOVERY_DIVIDE_ELSA_FILE`: harmonised ELSA file.
- `RECOVERY_DIVIDE_HRS_FILE`: RAND/HRS file used by the analysis.
- `RECOVERY_DIVIDE_MHAS_FILE`: harmonised MHAS file.
- `RECOVERY_DIVIDE_HRS_HARMONIZED_FILE`: optional harmonised HRS source used by schema audits.
- `RECOVERY_DIVIDE_CHARLS_RAW_ROOT` and `RECOVERY_DIVIDE_CHARLS2020_RAW_ROOT`: authorised CHARLS raw-release roots used by the 2020 adapter.
- `RECOVERY_DIVIDE_RLIB`: optional project R library.

The binding order is documented in `contracts/FROZEN_ANALYSIS_CONTRACT_v1.0.md` and `contracts/REVISION_ANALYSIS_CONTRACT_v1.1.md`. The numbered scripts preserve the executed sequence. Scripts in `analysis/core/` establish the source audit, Gate-0, primary models, bootstrap, life tables, and original sensitivities. Scripts 38–56 in `analysis/revision_v1_1/` implement the final estimator and robustness modules reported in the manuscript.

## Software and licences

Analysis code is released under the MIT License. Aggregate tables, figure source data, documentation, and other non-code materials are released under Creative Commons Attribution 4.0 International (CC BY 4.0). Third-party cohort data are not covered by these licences and are not included.

## Citation

Please cite the archived Zenodo record for this release: [https://doi.org/10.5281/zenodo.21938150](https://doi.org/10.5281/zenodo.21938150). The DOI resolves to the exact version deposited for peer review.

## Contact

Correspondence about the analysis may be directed to Xiaolong Liang (`204951@hospital.cqmu.edu.cn`), Fanghui Lu (`lufh@cqmu.edu.cn`), or Ziwei Wang (`drziweiwang@sina.com`).
