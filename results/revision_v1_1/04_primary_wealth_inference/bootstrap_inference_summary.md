# Household-cluster bootstrap inference summary

Generated: 2026-08-11 22:14:37 CDT
Analysis module: primary_wealth
Functional-difficulty definition: primary
Valid replicates: CHARLS=513; ELSA=500; HRS=500; MHAS=1045

## Key percentile intervals

   cohort                      metric point_estimate     ci_low    ci_high bootstrap_p_two_sided
   <char>                      <char>          <num>      <num>      <num>                 <num>
1: CHARLS        population_fdfle_gap     -2.6756095 -3.5208228 -1.8765387           0.003992016
2: CHARLS population_recovery_relapse     -1.0188424 -1.6304126 -0.3724336           0.003992016
3:   ELSA        population_fdfle_gap     -7.2089405 -7.9143574 -6.4943334           0.003992016
4:   ELSA population_recovery_relapse     -1.9899120 -2.3950810 -1.6255981           0.003992016
5:    HRS        population_fdfle_gap     -7.2410081 -7.6526471 -6.7992578           0.003992016
6:    HRS population_recovery_relapse     -1.6637101 -1.9379957 -1.3819794           0.003992016
7:   MHAS        population_fdfle_gap     -1.3567373 -2.6755239 -0.1420713           0.035964036
8:   MHAS population_recovery_relapse      0.1944029 -0.6088929  1.0459741           0.557442557

## Cross-cohort heterogeneity

                        metric cohorts inverse_variance_pooled cochran_q    df p_heterogeneity i2_percent
                        <char>   <int>                   <num>     <num> <int>           <num>      <num>
1:        population_fdfle_gap       4               -6.161792 159.62953     3    2.203334e-34   98.12065
2: population_recovery_relapse       4               -1.571302  26.68517     3    6.853451e-06   88.75780

## Monte Carlo stop rule

Endpoints triggering extension to 1000 replicates: 0
Empty data.table (0 rows and 11 cols): cohort,metric,half_replicates,full_replicates,half_ci_low,half_ci_high...
