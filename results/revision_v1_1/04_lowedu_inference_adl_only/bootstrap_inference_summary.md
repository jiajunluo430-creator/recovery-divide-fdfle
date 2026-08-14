# Household-cluster bootstrap inference summary

Generated: 2026-08-11 21:37:15 CDT
Analysis module: wealth_within_low_education
Functional-difficulty definition: adl_only
Valid replicates: CHARLS=500; ELSA=1000; HRS=1000; MHAS=1006

## Key percentile intervals

   cohort                      metric point_estimate     ci_low    ci_high bootstrap_p_two_sided
   <char>                      <char>          <num>      <num>      <num>                 <num>
1: CHARLS        population_fdfle_gap     -1.7630841 -2.7216004 -0.8582521           0.003992016
2: CHARLS population_recovery_relapse     -0.2919750 -0.8387163  0.3589201           0.331337325
3:   ELSA        population_fdfle_gap     -6.3475449 -7.6812837 -5.2358201           0.001998002
4:   ELSA population_recovery_relapse     -1.2458285 -2.1634710 -0.5828384           0.003996004
5:    HRS        population_fdfle_gap     -4.5690697 -5.5285065 -3.6124727           0.001998002
6:    HRS population_recovery_relapse     -1.5626430 -2.1882102 -0.8543607           0.001998002
7:   MHAS        population_fdfle_gap     -0.5687411 -1.9943688  0.8360539           0.361638362
8:   MHAS population_recovery_relapse      0.4427696 -0.3833106  1.3763503           0.259740260

## Cross-cohort heterogeneity

                        metric cohorts inverse_variance_pooled cochran_q    df p_heterogeneity i2_percent
                        <char>   <int>                   <num>     <num> <int>           <num>      <num>
1:        population_fdfle_gap       4              -3.2721252  56.24458     3    3.725200e-12   94.66615
2: population_recovery_relapse       4              -0.7333538  16.98015     3    7.134168e-04   82.33231

## Monte Carlo stop rule

Endpoints triggering extension to 1000 replicates: 0
Empty data.table (0 rows and 11 cols): cohort,metric,half_replicates,full_replicates,half_ci_low,half_ci_high...
