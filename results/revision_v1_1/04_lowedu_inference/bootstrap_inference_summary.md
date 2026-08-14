# Household-cluster bootstrap inference summary

Generated: 2026-08-11 20:41:27 CDT
Analysis module: wealth_within_low_education
Functional-difficulty definition: primary
Valid replicates: CHARLS=512; ELSA=500; HRS=500; MHAS=1039

## Key percentile intervals

   cohort                      metric point_estimate     ci_low    ci_high bootstrap_p_two_sided
   <char>                      <char>          <num>      <num>      <num>                 <num>
1: CHARLS        population_fdfle_gap     -2.1771735 -3.0702078 -1.2944612           0.003992016
2: CHARLS population_recovery_relapse     -0.9659235 -1.6936153 -0.2270713           0.015968064
3:   ELSA        population_fdfle_gap     -6.6341708 -7.8394493 -5.6362308           0.003992016
4:   ELSA population_recovery_relapse     -1.8694620 -2.7653398 -1.2210223           0.003992016
5:    HRS        population_fdfle_gap     -4.5634012 -5.6134515 -3.5529134           0.003992016
6:    HRS population_recovery_relapse     -1.1029230 -1.9190536 -0.4567403           0.003992016
7:   MHAS        population_fdfle_gap     -0.9364598 -2.3421585  0.3965118           0.191808192
8:   MHAS population_recovery_relapse      0.3529376 -0.5067114  1.3663588           0.375624376

## Cross-cohort heterogeneity

                        metric cohorts inverse_variance_pooled cochran_q    df p_heterogeneity i2_percent
                        <char>   <int>                   <num>     <num> <int>           <num>      <num>
1:        population_fdfle_gap       4               -3.699360  56.93039     3    2.659326e-12   94.73041
2: population_recovery_relapse       4               -1.004661  13.65661     3    3.411868e-03   78.03261

## Monte Carlo stop rule

Endpoints triggering extension to 1000 replicates: 0
Empty data.table (0 rows and 11 cols): cohort,metric,half_replicates,full_replicates,half_ci_low,half_ci_high...
