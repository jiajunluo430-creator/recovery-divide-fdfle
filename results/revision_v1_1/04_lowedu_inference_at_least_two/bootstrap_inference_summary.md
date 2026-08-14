# Household-cluster bootstrap inference summary

Generated: 2026-08-11 22:29:39 CDT
Analysis module: wealth_within_low_education
Functional-difficulty definition: at_least_two
Valid replicates: CHARLS=1049; ELSA=1020; HRS=514; MHAS=1023

## Key percentile intervals

   cohort                      metric point_estimate     ci_low     ci_high bootstrap_p_two_sided
   <char>                      <char>          <num>      <num>       <num>                 <num>
1: CHARLS        population_fdfle_gap     -2.2302121 -3.2397526 -1.24196421           0.001998002
2: CHARLS population_recovery_relapse     -0.1811774 -0.7942597  0.45256802           0.659340659
3:   ELSA        population_fdfle_gap     -6.6789783 -8.0216113 -5.48160373           0.001998002
4:   ELSA population_recovery_relapse     -1.1494236 -2.0109788 -0.48793751           0.001998002
5:    HRS        population_fdfle_gap     -4.8815140 -5.9116218 -3.92581738           0.003992016
6:    HRS population_recovery_relapse     -0.6126001 -1.3156361  0.01710645           0.059880240
7:   MHAS        population_fdfle_gap     -2.0752810 -3.5244124 -0.60487303           0.003996004
8:   MHAS population_recovery_relapse     -0.3195398 -0.9850594  0.51590703           0.551448551

## Cross-cohort heterogeneity

                        metric cohorts inverse_variance_pooled cochran_q    df p_heterogeneity i2_percent
                        <char>   <int>                   <num>     <num> <int>           <num>      <num>
1:        population_fdfle_gap       4              -3.9053997 37.506114     3    3.595754e-08   92.00130
2: population_recovery_relapse       4              -0.5386008  4.333059     3    2.276732e-01   30.76485

## Monte Carlo stop rule

Endpoints triggering extension to 1000 replicates: 2
   cohort                      metric half_replicates full_replicates half_ci_low half_ci_high full_ci_low full_ci_high max_endpoint_drift_years relative_halfwidth_change requires_1000
   <char>                      <char>           <int>           <int>       <num>        <num>       <num>        <num>                    <num>                     <num>        <lgcl>
1: CHARLS        population_fdfle_gap             500            1000  -3.1618998   -1.1026164  -3.2397526    -1.241964                0.1393478                0.03078155          TRUE
2: CHARLS population_recovery_relapse             500            1000  -0.6812978    0.4360643  -0.7942597     0.452568                0.1129620                0.10383603          TRUE
