# Household-cluster bootstrap inference summary

Generated: 2026-08-11 21:53:39 CDT
Analysis module: primary_education
Functional-difficulty definition: primary
Valid replicates: CHARLS=1049; ELSA=500; HRS=500; MHAS=1045

## Key percentile intervals

   cohort                      metric point_estimate    ci_low    ci_high bootstrap_p_two_sided
   <char>                      <char>          <num>     <num>      <num>                 <num>
1: CHARLS        population_fdfle_gap     -6.5584492 -9.840921 -3.6213315           0.001998002
2: CHARLS population_recovery_relapse     -2.0742461 -5.015732  0.4227745           0.125874126
3:   ELSA        population_fdfle_gap     -5.8497942 -6.983727 -4.8761960           0.003992016
4:   ELSA population_recovery_relapse     -1.8250683 -2.511221 -1.2877390           0.003992016
5:    HRS        population_fdfle_gap     -6.9876237 -7.550651 -6.3562715           0.003992016
6:    HRS population_recovery_relapse     -1.7572479 -2.122337 -1.3736823           0.003992016
7:   MHAS        population_fdfle_gap     -3.4988795 -6.324629 -1.2717483           0.001998002
8:   MHAS population_recovery_relapse     -0.9053895 -2.420535  0.3347169           0.159840160

## Cross-cohort heterogeneity

                        metric cohorts inverse_variance_pooled cochran_q    df p_heterogeneity i2_percent
                        <char>   <int>                   <num>     <num> <int>           <num>      <num>
1:        population_fdfle_gap       4               -6.608959  9.395935     3       0.0244646    68.0713
2: population_recovery_relapse       4               -1.735754  1.595805     3       0.6603414     0.0000

## Monte Carlo stop rule

Endpoints triggering extension to 1000 replicates: 6
   cohort                                           metric half_replicates full_replicates half_ci_low half_ci_high full_ci_low full_ci_high max_endpoint_drift_years relative_halfwidth_change requires_1000
   <char>                                           <char>           <int>           <int>       <num>        <num>       <num>        <num>                    <num>                     <num>        <lgcl>
1: CHARLS                            conditional_fdfle_gap             500            1000   -9.188408   -3.0380964   -9.195544   -3.1777761                0.1396796                0.02202548          TRUE
2: CHARLS                 population_contribution_recovery             500            1000   -2.798491    0.5605280   -2.852640    0.4480317                0.1124963                0.01767732          TRUE
3: CHARLS                             population_fdfle_gap             500            1000   -9.929170   -3.4845948   -9.840921   -3.6213315                0.1367367                0.03617368          TRUE
4:   MHAS                            conditional_fdfle_gap             500            1000   -6.392034   -0.9070778   -6.206072   -1.0236743                0.1859622                0.05838199          TRUE
5:   MHAS population_contribution_pre_difficulty_mortality             500            1000   -1.524021    1.4474810   -1.377950    1.4172581                0.1460710                0.06307002          TRUE
6:   MHAS                             population_fdfle_gap             500            1000   -6.308387   -1.1526828   -6.324629   -1.2717483                0.1190656                0.02034943          TRUE
