# Consolidated revision evidence snapshot

Generated: 2026-08-11 22:30:14 CDT

## Focused secondary binding household-bootstrap endpoints

   cohort                      metric point_estimate     ci_low    ci_high valid_replicates
   <char>                      <char>          <num>      <num>      <num>            <int>
1: CHARLS        population_fdfle_gap     -2.1771735 -3.0702078 -1.2944612              500
2: CHARLS population_recovery_relapse     -0.9659235 -1.6936153 -0.2270713              500
3:   ELSA        population_fdfle_gap     -6.6341708 -7.8394493 -5.6362308              500
4:   ELSA population_recovery_relapse     -1.8694620 -2.7653398 -1.2210223              500
5:    HRS        population_fdfle_gap     -4.5634012 -5.6134515 -3.5529134              500
6:    HRS population_recovery_relapse     -1.1029230 -1.9190536 -0.4567403              500
7:   MHAS        population_fdfle_gap     -0.9364598 -2.3421585  0.3965118             1000
8:   MHAS population_recovery_relapse      0.3529376 -0.5067114  1.3663588             1000

## Available primary education/wealth household-bootstrap endpoints

             analysis cohort                      metric point_estimate     ci_low    ci_high valid_replicates
               <char> <char>                      <char>          <num>      <num>      <num>            <int>
 1: primary_education CHARLS        population_fdfle_gap     -6.5584492 -9.8409212 -3.6213315             1000
 2: primary_education CHARLS population_recovery_relapse     -2.0742461 -5.0157319  0.4227745             1000
 3: primary_education   ELSA        population_fdfle_gap     -5.8497942 -6.9837269 -4.8761960              500
 4: primary_education   ELSA population_recovery_relapse     -1.8250683 -2.5112209 -1.2877390              500
 5: primary_education    HRS        population_fdfle_gap     -6.9876237 -7.5506508 -6.3562715              500
 6: primary_education    HRS population_recovery_relapse     -1.7572479 -2.1223374 -1.3736823              500
 7: primary_education   MHAS        population_fdfle_gap     -3.4988795 -6.3246293 -1.2717483             1000
 8: primary_education   MHAS population_recovery_relapse     -0.9053895 -2.4205351  0.3347169             1000
 9:    primary_wealth CHARLS        population_fdfle_gap     -2.6756095 -3.5208228 -1.8765387              500
10:    primary_wealth CHARLS population_recovery_relapse     -1.0188424 -1.6304126 -0.3724336              500
11:    primary_wealth   ELSA        population_fdfle_gap     -7.2089405 -7.9143574 -6.4943334              500
12:    primary_wealth   ELSA population_recovery_relapse     -1.9899120 -2.3950810 -1.6255981              500
13:    primary_wealth    HRS        population_fdfle_gap     -7.2410081 -7.6526471 -6.7992578              500
14:    primary_wealth    HRS population_recovery_relapse     -1.6637101 -1.9379957 -1.3819794              500
15:    primary_wealth   MHAS        population_fdfle_gap     -1.3567373 -2.6755239 -0.1420713             1000
16:    primary_wealth   MHAS population_recovery_relapse      0.1944029 -0.6088929  1.0459741             1000

## Continuous-time panel intensity ratios

    cohort                                         model                                contrast    process intensity_ratio    ci_low   ci_high
    <char>                                        <char>                                  <char>     <char>           <num>     <num>     <num>
 1: CHARLS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education      onset       1.2866064 1.1581150 1.4293538
 2: CHARLS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education  death_pre       1.0421919 0.5085490 2.1358099
 3: CHARLS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education   recovery       0.9938632 0.9005266 1.0968739
 4: CHARLS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education death_post       1.0066418 0.8698059 1.1650045
 5: CHARLS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education    relapse       1.3089243 1.1336987 1.5112330
 6:   ELSA continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education      onset       2.0567134 1.7695413 2.3904896
 7:   ELSA continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education  death_pre       1.3996267 0.7501226 2.6115129
 8:   ELSA continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education   recovery       0.7134740 0.6110992 0.8329992
 9:   ELSA continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education death_post       1.6471867 1.3155683 2.0623969
10:   ELSA continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education    relapse       1.4262449 1.1537347 1.7631217
11:    HRS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education      onset       1.5218749 1.3442819 1.7229298
12:    HRS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education  death_pre       1.6973538 1.2802435 2.2503611
13:    HRS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education   recovery       0.6728371 0.5849911 0.7738746
14:    HRS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education death_post       1.0797271 0.9535528 1.2225969
15:    HRS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education    relapse       1.0943956 0.9183574 1.3041781
16:   MHAS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education      onset       1.1913593 1.0736152 1.3220166
17:   MHAS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education  death_pre       1.2439525 0.9878329 1.5664772
18:   MHAS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education   recovery       0.9982096 0.8872841 1.1230028
19:   MHAS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education death_post       1.1506180 0.9989456 1.3253193
20:   MHAS continuous_time_panel_markov_age_sex_adjusted low_vs_high_wealth_within_low_education    relapse       1.0669411 0.8872865 1.2829717
    cohort                                         model                                contrast    process intensity_ratio    ci_low   ci_high
    <char>                                        <char>                                  <char>     <char>           <num>     <num>     <num>

## Missing evidence-file checks

Empty data.table (0 rows and 3 cols): manuscript_domain,source,exists
