# Single-transit calibration: score (pre-registration version 1.1)

Arm A: 109 stars scored of 120 drawn (skipped: {'under 3 clean transits': 11}; errors: 0).
Arm B: 515 trials scored on 129 of 150 hosts (errors: 0).

- A1 wrong period: 48.0 of 48 (100.0%, 95% upper bound 100.0%) - FAIL
- A2 false exclusion, SNR all: 7.0 of 109 (6.4%, 95% upper bound 11.7%) - FAIL
- A2 false exclusion, SNR 6-8: 0.0 of 14 (0.0%, 95% upper bound 19.3%) - pass
- A2 false exclusion, SNR 8-12: 2.0 of 17 (11.8%, 95% upper bound 32.6%) - FAIL
- A2 false exclusion, SNR 12-20: 0.0 of 17 (0.0%, 95% upper bound 16.2%) - pass
- A2 false exclusion, SNR >20: 2.0 of 38 (5.3%, 95% upper bound 15.7%) - FAIL
- A2 wrong period: 15.0 of 74 (20.3%, 95% upper bound 29.5%) - FAIL
- B false exclusion, SNR all: 4.0 of 515 (0.8%, 95% upper bound 1.8%) - pass
- B false exclusion, SNR 6-8: 0.0 of 26 (0.0%, 95% upper bound 10.9%) - pass
- B false exclusion, SNR 8-12: 0.0 of 47 (0.0%, 95% upper bound 6.2%) - pass
- B false exclusion, SNR 12-20: 1.0 of 46 (2.2%, 95% upper bound 9.9%) - pass
- B false exclusion, SNR >20: 0.0 of 248 (0.0%, 95% upper bound 1.2%) - pass
- B wrong period: 84.0 of 201 (41.8%, 95% upper bound 47.8%) - FAIL

**Outcome: fail**

## All rates

| arm                 | metric                              | snr_bin   |   n |       k |    rate |   upper95 |
|:--------------------|:------------------------------------|:----------|----:|--------:|--------:|----------:|
| A1                  | false_exclusion                     | all       | 109 |   0.000 |   0.000 |     0.027 |
| A1                  | false_exclusion                     | <6        |  23 |   0.000 |   0.000 |     0.122 |
| A1                  | false_exclusion                     | 6-8       |  13 |   0.000 |   0.000 |     0.206 |
| A1                  | false_exclusion                     | 8-12      |  17 |   0.000 |   0.000 |     0.162 |
| A1                  | false_exclusion                     | 12-20     |  18 |   0.000 |   0.000 |     0.153 |
| A1                  | false_exclusion                     | >20       |  38 |   0.000 |   0.000 |     0.076 |
| A1                  | false_exclusion_given_a_clean_epoch | all       |   0 |   0.000 | nan     |   nan     |
| A1                  | repeat_recovery                     | all       |   0 |   0.000 | nan     |   nan     |
| A1                  | repeat_recovery                     | <6        |   0 |   0.000 | nan     |   nan     |
| A1                  | repeat_recovery                     | 6-8       |   0 |   0.000 | nan     |   nan     |
| A1                  | repeat_recovery                     | 8-12      |   0 |   0.000 | nan     |   nan     |
| A1                  | repeat_recovery                     | 12-20     |   0 |   0.000 | nan     |   nan     |
| A1                  | repeat_recovery                     | >20       |   0 |   0.000 | nan     |   nan     |
| A1                  | wrong_period                        | all       |  48 |  48.000 |   1.000 |     1.000 |
| A1                  | wrong_period_resolved_only          | all       |   2 |   2.000 |   1.000 |     1.000 |
| A1                  | median_frac_lnP_excluded            | all       | 109 | nan     |   0.703 |   nan     |
| A1                  | median_frac_lnP_excluded            | <6        |  23 | nan     |   0.497 |   nan     |
| A1                  | median_frac_lnP_excluded            | 6-8       |  13 | nan     |   0.598 |   nan     |
| A1                  | median_frac_lnP_excluded            | 8-12      |  17 | nan     |   0.642 |   nan     |
| A1                  | median_frac_lnP_excluded            | 12-20     |  18 | nan     |   0.731 |   nan     |
| A1                  | median_frac_lnP_excluded            | >20       |  38 | nan     |   0.737 |   nan     |
| A2                  | false_exclusion                     | all       | 109 |   7.000 |   0.064 |     0.117 |
| A2                  | false_exclusion                     | <6        |  23 |   3.000 |   0.130 |     0.304 |
| A2                  | false_exclusion                     | 6-8       |  14 |   0.000 |   0.000 |     0.193 |
| A2                  | false_exclusion                     | 8-12      |  17 |   2.000 |   0.118 |     0.326 |
| A2                  | false_exclusion                     | 12-20     |  17 |   0.000 |   0.000 |     0.162 |
| A2                  | false_exclusion                     | >20       |  38 |   2.000 |   0.053 |     0.157 |
| A2                  | false_exclusion_given_a_clean_epoch | all       | 109 |   7.000 |   0.064 |     0.117 |
| A2                  | repeat_recovery                     | all       | 109 |  88.000 |   0.807 |     0.867 |
| A2                  | repeat_recovery                     | <6        |  23 |   6.000 |   0.261 |     0.451 |
| A2                  | repeat_recovery                     | 6-8       |  14 |  11.000 |   0.786 |     0.939 |
| A2                  | repeat_recovery                     | 8-12      |  17 |  16.000 |   0.941 |     0.997 |
| A2                  | repeat_recovery                     | 12-20     |  17 |  17.000 |   1.000 |     1.000 |
| A2                  | repeat_recovery                     | >20       |  38 |  38.000 |   1.000 |     1.000 |
| A2                  | wrong_period                        | all       |  74 |  15.000 |   0.203 |     0.295 |
| A2                  | wrong_period_resolved_only          | all       |   0 |   0.000 | nan     |   nan     |
| A2                  | median_frac_lnP_excluded            | all       | 109 | nan     |   0.711 |   nan     |
| A2                  | median_frac_lnP_excluded            | <6        |  23 | nan     |   0.492 |   nan     |
| A2                  | median_frac_lnP_excluded            | 6-8       |  14 | nan     |   0.579 |   nan     |
| A2                  | median_frac_lnP_excluded            | 8-12      |  17 | nan     |   0.676 |   nan     |
| A2                  | median_frac_lnP_excluded            | 12-20     |  17 | nan     |   0.734 |   nan     |
| A2                  | median_frac_lnP_excluded            | >20       |  38 | nan     |   0.747 |   nan     |
| B                   | false_exclusion                     | all       | 515 |   4.000 |   0.008 |     0.018 |
| B                   | false_exclusion                     | <6        | 145 |   3.000 |   0.021 |     0.053 |
| B                   | false_exclusion                     | 6-8       |  26 |   0.000 |   0.000 |     0.109 |
| B                   | false_exclusion                     | 8-12      |  47 |   0.000 |   0.000 |     0.062 |
| B                   | false_exclusion                     | 12-20     |  46 |   1.000 |   0.022 |     0.099 |
| B                   | false_exclusion                     | >20       | 248 |   0.000 |   0.000 |     0.012 |
| B                   | false_exclusion_given_a_clean_epoch | all       | 189 |   4.000 |   0.021 |     0.048 |
| B                   | repeat_recovery                     | all       | 188 | 131.000 |   0.697 |     0.752 |
| B                   | repeat_recovery                     | <6        |  49 |   6.000 |   0.122 |     0.227 |
| B                   | repeat_recovery                     | 6-8       |  10 |   4.000 |   0.400 |     0.696 |
| B                   | repeat_recovery                     | 8-12      |  17 |  12.000 |   0.706 |     0.876 |
| B                   | repeat_recovery                     | 12-20     |  11 |  10.000 |   0.909 |     0.995 |
| B                   | repeat_recovery                     | >20       | 100 |  98.000 |   0.980 |     0.996 |
| B                   | wrong_period                        | all       | 201 |  84.000 |   0.418 |     0.478 |
| B                   | wrong_period_resolved_only          | all       |   8 |   3.000 |   0.375 |     0.711 |
| B                   | median_frac_lnP_excluded            | all       | 515 | nan     |   0.612 |   nan     |
| B                   | median_frac_lnP_excluded            | <6        | 145 | nan     |   0.203 |   nan     |
| B                   | median_frac_lnP_excluded            | 6-8       |  26 | nan     |   0.454 |   nan     |
| B                   | median_frac_lnP_excluded            | 8-12      |  47 | nan     |   0.637 |   nan     |
| B                   | median_frac_lnP_excluded            | 12-20     |  46 | nan     |   0.654 |   nan     |
| B                   | median_frac_lnP_excluded            | >20       | 248 | nan     |   0.680 |   nan     |
| A1 single-TOI stars | false_exclusion                     | all       |  76 |   0.000 |   0.000 |     0.039 |
| A1 single-TOI stars | false_exclusion                     | <6        |  15 |   0.000 |   0.000 |     0.181 |
| A1 single-TOI stars | false_exclusion                     | 6-8       |  10 |   0.000 |   0.000 |     0.259 |
| A1 single-TOI stars | false_exclusion                     | 8-12      |   8 |   0.000 |   0.000 |     0.312 |
| A1 single-TOI stars | false_exclusion                     | 12-20     |  14 |   0.000 |   0.000 |     0.193 |
| A1 single-TOI stars | false_exclusion                     | >20       |  29 |   0.000 |   0.000 |     0.098 |
| A1 single-TOI stars | false_exclusion_given_a_clean_epoch | all       |   0 |   0.000 | nan     |   nan     |
| A1 single-TOI stars | repeat_recovery                     | all       |   0 |   0.000 | nan     |   nan     |
| A1 single-TOI stars | repeat_recovery                     | <6        |   0 |   0.000 | nan     |   nan     |
| A1 single-TOI stars | repeat_recovery                     | 6-8       |   0 |   0.000 | nan     |   nan     |
| A1 single-TOI stars | repeat_recovery                     | 8-12      |   0 |   0.000 | nan     |   nan     |
| A1 single-TOI stars | repeat_recovery                     | 12-20     |   0 |   0.000 | nan     |   nan     |
| A1 single-TOI stars | repeat_recovery                     | >20       |   0 |   0.000 | nan     |   nan     |
| A1 single-TOI stars | wrong_period                        | all       |  25 |  25.000 |   1.000 |     1.000 |
| A1 single-TOI stars | wrong_period_resolved_only          | all       |   0 |   0.000 | nan     |   nan     |
| A1 single-TOI stars | median_frac_lnP_excluded            | all       |  76 | nan     |   0.700 |   nan     |
| A1 single-TOI stars | median_frac_lnP_excluded            | <6        |  15 | nan     |   0.541 |   nan     |
| A1 single-TOI stars | median_frac_lnP_excluded            | 6-8       |  10 | nan     |   0.621 |   nan     |
| A1 single-TOI stars | median_frac_lnP_excluded            | 8-12      |   8 | nan     |   0.628 |   nan     |
| A1 single-TOI stars | median_frac_lnP_excluded            | 12-20     |  14 | nan     |   0.731 |   nan     |
| A1 single-TOI stars | median_frac_lnP_excluded            | >20       |  29 | nan     |   0.734 |   nan     |
| A2 single-TOI stars | false_exclusion                     | all       |  76 |   2.000 |   0.026 |     0.081 |
| A2 single-TOI stars | false_exclusion                     | <6        |  15 |   2.000 |   0.133 |     0.363 |
| A2 single-TOI stars | false_exclusion                     | 6-8       |  10 |   0.000 |   0.000 |     0.259 |
| A2 single-TOI stars | false_exclusion                     | 8-12      |   9 |   0.000 |   0.000 |     0.283 |
| A2 single-TOI stars | false_exclusion                     | 12-20     |  13 |   0.000 |   0.000 |     0.206 |
| A2 single-TOI stars | false_exclusion                     | >20       |  29 |   0.000 |   0.000 |     0.098 |
| A2 single-TOI stars | false_exclusion_given_a_clean_epoch | all       |  76 |   2.000 |   0.026 |     0.081 |
| A2 single-TOI stars | repeat_recovery                     | all       |  76 |  63.000 |   0.829 |     0.896 |
| A2 single-TOI stars | repeat_recovery                     | <6        |  15 |   4.000 |   0.267 |     0.511 |
| A2 single-TOI stars | repeat_recovery                     | 6-8       |  10 |   8.000 |   0.800 |     0.963 |
| A2 single-TOI stars | repeat_recovery                     | 8-12      |   9 |   9.000 |   1.000 |     1.000 |
| A2 single-TOI stars | repeat_recovery                     | 12-20     |  13 |  13.000 |   1.000 |     1.000 |
| A2 single-TOI stars | repeat_recovery                     | >20       |  29 |  29.000 |   1.000 |     1.000 |
| A2 single-TOI stars | wrong_period                        | all       |  53 |   6.000 |   0.113 |     0.211 |
| A2 single-TOI stars | wrong_period_resolved_only          | all       |   0 |   0.000 | nan     |   nan     |
| A2 single-TOI stars | median_frac_lnP_excluded            | all       |  76 | nan     |   0.714 |   nan     |
| A2 single-TOI stars | median_frac_lnP_excluded            | <6        |  15 | nan     |   0.542 |   nan     |
| A2 single-TOI stars | median_frac_lnP_excluded            | 6-8       |  10 | nan     |   0.623 |   nan     |
| A2 single-TOI stars | median_frac_lnP_excluded            | 8-12      |   9 | nan     |   0.676 |   nan     |
| A2 single-TOI stars | median_frac_lnP_excluded            | 12-20     |  13 | nan     |   0.727 |   nan     |
| A2 single-TOI stars | median_frac_lnP_excluded            | >20       |  29 | nan     |   0.747 |   nan     |
