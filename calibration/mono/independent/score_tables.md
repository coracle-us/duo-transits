## Arm A

Stars in the sample file with an output: 120. Finished: 109. Skipped or failed: 11.
- skipped: under 3 transits in clean data: 11 star(s): TIC 113921235, 156514476, 166648874, 257527578, 328934463, 365007485, 388804061, 405010127, 437444661, 4897275, 8260536

### A1 (hidden), all stars: 109 scored trials on 109 stars

- false_exclusion_rate overall: 0/109 = 0.0% (95% UB 2.7%)

| reference SNR | false exclusions | median frac ln P excluded | median sectors |
|---|---|---|---|
| <6 | 0/23 = 0.0% (95% UB 12.2%) | 0.497 | 9 |
| 6-8 | 0/13 = 0.0% (95% UB 20.6%) | 0.598 | 10 |
| 8-12 | 0/17 = 0.0% (95% UB 16.2%) | 0.642 | 5 |
| 12-20 | 0/18 = 0.0% (95% UB 15.3%) | 0.731 | 6 |
| >20 | 0/38 = 0.0% (95% UB 7.6%) | 0.737 | 6 |

- repeat_recovery_rate (trials with >= 1 other true transit at a clean epoch): 0/0 (n/a)
- wrong_period_rate (called resolved or narrowed, truth not in the set): 47/47 = 100.0% (95% UB 100.0%)
  - resolved: 2/2 = 100.0% (95% UB 100.0%)
  - narrowed: 45/45 = 100.0% (95% UB 100.0%)
  - of the 47 wrong: 47 had an alias event that is not a transit of the planet, 0 had a true repeat whose alias was excluded
- calls: no repeat 58, narrowed 45, all excluded 4, resolved 2
- median fraction of ln P (1-1000 d) excluded: 0.703 (quartiles 0.556-0.773)
- median frac ln P excluded by number of sectors: 1-2: 0.539 (n=7), 3-5: 0.618 (n=31), 6-10: 0.732 (n=43), 11-20: 0.762 (n=10), >20: 0.900 (n=18)
- reference measured by the pipeline's box (snr >= 5, depth within 0.33-3x): 91/109; catalogue depth and duration used for the rest

| multi_toi | false exclusions | repeat recovery | wrong period | median frac ln P |
|---|---|---|---|---|
| n_toi_on_star == 1 | 0/76 = 0.0% (95% UB 3.9%) | 0/0 (n/a) | 24/24 = 100.0% (95% UB 100.0%) | 0.700 |
| n_toi_on_star > 1 | 0/33 = 0.0% (95% UB 8.7%) | 0/0 (n/a) | 23/23 = 100.0% (95% UB 100.0%) | 0.707 |

| depth_bin | false exclusions | repeat recovery | wrong period | median frac ln P |
|---|---|---|---|---|
| (2000, 5000] | 0/33 = 0.0% (95% UB 8.7%) | 0/0 (n/a) | 18/18 = 100.0% (95% UB 100.0%) | 0.722 |
| (500, 2000] | 0/39 = 0.0% (95% UB 7.4%) | 0/0 (n/a) | 19/19 = 100.0% (95% UB 100.0%) | 0.606 |
| (5000, 25000] | 0/37 = 0.0% (95% UB 7.8%) | 0/0 (n/a) | 10/10 = 100.0% (95% UB 100.0%) | 0.761 |

### A1, n_toi_on_star == 1: 76 scored trials on 76 stars

- false_exclusion_rate overall: 0/76 = 0.0% (95% UB 3.9%)

| reference SNR | false exclusions | median frac ln P excluded | median sectors |
|---|---|---|---|
| <6 | 0/15 = 0.0% (95% UB 18.1%) | 0.541 | 8 |
| 6-8 | 0/10 = 0.0% (95% UB 25.9%) | 0.621 | 12 |
| 8-12 | 0/8 = 0.0% (95% UB 31.2%) | 0.628 | 6 |
| 12-20 | 0/14 = 0.0% (95% UB 19.3%) | 0.731 | 6 |
| >20 | 0/29 = 0.0% (95% UB 9.8%) | 0.734 | 6 |

- repeat_recovery_rate (trials with >= 1 other true transit at a clean epoch): 0/0 (n/a)
- wrong_period_rate (called resolved or narrowed, truth not in the set): 24/24 = 100.0% (95% UB 100.0%)
  - resolved: 0/0 (n/a)
  - narrowed: 24/24 = 100.0% (95% UB 100.0%)
  - of the 24 wrong: 24 had an alias event that is not a transit of the planet, 0 had a true repeat whose alias was excluded
- calls: no repeat 50, narrowed 24, all excluded 2
- median fraction of ln P (1-1000 d) excluded: 0.700 (quartiles 0.577-0.772)
- median frac ln P excluded by number of sectors: 1-2: 0.544 (n=6), 3-5: 0.625 (n=22), 6-10: 0.750 (n=29), 11-20: 0.729 (n=8), >20: 0.789 (n=11)
- reference measured by the pipeline's box (snr >= 5, depth within 0.33-3x): 64/76; catalogue depth and duration used for the rest

### A1, n_toi_on_star > 1: 33 scored trials on 33 stars

- false_exclusion_rate overall: 0/33 = 0.0% (95% UB 8.7%)

| reference SNR | false exclusions | median frac ln P excluded | median sectors |
|---|---|---|---|
| <6 | 0/8 = 0.0% (95% UB 31.2%) | 0.429 | 9 |
| 6-8 | 0/3 = 0.0% (95% UB 63.2%) | 0.556 | 7 |
| 8-12 | 0/9 = 0.0% (95% UB 28.3%) | 0.642 | 5 |
| 12-20 | 0/4 = 0.0% (95% UB 52.7%) | 0.731 | 6 |
| >20 | 0/9 = 0.0% (95% UB 28.3%) | 0.831 | 15 |

- repeat_recovery_rate (trials with >= 1 other true transit at a clean epoch): 0/0 (n/a)
- wrong_period_rate (called resolved or narrowed, truth not in the set): 23/23 = 100.0% (95% UB 100.0%)
  - resolved: 2/2 = 100.0% (95% UB 100.0%)
  - narrowed: 21/21 = 100.0% (95% UB 100.0%)
  - of the 23 wrong: 23 had an alias event that is not a transit of the planet, 0 had a true repeat whose alias was excluded
- calls: narrowed 21, no repeat 8, resolved 2, all excluded 2
- median fraction of ln P (1-1000 d) excluded: 0.707 (quartiles 0.548-0.782)
- median frac ln P excluded by number of sectors: 1-2: 0.437 (n=1), 3-5: 0.575 (n=9), 6-10: 0.661 (n=14), 11-20: 0.777 (n=2), >20: 0.952 (n=7)
- reference measured by the pipeline's box (snr >= 5, depth within 0.33-3x): 27/33; catalogue depth and duration used for the rest

### A2 (visible), all stars: 109 scored trials on 109 stars

- false_exclusion_rate overall: 7/109 = 6.4% (95% UB 11.7%)

| reference SNR | false exclusions | median frac ln P excluded | median sectors |
|---|---|---|---|
| <6 | 3/23 = 13.0% (95% UB 30.4%) | 0.492 | 9 |
| 6-8 | 0/14 = 0.0% (95% UB 19.3%) | 0.579 | 9 |
| 8-12 | 2/17 = 11.8% (95% UB 32.6%) | 0.676 | 5 |
| 12-20 | 0/17 = 0.0% (95% UB 16.2%) | 0.734 | 6 |
| >20 | 2/38 = 5.3% (95% UB 15.7%) | 0.747 | 6 |

- repeat_recovery_rate (trials with >= 1 other true transit at a clean epoch): 89/109 = 81.7% (95% UB 87.5%)
- wrong_period_rate (called resolved or narrowed, truth not in the set): 16/74 = 21.6% (95% UB 31.0%)
  - resolved: 1/2 = 50.0% (95% UB 97.5%)
  - narrowed: 15/72 = 20.8% (95% UB 30.3%)
  - of the 16 wrong: 12 had an alias event that is not a transit of the planet, 4 had a true repeat whose alias was excluded
- calls: narrowed 72, all excluded 23, no repeat 12, resolved 2
- median fraction of ln P (1-1000 d) excluded: 0.711 (quartiles 0.565-0.785)
- median frac ln P excluded by number of sectors: 1-2: 0.554 (n=7), 3-5: 0.621 (n=31), 6-10: 0.741 (n=43), 11-20: 0.768 (n=10), >20: 0.902 (n=18)
- reference measured by the pipeline's box (snr >= 5, depth within 0.33-3x): 91/109; catalogue depth and duration used for the rest

| multi_toi | false exclusions | repeat recovery | wrong period | median frac ln P |
|---|---|---|---|---|
| n_toi_on_star == 1 | 2/76 = 2.6% (95% UB 8.1%) | 64/76 = 84.2% (95% UB 90.6%) | 7/53 = 13.2% (95% UB 23.4%) | 0.714 |
| n_toi_on_star > 1 | 5/33 = 15.2% (95% UB 29.3%) | 25/33 = 75.8% (95% UB 87.3%) | 9/21 = 42.9% (95% UB 62.8%) | 0.711 |

| depth_bin | false exclusions | repeat recovery | wrong period | median frac ln P |
|---|---|---|---|---|
| (2000, 5000] | 1/33 = 3.0% (95% UB 13.6%) | 30/33 = 90.9% (95% UB 97.5%) | 6/22 = 27.3% (95% UB 46.8%) | 0.728 |
| (500, 2000] | 2/39 = 5.1% (95% UB 15.3%) | 24/39 = 61.5% (95% UB 74.6%) | 7/24 = 29.2% (95% UB 47.9%) | 0.609 |
| (5000, 25000] | 4/37 = 10.8% (95% UB 23.1%) | 35/37 = 94.6% (95% UB 99.0%) | 3/28 = 10.7% (95% UB 25.4%) | 0.763 |

### A2, n_toi_on_star == 1: 76 scored trials on 76 stars

- false_exclusion_rate overall: 2/76 = 2.6% (95% UB 8.1%)

| reference SNR | false exclusions | median frac ln P excluded | median sectors |
|---|---|---|---|
| <6 | 2/15 = 13.3% (95% UB 36.3%) | 0.542 | 8 |
| 6-8 | 0/10 = 0.0% (95% UB 25.9%) | 0.622 | 12 |
| 8-12 | 0/9 = 0.0% (95% UB 28.3%) | 0.676 | 6 |
| 12-20 | 0/13 = 0.0% (95% UB 20.6%) | 0.727 | 7 |
| >20 | 0/29 = 0.0% (95% UB 9.8%) | 0.747 | 6 |

- repeat_recovery_rate (trials with >= 1 other true transit at a clean epoch): 64/76 = 84.2% (95% UB 90.6%)
- wrong_period_rate (called resolved or narrowed, truth not in the set): 7/53 = 13.2% (95% UB 23.4%)
  - resolved: 1/2 = 50.0% (95% UB 97.5%)
  - narrowed: 6/51 = 11.8% (95% UB 21.9%)
  - of the 7 wrong: 5 had an alias event that is not a transit of the planet, 2 had a true repeat whose alias was excluded
- calls: narrowed 51, all excluded 17, no repeat 6, resolved 2
- median fraction of ln P (1-1000 d) excluded: 0.714 (quartiles 0.582-0.780)
- median frac ln P excluded by number of sectors: 1-2: 0.559 (n=6), 3-5: 0.635 (n=22), 6-10: 0.758 (n=29), 11-20: 0.734 (n=8), >20: 0.789 (n=11)
- reference measured by the pipeline's box (snr >= 5, depth within 0.33-3x): 64/76; catalogue depth and duration used for the rest

### A2, n_toi_on_star > 1: 33 scored trials on 33 stars

- false_exclusion_rate overall: 5/33 = 15.2% (95% UB 29.3%)

| reference SNR | false exclusions | median frac ln P excluded | median sectors |
|---|---|---|---|
| <6 | 1/8 = 12.5% (95% UB 47.1%) | 0.433 | 9 |
| 6-8 | 0/4 = 0.0% (95% UB 52.7%) | 0.553 | 6 |
| 8-12 | 2/8 = 25.0% (95% UB 60.0%) | 0.659 | 5 |
| 12-20 | 0/4 = 0.0% (95% UB 52.7%) | 0.734 | 6 |
| >20 | 2/9 = 22.2% (95% UB 55.0%) | 0.843 | 15 |

- repeat_recovery_rate (trials with >= 1 other true transit at a clean epoch): 25/33 = 75.8% (95% UB 87.3%)
- wrong_period_rate (called resolved or narrowed, truth not in the set): 9/21 = 42.9% (95% UB 62.8%)
  - resolved: 0/0 (n/a)
  - narrowed: 9/21 = 42.9% (95% UB 62.8%)
  - of the 9 wrong: 7 had an alias event that is not a transit of the planet, 2 had a true repeat whose alias was excluded
- calls: narrowed 21, no repeat 6, all excluded 6
- median fraction of ln P (1-1000 d) excluded: 0.711 (quartiles 0.549-0.785)
- median frac ln P excluded by number of sectors: 1-2: 0.441 (n=1), 3-5: 0.571 (n=9), 6-10: 0.664 (n=14), 11-20: 0.785 (n=2), >20: 0.955 (n=7)
- reference measured by the pipeline's box (snr >= 5, depth within 0.33-3x): 27/33; catalogue depth and duration used for the rest

### A2 stars where the true period is excluded (for human review): 7

| TIC | TOI | n TOI on star | P true (d) | ref SNR | ref depth (ppm) | ref dur (h) | TOI depth (ppm) | TOI dur (h) | clean true epochs | excluding | excluding epochs (BTJD: box depth) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 118798035 | 7510.01 | 3 | 11.53021 | 21.8 | 4780 | 4.29 | 3881 | 4.81 | 9.0 | 2.0 | 2371.2096:-721+/-301ppm;2382.7398:-48+/-321ppm |
| 149601126 | 2525.01 | 2 | 23.34480 | 8.3 | 11499 | 3.60 | 6474 | 6.41 | 31.0 | 21.0 | 1660.2922:790+/-702ppm;2057.1538:-1525+/-733ppm;2150.5330:852+/-751ppm;2197.2226:-648+/-733ppm;2220.5674:-407+/-733ppm;2243.9122:-386+/-733ppm;2290.6018:-56+/-751ppm;2313.9466:461+/-733ppm;2383.9810:611+/-733ppm;2967.6009:145+/-751ppm;2990.9457:546+/-733ppm;3060.9801:322+/-751ppm;3084.3249:-352+/-751ppm;3131.0145:679+/-733ppm;3201.0489:173+/-742ppm;3667.9449:144+/-751ppm;3714.6345:514+/-733ppm;3737.9793:231+/-733ppm;3878.0481:-61+/-733ppm;3924.7377:323+/-733ppm;3971.4273:507+/-733ppm |
| 170729775 | 2449.01 | 1 | 16.65742 | 1.6 | 10168 | 7.79 | 10168 | 7.79 | 5.0 | 4.0 | 4002.5609:-128+/-65ppm;4035.8758:57+/-65ppm;4285.7371:49+/-65ppm;4302.3945:-35+/-66ppm |
| 349972412 | 4504.01 | 2 | 82.93544 | 0.9 | 13068 | 7.50 | 13068 | 7.50 | 11.0 | 8.0 | 1566.8022:196+/-277ppm;1649.7377:-549+/-277ppm;2064.4149:-160+/-277ppm;2147.3503:-265+/-277ppm;2230.2858:-158+/-277ppm;3723.1237:3+/-277ppm;3888.9946:140+/-277ppm;3971.9301:-279+/-277ppm |
| 459970307 | 1154.01 | 1 | 21.74551 | 2.6 | 659 | 4.20 | 659 | 4.20 | 36.0 | 1.0 | 3339.3325:-487+/-242ppm |
| 55652896 | 216.01 | 2 | 34.50812 | 72.8 | 14492 | 5.40 | 17388 | 3.82 | 24.0 | 13.0 | 1538.4059:1778+/-182ppm;1607.4222:-136+/-182ppm;1676.4384:-170+/-183ppm;2056.0278:-129+/-185ppm;2090.5359:-243+/-183ppm;2125.0440:-73+/-185ppm;2159.5522:2310+/-185ppm;2366.6009:1458+/-182ppm;3712.4177:83+/-182ppm;3850.4502:-215+/-182ppm;3884.9583:27+/-182ppm;3919.4665:154+/-185ppm;4022.9908:-52+/-182ppm |
| 63898957 | 261.02 | 2 | 13.04931 | 8.6 | 543 | 4.92 | 626 | 4.92 | 5.0 | 1.0 | 2133.5487:-94+/-74ppm |

## Arm B

Hosts with an output: 150. Trials scored: 516 on 129 hosts. Trials or hosts not scored: 21.
- error: RuntimeError: no light curve: : 21: TIC 120612828, 125624942, 155216235, 160530638, 161173780, 23547705, 282348445, 284401476, 300293155, 306894786, 30982513, 33494413, 33765894, 350446770, 376796913, 38315113, 390008097, 434315910, 55557788, 77231368, 94487756
- hosts on solar mass and radius: 11; hosts cooler than 3500 K (3500 K coefficients): 3; hosts without a TIC Teff (5772 K used): 1

### B (injections): 516 scored trials on 129 stars

- false_exclusion_rate overall: 4/516 = 0.8% (95% UB 1.8%)

| reference SNR | false exclusions | median frac ln P excluded | median sectors |
|---|---|---|---|
| <6 | 3/148 = 2.0% (95% UB 5.2%) | 0.206 | 4 |
| 6-8 | 0/26 = 0.0% (95% UB 10.9%) | 0.454 | 4 |
| 8-12 | 0/47 = 0.0% (95% UB 6.2%) | 0.637 | 5 |
| 12-20 | 1/47 = 2.1% (95% UB 9.7%) | 0.650 | 4 |
| >20 | 0/248 = 0.0% (95% UB 1.2%) | 0.680 | 5 |

- repeat_recovery_rate (trials with >= 1 other true transit at a clean epoch): 131/188 = 69.7% (95% UB 75.2%)
- wrong_period_rate (called resolved or narrowed, truth not in the set): 83/200 = 41.5% (95% UB 47.5%)
  - resolved: 4/9 = 44.4% (95% UB 74.9%)
  - narrowed: 79/191 = 41.4% (95% UB 47.6%)
  - of the 83 wrong: 82 had an alias event that is not a transit of the planet, 1 had a true repeat whose alias was excluded
- calls: no repeat 286, narrowed 191, all excluded 30, resolved 9
- median fraction of ln P (1-1000 d) excluded: 0.610 (quartiles 0.322-0.712)
- median frac ln P excluded by number of sectors: 1-2: 0.425 (n=92), 3-5: 0.623 (n=236), 6-10: 0.735 (n=140), 11-20: 0.762 (n=16), >20: 0.868 (n=32)
- reference measured by the pipeline's box (snr >= 5, depth within 0.33-3x): 382/516; catalogue depth and duration used for the rest

| inj_depth_bin | false exclusions | repeat recovery | wrong period | median frac ln P |
|---|---|---|---|---|
| 500-2000 | 2/191 = 1.0% (95% UB 3.3%) | 19/64 = 29.7% (95% UB 40.5%) | 33/47 = 70.2% (95% UB 81.0%) | 0.323 |
| 2000-5000 | 0/122 = 0.0% (95% UB 2.4%) | 34/43 = 79.1% (95% UB 88.6%) | 28/60 = 46.7% (95% UB 58.0%) | 0.624 |
| 5000-20000 | 2/203 = 1.0% (95% UB 3.1%) | 78/81 = 96.3% (95% UB 99.0%) | 22/93 = 23.7% (95% UB 32.0%) | 0.680 |

| ref_cadence | false exclusions | repeat recovery | wrong period | median frac ln P |
|---|---|---|---|---|
| QLP 1800 s | 4/25 = 16.0% (95% UB 33.0%) | 7/14 = 50.0% (95% UB 73.6%) | 7/10 = 70.0% (95% UB 91.3%) | 0.535 |
| QLP 200 s | 0/63 = 0.0% (95% UB 4.6%) | 2/5 = 40.0% (95% UB 81.1%) | 9/12 = 75.0% (95% UB 92.8%) | 0.427 |
| QLP 600 s | 0/24 = 0.0% (95% UB 11.7%) | 3/9 = 33.3% (95% UB 65.5%) | 7/8 = 87.5% (95% UB 99.4%) | 0.299 |
| SPOC 120 s | 0/296 = 0.0% (95% UB 1.0%) | 81/112 = 72.3% (95% UB 79.2%) | 39/117 = 33.3% (95% UB 41.2%) | 0.633 |
| TESS-SPOC 1800 s | 0/51 = 0.0% (95% UB 5.7%) | 27/32 = 84.4% (95% UB 93.6%) | 11/33 = 33.3% (95% UB 49.1%) | 0.662 |
| TESS-SPOC 200 s | 0/17 = 0.0% (95% UB 16.2%) | 3/3 = 100.0% (95% UB 100.0%) | 0/4 = 0.0% (95% UB 52.7%) | 0.517 |
| TESS-SPOC 600 s | 0/40 = 0.0% (95% UB 7.2%) | 8/13 = 61.5% (95% UB 83.4%) | 10/16 = 62.5% (95% UB 82.2%) | 0.550 |

| tmag_bin | false exclusions | repeat recovery | wrong period | median frac ln P |
|---|---|---|---|---|
| 10-12 | 0/200 = 0.0% (95% UB 1.5%) | 51/72 = 70.8% (95% UB 79.6%) | 32/74 = 43.2% (95% UB 53.5%) | 0.586 |
| <10 | 0/200 = 0.0% (95% UB 1.5%) | 63/69 = 91.3% (95% UB 96.1%) | 26/89 = 29.2% (95% UB 38.2%) | 0.662 |
| >12 | 4/116 = 3.4% (95% UB 7.7%) | 17/47 = 36.2% (95% UB 49.2%) | 25/37 = 67.6% (95% UB 80.1%) | 0.307 |

| inj_P_bin | false exclusions | repeat recovery | wrong period | median frac ln P |
|---|---|---|---|---|
| 20-50 | 1/115 = 0.9% (95% UB 4.1%) | 52/83 = 62.7% (95% UB 71.5%) | 12/58 = 20.7% (95% UB 31.4%) | 0.576 |
| 50-100 | 2/97 = 2.1% (95% UB 6.3%) | 37/48 = 77.1% (95% UB 86.6%) | 21/52 = 40.4% (95% UB 52.7%) | 0.625 |
| 100-300 | 1/146 = 0.7% (95% UB 3.2%) | 29/39 = 74.4% (95% UB 85.4%) | 21/48 = 43.8% (95% UB 56.6%) | 0.635 |
| 300-1000 | 0/158 = 0.0% (95% UB 1.9%) | 13/18 = 72.2% (95% UB 88.4%) | 29/42 = 69.0% (95% UB 80.6%) | 0.587 |

### B trials where the true period is excluded: 4

| TIC | draw | P (d) | depth (ppm) | b | T14 (h) | ref SNR | ref used | ref depth (ppm) | ref dur (h) | excluding epochs |
|---|---|---|---|---|---|---|---|---|---|---|
| 339672790 | 1.0 | 157.704 | 1476 | 0.47 | 12.08 | 4.4 | catalogue | 1476 | 12.08 | 3757.5565:-790+/-393ppm |
| 382576217 | 0.0 | 79.223 | 7141 | 0.29 | 8.14 | 12.6 | box | 8371 | 9.14 | 3748.0965:-37963+/-940ppm |
| 382576217 | 2.0 | 29.157 | 1707 | 0.66 | 4.52 | 1.4 | catalogue | 1707 | 4.52 | 2043.7490:-10913+/-1075ppm;3764.0111:-51032+/-1113ppm |
| 407663579 | 0.0 | 62.162 | 17560 | 0.31 | 20.40 | 3.8 | catalogue | 17560 | 20.40 | 1603.0156:-44535+/-1203ppm;2348.9630:-160526+/-1096ppm |

## Verdict under pass_criteria and outcomes

- completeness: 109 arm-A stars (need >= 60), 516 arm-B trials (need >= 400): complete
- false exclusion, A2, overall: 7/109 = 6.4% (95% UB 11.7%) -> MISSED
  - SNR 6-8: 0/14 = 0.0% (95% UB 19.3%) -> met
  - SNR 8-12: 2/17 = 11.8% (95% UB 32.6%) -> MISSED
  - SNR 12-20: 0/17 = 0.0% (95% UB 16.2%) -> met
  - SNR >20: 2/38 = 5.3% (95% UB 15.7%) -> MISSED
- false exclusion, B, overall: 4/516 = 0.8% (95% UB 1.8%) -> met
  - SNR 6-8: 0/26 = 0.0% (95% UB 10.9%) -> met
  - SNR 8-12: 0/47 = 0.0% (95% UB 6.2%) -> met
  - SNR 12-20: 1/47 = 2.1% (95% UB 9.7%) -> met
  - SNR >20: 0/248 = 0.0% (95% UB 1.2%) -> met
- wrong period, A1: 47/47 = 100.0% (95% UB 100.0%) -> MISSED
- wrong period, A2: 16/74 = 21.6% (95% UB 31.0%) -> MISSED
- wrong period, B: 83/200 = 41.5% (95% UB 47.5%) -> MISSED

**Outcome: FAIL**: false exclusion A2 overall 7/109; false exclusion A2 SNR 8-12 2/17; false exclusion A2 SNR >20 2/38; wrong period A1 47/47; wrong period A2 16/74; wrong period B 83/200
