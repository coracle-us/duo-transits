# Round 2 test: score (pre-registration 1.1)

Arm A: 123 one-TOI stars scored (146 stars in all). Arm B: 600 trials on 150 hosts.

**Outcome: partial**

| criterion                                                             |       k |   n |   rate |   upper95 |   limit | verdict              |
|:----------------------------------------------------------------------|--------:|----:|-------:|----------:|--------:|:---------------------|
| false exclusion, A2 one-TOI, all                                      |   0.000 | 123 |  0.000 |     0.024 |   0.050 | pass                 |
| false exclusion, A2 one-TOI, SNR 6-8                                  |   0.000 |  17 |  0.000 |     0.162 |   0.050 | reported, not judged |
| false exclusion, A2 one-TOI, SNR 8-12                                 |   0.000 |  15 |  0.000 |     0.181 |   0.050 | reported, not judged |
| false exclusion, A2 one-TOI, SNR 12-20                                |   0.000 |  15 |  0.000 |     0.181 |   0.050 | reported, not judged |
| false exclusion, A2 one-TOI, SNR >20                                  |   0.000 |  41 |  0.000 |     0.070 |   0.050 | reported, not judged |
| false exclusion, B, all                                               |   0.000 | 600 |  0.000 |     0.005 |   0.050 | pass                 |
| false exclusion, B, SNR 6-8                                           |   0.000 |  46 |  0.000 |     0.063 |   0.050 | reported, not judged |
| false exclusion, B, SNR 8-12                                          |   0.000 |  52 |  0.000 |     0.056 |   0.050 | reported, not judged |
| false exclusion, B, SNR 12-20                                         |   0.000 |  71 |  0.000 |     0.041 |   0.050 | pass                 |
| false exclusion, B, SNR >20                                           |   0.000 | 254 |  0.000 |     0.012 |   0.050 | pass                 |
| false exclusion, A2 multi-TOI                                         |   0.000 |  23 |  0.000 |     0.122 |   0.050 | reported, not judged |
| median ln P excluded, A1 one-TOI, SNR >= 8                            | nan     |  71 |  0.654 |   nan     |   0.500 | pass                 |
| median ln P excluded, B, SNR >= 8                                     | nan     | 377 |  0.670 |   nan     |   0.500 | pass                 |
| false event, A1 one-TOI                                               |   7.000 | 123 |  0.057 |     0.104 |   0.050 | FAIL                 |
| false event, B, own events removed                                    |  11.000 | 600 |  0.018 |     0.030 |   0.050 | pass                 |
| false event, B, all events                                            |  84.000 | 600 |  0.140 |     0.165 |   0.050 | reported, not judged |
| passing event in B null runs                                          |  85.000 | 600 |  0.142 |     0.167 |   0.050 | reported, not judged |
| wrong period, A2 one-TOI                                              |   6.000 |  76 |  0.079 |     0.150 |   0.030 | FAIL                 |
| repeat recovery, A2 one-TOI                                           |  81.000 | 123 |  0.659 |     0.729 |   0.600 | pass                 |
| wrong period, B, own events removed                                   |   6.000 | 151 |  0.040 |     0.077 |   0.030 | FAIL                 |
| repeat recovery, B, own events removed                                | 150.000 | 230 |  0.652 |     0.704 |   0.600 | pass                 |
| wrong period, B, all events                                           |  46.000 | 178 |  0.258 |     0.318 |   0.030 | reported, not judged |
| repeat recovery, B, all events                                        | 150.000 | 230 |  0.652 |     0.704 |   0.600 | reported, not judged |
| probability calibration: observed 0, predicted 1.63, 95% interval 0-5 |   0.000 | 723 |  1.628 |   nan     | nan     | pass                 |
