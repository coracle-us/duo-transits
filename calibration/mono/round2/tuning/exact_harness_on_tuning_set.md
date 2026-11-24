# Round 2 test: score (pre-registration 1.1)

Arm A: 76 one-TOI stars scored (109 stars in all). Arm B: 515 trials on 129 hosts.

**Outcome: incomplete**

| criterion                                                             |       k |   n |   rate |   upper95 |   limit | verdict              |
|:----------------------------------------------------------------------|--------:|----:|-------:|----------:|--------:|:---------------------|
| false exclusion, A2 one-TOI, all                                      |   0.000 |  76 |  0.000 |     0.039 |   0.050 | pass                 |
| false exclusion, A2 one-TOI, SNR 6-8                                  |   0.000 |  10 |  0.000 |     0.259 |   0.050 | reported, not judged |
| false exclusion, A2 one-TOI, SNR 8-12                                 |   0.000 |   9 |  0.000 |     0.283 |   0.050 | reported, not judged |
| false exclusion, A2 one-TOI, SNR 12-20                                |   0.000 |  13 |  0.000 |     0.206 |   0.050 | reported, not judged |
| false exclusion, A2 one-TOI, SNR >20                                  |   0.000 |  29 |  0.000 |     0.098 |   0.050 | reported, not judged |
| false exclusion, B, all                                               |   0.000 | 515 |  0.000 |     0.006 |   0.050 | pass                 |
| false exclusion, B, SNR 6-8                                           |   0.000 |  26 |  0.000 |     0.109 |   0.050 | reported, not judged |
| false exclusion, B, SNR 8-12                                          |   0.000 |  47 |  0.000 |     0.062 |   0.050 | reported, not judged |
| false exclusion, B, SNR 12-20                                         |   0.000 |  46 |  0.000 |     0.063 |   0.050 | reported, not judged |
| false exclusion, B, SNR >20                                           |   0.000 | 248 |  0.000 |     0.012 |   0.050 | pass                 |
| false exclusion, A2 multi-TOI                                         |   2.000 |  33 |  0.061 |     0.179 |   0.050 | reported, not judged |
| median ln P excluded, A1 one-TOI, SNR >= 8                            | nan     |  51 |  0.719 |   nan     |   0.500 | pass                 |
| median ln P excluded, B, SNR >= 8                                     | nan     | 341 |  0.666 |   nan     |   0.500 | pass                 |
| false event, A1 one-TOI                                               |   2.000 |  76 |  0.026 |     0.081 |   0.050 | pass                 |
| false event, B, own events removed                                    |   8.000 | 515 |  0.016 |     0.028 |   0.050 | pass                 |
| false event, B, all events                                            |  56.000 | 515 |  0.109 |     0.134 |   0.050 | reported, not judged |
| passing event in B null runs                                          |  55.000 | 515 |  0.107 |     0.132 |   0.050 | reported, not judged |
| wrong period, A2 one-TOI                                              |   0.000 |  52 |  0.000 |     0.056 |   0.030 | pass                 |
| repeat recovery, A2 one-TOI                                           |  61.000 |  76 |  0.803 |     0.874 |   0.600 | pass                 |
| wrong period, B, own events removed                                   |   1.000 | 115 |  0.009 |     0.041 |   0.030 | pass                 |
| repeat recovery, B, own events removed                                | 120.000 | 188 |  0.638 |     0.697 |   0.600 | pass                 |
| wrong period, B, all events                                           |  26.000 | 134 |  0.194 |     0.259 |   0.030 | reported, not judged |
| repeat recovery, B, all events                                        | 120.000 | 188 |  0.638 |     0.697 |   0.600 | reported, not judged |
| probability calibration: observed 0, predicted 1.37, 95% interval 0-4 |   0.000 | 591 |  1.366 |   nan     | nan     | pass                 |
