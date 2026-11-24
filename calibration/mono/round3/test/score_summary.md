# Round 3 test: score (pre-registration 1.0)

Arm A: 90 one-TOI stars scored (116 stars in all). Arm B: 1072 trials on 269 hosts.

**Outcome: not defined by the pre-registration: wrong period (true event) and alias power met, repeat recovery or false exclusion missed**

| criterion                                                                       |       k |    n |   rate |   upper95 |   limit | verdict              |
|:--------------------------------------------------------------------------------|--------:|-----:|-------:|----------:|--------:|:---------------------|
| wrong period, A2 one-TOI                                                        |   1.000 |   49 |  0.020 |     0.093 |   0.030 | pass                 |
| wrong period (true event), A2 one-TOI                                           |   0.000 |   48 |  0.000 |     0.061 |   0.030 | pass                 |
| repeat recovery, A2 one-TOI                                                     |  49.000 |   90 |  0.544 |     0.634 |   0.600 | FAIL                 |
| alias power, A2 one-TOI: median excluded 0.9757, round-2 rule 0.9846            | nan     |   48 |  0.991 |   nan     |   0.800 | pass                 |
| wrong period, B, own events removed                                             |   3.000 |  208 |  0.014 |     0.037 |   0.030 | pass                 |
| wrong period (true event), B, own events removed                                |   0.000 |  205 |  0.000 |     0.015 |   0.030 | pass                 |
| repeat recovery, B, own events removed                                          | 203.000 |  365 |  0.556 |     0.600 |   0.600 | FAIL                 |
| alias power, B, own events removed: median excluded 0.9838, round-2 rule 0.9848 | nan     |  205 |  0.999 |   nan     |   0.800 | pass                 |
| wrong period, A2 multi-TOI                                                      |   6.000 |   14 |  0.429 |     0.675 |   0.030 | reported, not judged |
| wrong period (true event), A2 multi-TOI                                         |   0.000 |    8 |  0.000 |     0.312 |   0.030 | reported, not judged |
| repeat recovery, A2 multi-TOI                                                   |  12.000 |   26 |  0.462 |     0.638 |   0.600 | reported, not judged |
| alias power, A2 multi-TOI: median excluded 0.9241, round-2 rule 0.9756          | nan     |    8 |  0.947 |   nan     |   0.800 | reported, not judged |
| wrong period, B, all events                                                     |  54.000 |  242 |  0.223 |     0.272 |   0.030 | reported, not judged |
| wrong period (true event), B, all events                                        |   0.000 |  188 |  0.000 |     0.016 |   0.030 | reported, not judged |
| repeat recovery, B, all events                                                  | 203.000 |  365 |  0.556 |     0.600 |   0.600 | reported, not judged |
| alias power, B, all events: median excluded 0.9839, round-2 rule 0.9848         | nan     |  188 |  0.999 |   nan     |   0.800 | reported, not judged |
| false event, A1 one-TOI                                                         |   3.000 |   90 |  0.033 |     0.084 |   0.050 | pass                 |
| false event, B, own events removed                                              |  10.000 | 1072 |  0.009 |     0.016 |   0.050 | pass                 |
| false event, B, all events                                                      |  93.000 | 1072 |  0.087 |     0.102 |   0.050 | reported, not judged |
| false exclusion, A2 one-TOI                                                     |   0.000 |   90 |  0.000 |     0.033 |   0.050 | pass                 |
| false exclusion, B                                                              |   0.000 | 1072 |  0.000 |     0.003 |   0.050 | pass                 |
