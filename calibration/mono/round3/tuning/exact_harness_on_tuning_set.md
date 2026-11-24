# Round 3 test: score (pre-registration 1.0)

Arm A: 199 one-TOI stars scored (255 stars in all). Arm B: 1115 trials on 279 hosts.

**Outcome: pass**

| criterion                                                                       |       k |    n |   rate |   upper95 |   limit | verdict              |
|:--------------------------------------------------------------------------------|--------:|-----:|-------:|----------:|--------:|:---------------------|
| wrong period, A2 one-TOI                                                        |   2.000 |  133 |  0.015 |     0.047 |   0.030 | pass                 |
| wrong period (true event), A2 one-TOI                                           |   0.000 |  131 |  0.000 |     0.023 |   0.030 | pass                 |
| repeat recovery, A2 one-TOI                                                     | 133.000 |  199 |  0.668 |     0.724 |   0.600 | pass                 |
| alias power, A2 one-TOI: median excluded 0.9731, round-2 rule 0.9807            | nan     |  131 |  0.992 |   nan     |   0.800 | pass                 |
| wrong period, B, own events removed                                             |   3.000 |  260 |  0.012 |     0.030 |   0.030 | pass                 |
| wrong period (true event), B, own events removed                                |   0.000 |  257 |  0.000 |     0.012 |   0.030 | pass                 |
| repeat recovery, B, own events removed                                          | 259.000 |  420 |  0.617 |     0.656 |   0.600 | pass                 |
| alias power, B, own events removed: median excluded 0.9844, round-2 rule 0.9848 | nan     |  257 |  1.000 |   nan     |   0.800 | pass                 |
| wrong period, A2 multi-TOI                                                      |   6.000 |   31 |  0.194 |     0.347 |   0.030 | reported, not judged |
| wrong period (true event), A2 multi-TOI                                         |   1.000 |   26 |  0.038 |     0.170 |   0.030 | reported, not judged |
| repeat recovery, A2 multi-TOI                                                   |  31.000 |   56 |  0.554 |     0.667 |   0.600 | reported, not judged |
| alias power, A2 multi-TOI: median excluded 0.9682, round-2 rule 0.9830          | nan     |   26 |  0.985 |   nan     |   0.800 | reported, not judged |
| wrong period, B, all events                                                     |  65.000 |  309 |  0.210 |     0.252 |   0.030 | reported, not judged |
| wrong period (true event), B, all events                                        |   0.000 |  244 |  0.000 |     0.012 |   0.030 | reported, not judged |
| repeat recovery, B, all events                                                  | 259.000 |  420 |  0.617 |     0.656 |   0.600 | reported, not judged |
| alias power, B, all events: median excluded 0.9842, round-2 rule 0.9846         | nan     |  244 |  1.000 |   nan     |   0.800 | reported, not judged |
| false event, A1 one-TOI                                                         |   6.000 |  199 |  0.030 |     0.059 |   0.050 | pass                 |
| false event, B, own events removed                                              |  13.000 | 1115 |  0.012 |     0.018 |   0.050 | pass                 |
| false event, B, all events                                                      | 134.000 | 1115 |  0.120 |     0.137 |   0.050 | reported, not judged |
| false exclusion, A2 one-TOI                                                     |   0.000 |  199 |  0.000 |     0.015 |   0.050 | pass                 |
| false exclusion, B                                                              |   0.000 | 1115 |  0.000 |     0.003 |   0.050 | pass                 |
