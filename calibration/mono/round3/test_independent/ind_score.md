# Round 3 test: independent score

From `ind_a_rows.csv` and `ind_b_rows.csv` (`ind_test3.py`), scored by `ind_score.py`, written from
`prereg_round3.yaml` only. Made before the primary rows and score were read.

Scored: 90 one-TOI arm-A stars (A1 and A2 both run; 26 multi-TOI stars reported, not judged), 1072 arm-B trials on 269 hosts. Failed runs: 0 arm-A rows, 0 arm-B rows.

## Outcome: **NOT DEFINED by the pre-registration**

The outcomes section of prereg_round3.yaml does not define this combination: missed = repeat_recovery_rate; wrong_period_rate_true_event and alias_power are met, so it is not `fail`; not every criterion is met, so it is not `pass`; `calls_after_vetting` requires repeat_recovery_rate and false_exclusion_rate to be met, so it is not that either.

| Criterion | Limit | Verdict |
|---|---|---|
| wrong_period_rate | <= 3% (A2 needs >= 40 called) | met |
| wrong_period_rate_true_event | <= 3% (same minimum) | met |
| false_event_rate | <= 5% | met |
| repeat_recovery_rate | >= 60% | MISSED |
| alias_power | round-3 median >= 0.8 x round-2 median | met |
| false_exclusion_rate | <= 5% | met |

## Per arm

| Criterion | Arm | k | n | rate | verdict | note |
|---|---|---|---|---|---|---|
| wrong_period_rate | A2 one-TOI | 1 | 49 | 2.0% | met |  |
| wrong_period_rate | B own removed | 3 | 208 | 1.4% | met |  |
| wrong_period_rate | B all events | 54 | 242 | 22.3% | reported |  |
| wrong_period_rate | A2 multi-TOI | 6 | 14 | 42.9% | reported |  |
| wrong_period_rate_true_event | A2 one-TOI | 0 | 48 | 0.0% | met |  |
| wrong_period_rate_true_event | B own removed | 0 | 205 | 0.0% | met |  |
| wrong_period_rate_true_event | B all events | 0 | 188 | 0.0% | reported |  |
| wrong_period_rate_true_event | A2 multi-TOI | 0 | 8 | 0.0% | reported |  |
| false_event_rate | A1 one-TOI | 3 | 90 | 3.3% | met | stars with any passing event |
| false_event_rate | A1 one-TOI, passing event not at a true transit | 3 | 90 | 3.3% | reported |  |
| false_event_rate | B own removed | 10 | 1072 | 0.9% | met |  |
| false_event_rate | B all events | 93 | 1072 | 8.7% | reported |  |
| false_event_rate | A1 multi-TOI | 8 | 26 | 30.8% | reported |  |
| false_event_rate | A2 one-TOI (not a criterion) | 4 | 90 | 4.4% | reported |  |
| null runs with a passing event | B | 100 | 1072 | 9.3% | reported |  |
| repeat_recovery_rate | A2 one-TOI | 49 | 90 | 54.4% | MISSED |  |
| repeat_recovery_rate | B own removed | 203 | 365 | 55.6% | MISSED |  |
| repeat_recovery_rate | B all events | 203 | 365 | 55.6% | reported |  |
| repeat_recovery_rate | A2 multi-TOI | 12 | 26 | 46.2% | reported |  |
| false_exclusion_rate | A2 one-TOI | 0 | 90 | 0.0% | met |  |
| false_exclusion_rate | B | 0 | 1072 | 0.0% | met |  |
| false_exclusion_rate | A2 multi-TOI | 0 | 26 | 0.0% | reported |  |

## Alias power (called trials whose alias event is a true transit)

| Arm | n | median fraction excluded, round 3 | round 2 rule, same trials | ratio | verdict |
|---|---|---|---|---|---|
| A2 one-TOI | 48 | 0.9757 | 0.9846 | 0.991 | met |
| B own removed | 205 | 0.9838 | 0.9848 | 0.999 | met |
| A2 all stars | 56 | 0.9751 | 0.9846 | 0.990 | reported |
| B all events | 188 | 0.9839 | 0.9848 | 0.999 | reported |

## Calls

| Arm | resolved | narrowed | all excluded | no repeat |
|---|---|---|---|---|
| A1 one-TOI | 0 | 3 | 0 | 87 |
| A2 one-TOI | 5 | 44 | 0 | 41 |
| B own removed | 17 | 191 | 3 | 861 |
| B all events | 21 | 221 | 22 | 808 |

## Trials behind every miss

**A2 one-TOI wrong calls** (1)

| tic | toi | ref_snr | truth_P | t2 | t2_snr | call | n_alias_ok | n_alias | event_true |
|---|---|---|---|---|---|---|---|---|---|
| 192790476 | 2459 | 9.911 | 19.1 | 4024 | 14.06 | narrowed | 49 | 2571 | False |

**B (own removed) wrong calls** (3)

| tic | draw | ref_snr | truth_P | t2@ | t2_snr@ | call@ | n_alias_ok@ | n_alias@ | event_true@ |
|---|---|---|---|---|---|---|---|---|---|
| 339399841 | 1 | 5.614 | 751.2 | 4196 | 14.34 | narrowed | 16 | 1064 | False |
| 143355733 | 1 | 8.093 | 744.4 | 3146 | 17.01 | narrowed | 9 | 1038 | False |
| 139695717 | 1 | 6.51 | 157.5 | 3746 | 8.213 | narrowed | 1329 | 1510 | False |

**A1 one-TOI stars with a passing event** (3)

| tic | toi | ref_snr | truth_P | n_vetted | t2 | t2_snr | call |
|---|---|---|---|---|---|---|---|
| 364395234 | 1232 | 11.45 | 14.25 | 34 | 1483 | 11.87 | narrowed |
| 237232044 | 1443 | 6.752 | 23.54 | 1 | 3486 | 17.21 | narrowed |
| 192790476 | 2459 | 9.911 | 19.1 | 3 | 4024 | 14.05 | narrowed |

**B (own removed) trials with a false event** (10)

| tic | draw | ref_snr | truth_P | t2@ | t2_snr@ | call@ | n_alias_ok@ | n_alias@ | event_true@ |
|---|---|---|---|---|---|---|---|---|---|
| 325081244 | 1 | 535.8 | 97.96 | 4202 | 655 | narrowed | 11 | 2155 | True |
| 155157963 | 2 | 11.73 | 22.82 | 3818 | 10.14 | all excluded | 0 | 0 | False |
| 164365027 | 0 | 15.97 | 20.67 | 2748 | 16.13 | narrowed | 13 | 744 | True |
| 278990954 | 2 | 136.3 | 315.8 | 3042 | 124.4 | narrowed | 3 | 1579 | True |
| 209217764 | 3 | 29.92 | 217.3 | 3806 | 14.6 | all excluded | 0 | 67 | False |
| 388909695 | 0 | 45.24 | 731 | 3547 | 145.8 | narrowed | 7 | 730 | True |
| 339399841 | 1 | 5.614 | 751.2 | 4196 | 14.34 | narrowed | 16 | 1064 | False |
| 143355733 | 1 | 8.093 | 744.4 | 3146 | 17.01 | narrowed | 9 | 1038 | False |
| 189568918 | 2 | 36.97 | 48.15 | 4068 | 45.62 | narrowed | 49 | 2503 | True |
| 139695717 | 1 | 6.51 | 157.5 | 3746 | 8.213 | narrowed | 1329 | 1510 | False |

**False exclusions, A2 one-TOI** (0)

**False exclusions, B** (0)

**False exclusions, A2 multi-TOI (not judged)** (0)

