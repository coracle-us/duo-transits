# Round 2 test: independent score (ind_test.py rows, ind_score.py)

Rows: `ind_a_rows.csv`, `ind_b_rows.csv`

Scored: 123 one-TOI arm-A stars (both A1 and A2 scored), 23 multi-TOI stars, 600 of 600 arm-B trials; 292 of 292 arm-A rows scored.

## Outcome: **partial**

| Criterion | met |
|---|---|
| false_exclusion_rate | yes |
| log_period_excluded | yes |
| false_event_rate | NO |
| wrong_period_rate | NO |
| repeat_recovery_rate | yes |
| probability_calibration | yes |

## Criteria table

| Criterion | Sample | k | n | value | 95% upper | limit | verdict |
|---|---|---|---|---|---|---|---|
| false_exclusion | A2 one-TOI, overall | 0 | 123 | 0.0000 | 0.0241 | <= 0.05 | met |
| false_exclusion | A2 one-TOI, ref SNR 6-8 (under 60 trials) | 0 | 17 | 0.0000 | 0.1616 | <= 0.05 | not judged |
| false_exclusion | A2 one-TOI, ref SNR 8-12 (under 60 trials) | 0 | 15 | 0.0000 | 0.1810 | <= 0.05 | not judged |
| false_exclusion | A2 one-TOI, ref SNR 12-20 (under 60 trials) | 0 | 15 | 0.0000 | 0.1810 | <= 0.05 | not judged |
| false_exclusion | A2 one-TOI, ref SNR >20 (under 60 trials) | 0 | 41 | 0.0000 | 0.0705 | <= 0.05 | not judged |
| false_exclusion | A2 one-TOI, ref SNR under 6 or none | 0 | 35 | 0.0000 | 0.0820 | <= 0.05 | not judged |
| false_exclusion | B all trials, overall | 0 | 600 | 0.0000 | 0.0050 | <= 0.05 | met |
| false_exclusion | B all trials, ref SNR 6-8 (under 60 trials) | 0 | 46 | 0.0000 | 0.0630 | <= 0.05 | not judged |
| false_exclusion | B all trials, ref SNR 8-12 (under 60 trials) | 0 | 52 | 0.0000 | 0.0560 | <= 0.05 | not judged |
| false_exclusion | B all trials, ref SNR 12-20 | 0 | 71 | 0.0000 | 0.0413 | <= 0.05 | met |
| false_exclusion | B all trials, ref SNR >20 | 0 | 254 | 0.0000 | 0.0117 | <= 0.05 | met |
| false_exclusion | B all trials, ref SNR under 6 or none | 0 | 177 | 0.0000 | 0.0168 | <= 0.05 | not judged |
| log_period_excluded | A1 one-TOI, ref SNR >= 8, median |  | 71 | 0.6541 |  | >= 0.5 | met |
| log_period_excluded | B all trials, ref SNR >= 8, median |  | 377 | 0.6702 |  | >= 0.5 | met |
| false_event | A1 one-TOI stars | 7 | 123 | 0.0569 | 0.1042 | <= 0.05 | MISSED |
| false_event | B, own events removed | 11 | 600 | 0.0183 | 0.0302 | <= 0.05 | met |
| false_event | B, all events (reported beside) | 84 | 600 | 0.1400 | 0.1655 | <= 0.05 | not judged |
| false_event | B null runs with a passing event (reported, not judged) | 85 | 600 | 0.1417 | 0.1672 | <= 0.05 | not judged |
| wrong_period | A2 one-TOI | 6 | 76 | 0.0789 | 0.1499 | <= 0.03 | MISSED |
| wrong_period | B, own events removed | 6 | 151 | 0.0397 | 0.0769 | <= 0.03 | MISSED |
| wrong_period | B, all events (reported beside) | 46 | 178 | 0.2584 | 0.3180 | <= 0.03 | not judged |
| repeat_recovery | A2 one-TOI | 81 | 123 | 0.6585 |  | >= 0.6 | met |
| repeat_recovery | B, own events removed | 150 | 230 | 0.6522 |  | >= 0.6 | met |
| repeat_recovery | B, all events (reported beside) | 150 | 230 | 0.6522 |  | >= 0.6 | not judged |
| probability_calibration | A2 one-TOI + B: observed 0, predicted 1.63, 95% Poisson interval 0-5 |  | 723 | 0.0000 |  | in 0-5 | met |
| multi-TOI (not judged) | false exclusion, A2 | 0 | 23 | 0.0000 | 0.1221 | <= 0.05 | not judged |
| multi-TOI (not judged) | false event, A1 | 5 | 23 | 0.2174 | 0.4039 | <= 0.05 | not judged |
| multi-TOI (not judged) | wrong period, A2 | 4 | 11 | 0.3636 | 0.6502 | <= 0.03 | not judged |
| multi-TOI (not judged) | repeat recovery, A2 | 10 | 23 | 0.4348 |  | >= 0.6 | not judged |
| multi-TOI (not judged) | ln P excluded, A1, ref SNR >= 8, median |  | 9 | 0.6828 |  | >= 0.5 | not judged |

'95% upper' is the one-sided 95% Clopper-Pearson bound (reported, not part of any criterion).

## Probability calibration by R5 bin

| R5 bin | p | trials (A2 one-TOI + B) | predicted | observed |
|---|---|---|---|---|
| unmeasured | 0.0010 | 198 (27 + 171) | 0.195 | 0 |
| 5-8 | 0.0024 | 80 (25 + 55) | 0.189 | 0 |
| 8-12 | 0.0025 | 67 (15 + 52) | 0.165 | 0 |
| 12-20 | 0.0026 | 85 (15 + 70) | 0.220 | 0 |
| >20 | 0.0029 | 293 (41 + 252) | 0.858 | 0 |
| total | | 723 | 1.627 | 0 |
