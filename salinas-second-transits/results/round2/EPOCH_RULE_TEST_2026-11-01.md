# Epoch exclusion rule test, 2026-11-01

Script: `round2/s10_epoch_rule_test.py`. Full tables: `epoch_rule_test_summary.txt` and the
`epoch_rule_*.csv` files here. After this test the pipeline switched to half3s; "current" below means the
rule it used until then.

With sigma = hypot(epoch noise, reference depth error), an epoch excludes a period when its box depth is:

| Rule | Below |
|---|---|
| current | max(depth0 - 3 sigma, 0.5 depth0) |
| both | min(depth0 - 3 sigma, 0.5 depth0) |
| half3s (proposed) | 0.5 depth0 - 3 sigma |

The comparison code reproduced all 53 stored alias tables and all 39 quoted-period statuses under the
current rule before the switch, and under half3s after it.

## Confirmed planets: the main test

90 TOIs dispositioned CP or KP on ExoFOP, drawn at random before looking at their light curves: 30 in
each depth range (500-2000, 2000-5000, 5000-25000 ppm), P = 3-40 d, none in `candidates.csv`. 83 had at
least 3 clean transits over at least 300 d. T1 and T2 are the first and last clean transits, each with
15 min of scatter added. Depth, duration and starting ephemeris come from the TOI table.

| Epoch SNR | Stars | True period excluded: current | both | half3s | Median false aliases left: current | both | half3s |
|---|---|---|---|---|---|---|---|
| under 6 | 26 | 21 | 10 | 1 | 8.5 | 22 | 101 |
| 6 to 12 | 19 | 12 | 6 | 2 | 20 | 20 | 29 |
| over 12 | 38 | 35 | 2 | 1 | 21 | 21 | 21 |
| all | 83 | 68 (82%) | 18 (22%) | 4 (5%) | 16 | 20 | 30 |

The table above uses the TOI table's fitted depth as the reference, which runs deeper than a box mean, so
it overstates how often the current rule excludes the true period in the pipeline. The pipeline compares
box depth with box depth, as below.

Same test with the reference depth measured the way the pipeline does (box depth at T1):

| Epoch SNR | Stars | True period excluded: current | both | half3s | Median false aliases left: current | both | half3s |
|---|---|---|---|---|---|---|---|
| under 6 | 35 | 24 | 8 | 1 | 9 | 18 | 102 |
| 6 to 12 | 16 | 5 | 3 | 1 | 23 | 24 | 28 |
| over 12 | 32 | 17 | 2 | 1 | 20.5 | 21 | 21 |
| all | 83 | 46 (55%) | 13 (16%) | 3 (4%) | 16 | 20 | 32 |

`round2/s10_epoch_rule_test.py known --set toi --ref box`, outputs `epoch_rule_known_*_toi_boxref.csv`.
Same stars, seed and T1/T2 as above. As in `st2_pipeline.search_one`, `core.refine_event` at T1 gives the
reference depth, its error and the duration, kept only if its SNR is at least 5 and it is 0.33-3 times the
catalogue depth. That held for 63 of 83 stars (median box/TOI depth 0.94, quartiles 0.86-0.99); the other
20, all with a T1 box SNR under 5, keep the TOI values, as the pipeline would. On the 63 alone the true
period is excluded for 30 / 7 / 2 under current / both / half3s. Epoch SNR bins use each run's own
reference, so stars move between bins.

Each star has a median of 1853 aliases. False aliases exclude multiples of the true period, which the
data can't rule out.

Share of clean epochs that exclude, by epoch SNR:

| Epoch SNR | Real transit: current | both | half3s | Flat data: current | both | half3s |
|---|---|---|---|---|---|---|
| under 4 | 25% | 1% | 0% | 92% | 46% | 6% |
| 4 to 6 | 10% | 4% | 1% | 98% | 95% | 27% |
| 6 to 8 | 3% | 1% | 0% | 100% | 100% | 70% |
| 8 to 10 | 20% | 7% | 1% | 100% | 100% | 98% |
| 10 to 15 | 13% | 1% | 1% | 100% | 100% | 100% |
| 15 to 20 | 45% | 0% | 0% | 100% | 100% | 100% |
| 20 to 30 | 63% | 0% | 0% | 100% | 100% | 100% |
| over 30 | 96% | 1% | 1% | 100% | 100% | 100% |

The median real transit has a box depth of 0.73-0.93 of the TOI depth, because a box mean over the full
duration is shallower than a fitted depth. Under the current rule, the stronger the signal, the more
often that normal offset excludes the true period.

The 4 stars still excluded under half3s:

- TIC 343628284 (TOI-1448.01): the period hit the edge of the refinement grid and the stacked depth is
  0.27 of the TOI depth, so this is probably a bad ephemeris, not a rule failure.
- TIC 142276270 (TOI-1136.01), 333473672 (TOI-1422.01), 447061717 (TOI-1231.01): one flat clean epoch
  each. Not explained. Timing variations or a data problem at that epoch are possible.

## Injections

Box transits added to the cached light curves of 89 candidate stars, 141 trials per SNR, with the
reference depth error set to 5%. The share of trials where the true period was excluded:

| Epoch SNR | Other transits at full depth: current | both | half3s | At 0.7 depth: current | both | half3s |
|---|---|---|---|---|---|---|
| 4 | 29% | 9% | 4% | 67% | 32% | 6% |
| 6 | 11% | 9% | 4% | 53% | 49% | 4% |
| 8 | 9% | 6% | 1% | 66% | 41% | 4% |
| 12 | 9% | 4% | 1% | 89% | 23% | 4% |
| 20 | 6% | 1% | 0% | 99% | 6% | 1% |
| 50 | 4% | 0% | 0% | 100% | 1% | 0% |

Median false aliases left at full depth, current / both / half3s: 21 / 25 / 147 at SNR 4, 21 / 21 / 47 at
SNR 6, 20 / 21 / 27 at SNR 8, and 20 / 20 / 21 at SNR 12. Injections go in after detrending, so they don't
include any depth lost to it.

## Real candidates

53 alias events, epoch SNR 4.3 to 27.4, median 8.8. Aliases allowed or supported:

| Rule | Median | Mean |
|---|---|---|
| current | 17 | 19.4 |
| both | 17 | 19.6 |
| half3s | 21 | 27.0 |

Under half3s, 23 events don't change, the median grows by 1, and the largest growth is TIC 302377820b at
SNR 4.8, from 29 to 141. The other large increases are also at SNR under 6.5.

Quoted periods that change under half3s (`epoch_rule_candidate_periods.csv`):

| TIC | Period (d) | current | both | half3s |
|---|---|---|---|---|
| 230133616 | 23.7533, 23.7535 | excluded | supported | supported |
| 275647095 | 17.0640, 17.0643 | excluded | supported | supported |
| 241550976 | 21.1968, 21.1971 | excluded | supported | supported |
| 317132251 | 25.2103 | excluded | allowed | allowed |
| 404736082b | 8.9212 | excluded | excluded | supported |

"Supported" means a matching depth at some epoch. It doesn't confirm a period.

## Limits

- The 57 `multi` rows in `candidates.csv` turned out not to be known periods: only 3 are on the TOI list,
  and one of those (TIC 380255458, K2-140 b) lists 13.10 d against the planet's 6.57 d. Results for that
  set are in the summary but carry no weight here (35 of 57 stars usable; true period excluded under
  current / both / half3s: 24 / 11 / 7).
- The planets have shorter periods than the single-transit candidates and aren't a random sample of them.
- For the candidates, the pipeline's reference depth is the box depth re-measured on the paper event
  (median 0.83 of the paper's depth), not the paper's fitted depth. The box-reference table is the
  like-for-like comparison.
- Under half3s an epoch with SNR under about 6 can barely exclude anything (0.5 depth0 - 3 sigma is near
  zero), so a period marked supported or allowed at low epoch SNR carries little weight. 302377820b, at
  alias-event SNR 4.8, goes from 29 to 141 allowed periods.
- The common-period test doesn't use this rule, so its flags are unaffected.
