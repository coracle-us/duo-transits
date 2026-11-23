# Round 3 tuning log

Tuning used only stars already seen (pre-registration 1.0, `tuning.data`): the 255 arm-A stars and
279 arm-B hosts of round 2's tuning and test sets, 1,625 trials, stored by `dump3.py`; and the same
arm-B trials with a 10% error on the catalogued depth and duration (`--noise 0.1`, seed 20261205).
"A test2" and "B test2" below are round 2's test stars; they are tuning data now.

Variant tables: `tuning/base.csv`, `t1.csv`, `t1ref.csv`, `t5.csv`, `t3.csv`, `final_noisy.csv`
(`tune3.py <family>`). Exact numbers for the proposed rule: `tuning/exact_harness_on_tuning_set.md`
(`harness3.py` on every tuning trial, arm B with the 10% error).

## Baseline

`tune3.py base` with every round-3 switch off reproduces round 2 on its own test set: false events
7 of 123, wrong period 6 of 76 and 6 of 151, recovery 81 of 123 and 150 of 230.

## T1, period tolerance (`t1.csv`, `t1ref.csv`; 24 variants each)

Growth factor c of 1, 2, 3; cap of 1, 2, 4 durations or none; timing error fixed at 15 min or also
0.5 duration/SNR; catalogued or fitted reference time.

- Wrong periods with a real second transit, all tuning trials: 5 of 127 (arm A) and 1 of 260 (arm
  B) under round 2; 0 of 140 and 0 of 265 with c = 2 and a cap of 1 or 2 durations.
- c = 1 leaves 1 in arm B. A cap of 4 or none leaves 1 in arm A (an epoch whose transit lies
  further off than the window predicts still excludes; a cap makes such far epochs unconstrained).
- The duration/SNR term changes nothing. The fitted reference time is not better (1 of 140).
- Alias power, median share of aliases excluded: 0.9747 to 0.9701 (arm A), 0.9842 to 0.9838 (arm B)
  for c = 2, cap 2. Cap 1 costs more (0.9663).
- Chosen: c = 2, cap 2 durations, 15 min per transit time, catalogued reference time.

## T3, minimum period (`t3.csv`)

Off, catalogued duration, half the limit: no change in any wrong-period count. Left as round 2.

## T5, event vetting (`t5.csv`; 48 variants on top of the chosen T1)

- `sec_box_snr` (event depth over the scatter of same-duration boxes in its sector) at 8, 9, 10:
  arm-A false events 9, 7, 6, 6 of 199 (none, 8, 9, 10); arm-B wrong periods 7, 5, 4, 4;
  recovery in arm B 270, 262, 259, 250 of 418.
- Common-period veto: removes real repeats (arm-A recovery 142 to 55 of 199, because the same
  drift breaks the common period). Not adopted.
- Higher event SNR (10 or 12) when the reference SNR is below 7: one or two fewer false events,
  9 to 13 fewer recoveries in arm A. Not adopted.
- Chosen: `sec_box_snr >= 9`.

## Not tried

T2 (refitting the period), the two-product agreement check of T5. T4 changes labels only
("resolved" or "narrowed"), not whether the truth is in the alias set, so it cannot move a criterion.

## The proposed rule on the tuning set (exact harness; not a result)

| Criterion | Arm A, one-TOI | Arm B, own removed, 10% catalogue error | Limit |
|---|---|---|---|
| Wrong period | 2 of 133 (1.5%) | 3 of 260 (1.2%) | <= 3% |
| Wrong period, true event | 0 of 131 | 0 of 257 | <= 3% |
| False event | 6 of 199 (3.0%) | 13 of 1115 (1.2%) | <= 5% |
| Repeat recovery | 133 of 199 (67%) | 259 of 420 (62%) | >= 60% |
| Alias power against round 2 | 0.992 | 1.000 | >= 0.8 |
| False exclusion | 0 of 199 | 0 of 1115 | <= 5% |

## Margins (told to Jon before the freeze)

- Arm-B recovery is 62% against a 60% guard. With about 380 test trials that have a clean repeat
  its standard error is about 2.5 points.
- Arm A in the test has about 98 one-TOI stars and about 55 called trials. One wrong call passes,
  two do not; four false-event stars pass, five do not. At the tuning-set rates the chance of
  missing each is roughly 1 in 5.
- The wrong periods left are false events, which T1 cannot touch.
- The 13 new arm-A stars above 100 d have few transits; some may have under 3 clean ones and be
  skipped, as round 1's rule says.
