# Round 2 test: independent re-derivation compared with the primary run

Independent rows: `ind_a_rows.csv`, `ind_b_rows.csv` (`ind_test.py`); independent score: `ind_score.md`
(`ind_score.py`); ambiguities: `AMBIGUITIES.md`. Their SHA-256 sums were written to
`SHA256_before_compare.txt` before `../test/arm_a_rows.csv`, `arm_b_rows.csv`, `score_summary.md`
and `score_rates.csv` were opened. This table is printed by `compare.py`.

## Result

- **Per trial: no disagreement.** All 292 arm-A and all 600 arm-B trials agree in every compared
  column; the numeric columns agree to the last digit (largest absolute difference 0).
- **Criteria: no disagreement.** All 23 rate rows have the same k, n, rate and verdict. Both
  outcomes are **partial** (false event in A1 and wrong period in A2 and in B missed).
- One difference in a reported number, not in a verdict: predicted false-exclusion count 1.6282
  (primary) against 1.6272 (independent). The independent scorer uses p = (k + 0.5)/(n + 1)
  unrounded; 1.6282 is what the four-decimal p of the FROZEN_RULE table (and FREEZE.json) gives.
  Both round to 1.63 and give the interval 0-5. Inferred from the numbers; `score2.py` was not read.

## Per-column agreement

Booleans and calls compared exactly; numbers with relative tolerance 1e-9. Empty against empty
counts as agreement (t2 where there is no alias event).

### Arm A: 292 trials, matched one to one on ['tic', 'toi', 'arm']; primary status: {'ok': 292}

| column | agree | disagree | max abs difference |
|---|---|---|---|
| truth_excluded | 292 | 0 |  |
| false_event | 292 | 0 |  |
| repeat_detected | 292 | 0 |  |
| truth_in_alias_set | 292 | 0 |  |
| ref_measured | 292 | 0 |  |
| call | 292 | 0 |  |
| frac_lnP_excluded | 292 | 0 | 0 |
| ref_snr | 292 | 0 | 0 |
| n_clean_repeats | 292 | 0 | 0 |
| n_vetted | 292 | 0 | 0 |
| n_alias_ok | 292 | 0 | 0 |
| t2 | 292 | 0 | 0 |
| depth_ref | 292 | 0 | 0 |
| dur_ref_d | 292 | 0 | 0 |

### Arm B: 600 trials, matched one to one on ['tic', 'draw']; primary status: {'ok': 600}

| column | agree | disagree | max abs difference |
|---|---|---|---|
| truth_excluded | 600 | 0 |  |
| false_event | 600 | 0 |  |
| repeat_detected | 600 | 0 |  |
| truth_in_alias_set | 600 | 0 |  |
| ref_measured | 600 | 0 |  |
| call | 600 | 0 |  |
| frac_lnP_excluded | 600 | 0 | 0 |
| ref_snr | 600 | 0 | 0 |
| n_clean_repeats | 600 | 0 | 0 |
| n_vetted | 600 | 0 | 0 |
| n_alias_ok | 600 | 0 | 0 |
| t2 | 600 | 0 | 0 |
| depth_ref | 600 | 0 | 0 |
| dur_ref_d | 600 | 0 | 0 |
| false_event_own | 600 | 0 |  |
| repeat_detected_own | 600 | 0 |  |
| call_own | 600 | 0 |  |
| truth_in_alias_set_own | 600 | 0 |  |
| n_vetted_own | 600 | 0 | 0 |
| n_alias_ok_own | 600 | 0 | 0 |
| t2_own | 600 | 0 | 0 |
| null_event | 600 | 0 |  |
| n_null_events | 600 | 0 | 0 |

### Disagreeing trials

None.

Also equal, beyond the columns asked for: n_vetted, n_alias_ok, t2, depth_ref, dur_ref_d,
n_clean_repeats, null_event, n_null_events, and the number of sectors per star (`n_sectors`).

## Criteria tables side by side

| criterion | primary k | primary n | primary rate | primary verdict | independent k | independent n | independent rate | independent verdict | same |
|---|---|---|---|---|---|---|---|---|---|
| false exclusion, A2 one-TOI, all | 0 | 123 | 0.0000 | pass | 0 | 123 | 0.0000 | met | yes |
| false exclusion, A2 one-TOI, SNR 6-8 | 0 | 17 | 0.0000 | reported, not judged | 0 | 17 | 0.0000 | not judged | yes |
| false exclusion, A2 one-TOI, SNR 8-12 | 0 | 15 | 0.0000 | reported, not judged | 0 | 15 | 0.0000 | not judged | yes |
| false exclusion, A2 one-TOI, SNR 12-20 | 0 | 15 | 0.0000 | reported, not judged | 0 | 15 | 0.0000 | not judged | yes |
| false exclusion, A2 one-TOI, SNR >20 | 0 | 41 | 0.0000 | reported, not judged | 0 | 41 | 0.0000 | not judged | yes |
| false exclusion, B, all | 0 | 600 | 0.0000 | pass | 0 | 600 | 0.0000 | met | yes |
| false exclusion, B, SNR 6-8 | 0 | 46 | 0.0000 | reported, not judged | 0 | 46 | 0.0000 | not judged | yes |
| false exclusion, B, SNR 8-12 | 0 | 52 | 0.0000 | reported, not judged | 0 | 52 | 0.0000 | not judged | yes |
| false exclusion, B, SNR 12-20 | 0 | 71 | 0.0000 | pass | 0 | 71 | 0.0000 | met | yes |
| false exclusion, B, SNR >20 | 0 | 254 | 0.0000 | pass | 0 | 254 | 0.0000 | met | yes |
| false exclusion, A2 multi-TOI | 0 | 23 | 0.0000 | reported, not judged | 0 | 23 | 0.0000 | not judged | yes |
| median ln P excluded, A1 one-TOI, SNR >= 8 |  | 71 | 0.6541 | pass |  | 71 | 0.6541 | met | yes |
| median ln P excluded, B, SNR >= 8 |  | 377 | 0.6702 | pass |  | 377 | 0.6702 | met | yes |
| false event, A1 one-TOI | 7 | 123 | 0.0569 | FAIL | 7 | 123 | 0.0569 | MISSED | yes |
| false event, B, own events removed | 11 | 600 | 0.0183 | pass | 11 | 600 | 0.0183 | met | yes |
| false event, B, all events | 84 | 600 | 0.1400 | reported, not judged | 84 | 600 | 0.1400 | not judged | yes |
| passing event in B null runs | 85 | 600 | 0.1417 | reported, not judged | 85 | 600 | 0.1417 | not judged | yes |
| wrong period, A2 one-TOI | 6 | 76 | 0.0789 | FAIL | 6 | 76 | 0.0789 | MISSED | yes |
| repeat recovery, A2 one-TOI | 81 | 123 | 0.6585 | pass | 81 | 123 | 0.6585 | met | yes |
| wrong period, B, own events removed | 6 | 151 | 0.0397 | FAIL | 6 | 151 | 0.0397 | MISSED | yes |
| repeat recovery, B, own events removed | 150 | 230 | 0.6522 | pass | 150 | 230 | 0.6522 | met | yes |
| wrong period, B, all events | 46 | 178 | 0.2584 | reported, not judged | 46 | 178 | 0.2584 | not judged | yes |
| repeat recovery, B, all events | 150 | 230 | 0.6522 | reported, not judged | 150 | 230 | 0.6522 | not judged | yes |
| probability calibration | observed 0, predicted 1.63, 95% interval 0-5 (predicted 1.6282) | 723 | | pass | observed 0, predicted 1.63, 95% Poisson interval 0-5 (predicted 1.6272) | 723 | | met | yes (see note) |

Outcome: primary **partial**, independent **partial**. 123 one-TOI arm-A stars and 600 arm-B
trials scored in both (thresholds 80 and 400).

## Causes of disagreement

None to trace. What the agreement does and does not show is in the next section.

## What this agreement does not verify

- Both harnesses call the same pipeline functions (`st2core`, `st2_pipeline`) on the same cached
  light curves, and `ind_test.py` descends from `ind2.py`, which had already been reconciled with
  the primary harness on 20 tuning stars. A fault in the pipeline, in FROZEN_RULE.md itself, or in
  a reading of it shared by both would not show here. Bit-level equality reflects the shared code.
- The trial definitions were taken as given: the sample draws (seeds, pools, exclusions), the
  refined "true" ephemerides of arm A, the injection parameters and limb-darkening coefficients,
  the 15-min time errors and the `multi` flags were not re-derived.
- The R5 table (tuning-set k and n) was taken from FROZEN_RULE.md, not recomputed.
- Light-curve product lists come from a live MAST query in both runs; that they returned the same
  sectors is shown only by the equal `n_sectors`, not by a file-level comparison.
- Whether the freeze preceded the loading of test stars was not checked.

## The depth-error choice (AMBIGUITIES 1-2) on the test set

`ind2.py` used the reference depth error in the alias table; FROZEN_RULE says the event's own.
With the reference error instead: `call` changes in 2 arm-A trials and 3 arm-B trials (0 with own
events removed), `truth_in_alias_set` in none, `n_alias_ok` in 36 (A), 37 (B), 15 (B, own
removed). The primary rows match the event-error reading in every trial.

## Findings about the rule, the spec and the data (not disagreements)

1. **Five of the six A2 wrong periods are true second transits whose true alias the alias table
   excludes.** In TIC 232967440, 231736113, 21744120, 23769326 and 203214138 the alias event is
   within 0.2 h of a true transit and |dT - n P_true| is far below dur/2, but alias n is `excluded`:
   P_n = dT/n differs from P_true by the timing error of the two events divided by n, and over
   hundreds of cycles the predicted epochs drift 3 to 13 h off the true transits (durations here
   are 1.5 to 4 h), onto flat data, where they "exclude". R3's half-duration tolerance does not
   grow with the cycle count. Four of the five are then called `resolved` at n = 1 (for example
   P = 14.13 d for a 7.06 d planet). The same mechanism gives arm-B TIC 141522677 draw 3 (true
   alias n = 2 excluded by one epoch that has drifted 14.6 h). These are false exclusions of the
   true period inside the alias table; the false-exclusion criterion (epochs at the exact true
   period) does not count them, and the outcome `partial` ("second-transit calls need vetting by a
   person") suggests the fault is in the events, which it is not for these six. The sixth A2 case
   (TIC 123482865) is a real false event 62 h from a true epoch.
2. **Arm-B TIC 82680046 draw 1**: the alias event is the injected transit recovered 11.8 h from its
   true time with dur/2 = 10.1 h (duration 0.84 d), so it scores as a false event and a wrong
   period; and the true alias (62.1 d) is labelled `allowed_below_pmin`, so the trial is
   "resolved" at n = 1. `core.pmin_from_duration` is fed the box-search duration, which
   `refine_event` can stretch to twice the catalogued one; the pipeline's docstring calls pmin a
   soft limit, the call rule treats it as hard.
3. **The other four arm-B wrong periods** (265465927 draw 1; 76804724 draws 0 and 3; 271693283
   draw 3) are false events far from any injected transit that the null run did not flag.
   TIC 76804724 and 265465927 each carry two of the eleven own-removed false events.
4. **Round 1's pre-registration and the pipeline disagree on light-curve choice** ("QLP where
   available, else TESS-SPOC" against SPOC, then TESS-SPOC, then QLP in `best_per_sector`).
   FROZEN_RULE says "as in round 1". Both harnesses follow the pipeline (sector products used on
   the 296 test stars: 1684 SPOC, 304 TESS-SPOC, 396 QLP).
5. **The pre-registration has no outcome** for "false exclusion and ln P met, probability
   calibration missed". Not reached here.
6. **Judged SNR bins.** Only B 12-20 (71 trials) and B over 20 (254) reach 60 trials; no arm-A bin
   does (17, 15, 15, 41), and B 6-8 (46) and 8-12 (52) do not. The false-exclusion criterion below
   SNR 12 therefore rests on the overall rates alone. FROZEN_RULE predicted only the over-20 bin
   would be judged.
7. **Zero false exclusions in 723 trials** makes the calibration test weak: with 1.63 predicted,
   any count from 0 to 5 passes.
8. **Arm A reaches only 38.5 d** (pool 5 to 100 d) and has 123 one-TOI stars, 46 rather than 50
   in the 2000-5000 ppm bin; 8 arm-A rows have no stellar mass and 2 no radius, so no alias of
   theirs can be `allowed_below_pmin`.
9. **Arm B's catalogued depth is exactly inj_k squared** with no depth error (20% assumed when
   the reference is unmeasured), and the catalogued duration carries no error; only the time has
   the 15-min scatter. Real catalogue depths and durations are less exact than that.
