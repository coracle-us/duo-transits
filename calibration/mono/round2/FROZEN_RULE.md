# Round 2: the frozen rule

**Status: FROZEN. Jon approved step 2 on 2026-11-23 (`../APPROVALS.md`), before any test-set star
was loaded. Nothing in this file, the code or the samples listed in `FREEZE.json` may change.** Chosen on the
tuning set only (`TUNING_LOG.md`). In code: `harness.FROZEN` (a `rule2.Rule`; the name is what it
will be once approved).

Everything not listed here is as in round 1 (pre-registration 1.1 of round 1, with its three
harness corrections): light-curve choice, 10-min binning, detrending, the reference-depth rule, the
event search, the clean-epoch rule, the period grid.

## Definitions used below

- `dur`, `depth0`, `depth0_err`: the reference duration, depth and depth error of the trial.
  The reference is *measured* if `core.refine_event` at the catalogued time gives SNR >= 5 and a
  depth strictly between 0.33 and 3 times the catalogued depth; then they are its duration, depth and depth
  error. Otherwise they are the catalogued duration and depth, with the catalogued depth error or,
  if there is none, 20% of the depth.
- `T, F`: all detrended data of the star, sorted. `sig = core.duration_scale_noise(T, F, dur)`.
- A *box* at time c: `core.box_stats(T, F, c, dur)` gives its depth and coverage `frac`; its noise
  is `err = sig / sqrt(max(frac, 0.001))`; its edge distance is the distance from c to the nearest
  of `core.gap_edges(T)`. Its threshold is `thr = 0.5 depth0 - 3 hypot(err, depth0_err)`.
- A box is *clean* if `frac >= 0.8` and edge distance `>= 0.5 d`. It is *covered* if `frac >= 0.5`.

## Exclusion (R1, R2, R3)

A predicted epoch at time c **excludes** if all of these hold:

1. the reference is measured (R2: an unmeasured reference excludes nothing);
2. the box at c is clean, and `-3 err <= depth < thr` (R1: the 3-sigma floor);
3. for every offset j = -w .. w in steps of dur/8, the box at c + j dur/8, if covered, has
   `depth < thr` with its own err and thr (R3). w = 4 (half a duration); w = 8 (one duration) on a
   star with more than one planet candidate (more than one TOI in arm A).

Uses: (i) the true period is falsely excluded if any epoch t_cat + k P_true (k != 0) inside the data
span excludes; (ii) ln P excluded: periods log-spaced from 1 to 1000 d with step
dur / (4 max(baseline, dur)) in ln P, baseline = max(T) - min(T) (as round 1's code); a period is excluded if any
epoch t_cat + kP (k != 0), rounded to a time grid of step dur/8 starting at min(T), is an excluding
grid point; (iii) alias tables, below.

## Second events (R4)

Events are the pipeline's (`pl.sector_events` per detrended sector, top 3 per sector; in the
reference sector an event within one `dur` of the catalogued time is dropped). An event is used
only if all of these hold:

- SNR >= 8 (round 1: 7);
- depth ratio to `depth0` between 1/1.7 and 1.7 (round 1: 0.5 to 2);
- duration ratio to `dur` between 0.5 and 2 (unchanged);
- SNR >= 1.5 times the sector's strongest brightening, `anti_snr` (unchanged);
- the event is clean: `frac >= 0.8` and `edge_dist_d >= 0.5` (new, pre-registration 1.1 (a)).

The alias event is the passing event with the highest SNR.

## Calls

With an alias event at t2: t1, t2 = the catalogued time and the event time in order, dT = t2 - t1.
For n = 1, 2, ... while P_n = dT / n >= 1 d (n <= 3000): the epochs are t1 + k P_n for
k from floor((min(T) - t1) / P_n) - 1 to ceil((max(T) - t1) / P_n) + 1, without k = 0 and k = n.
The alias is `excluded` if any epoch excludes; else `supported` if any epoch's box is clean with
`|depth - depth0| < 3 hypot(err, depth0_err)`; else `allowed_below_pmin` if P_n is below
`core.pmin_from_duration(dur, mstar, rstar, depth0)`; else `allowed`.
Call: `resolved` if exactly one alias is allowed or supported, `narrowed` if more than one,
`all excluded` if none; `no repeat` without an alias event.
In the alias table only, `depth0_err` is the alias event's own depth error, not the reference's
(as `st2_pipeline.step_aliases`; round 1's third harness correction). The independent check of
2026-11-22 found this sentence missing.
The truth is in the alias set if, for n = round(dT / P_true), 1 <= n <= number of aliases,
`|dT - n P_true| < dur / 2`, and alias n is not `excluded`.

An event is *at a true transit* if it is within dur/2 of an epoch truth_t0 + k P_true other than
the reference epoch (epochs further than dur/2 from the catalogued time). A trial has a *false
event* if any passing event is not at a true transit; a *repeat is detected* if any passing event
is at one. `n_clean_repeats` counts true epochs, other than the reference, that are clean in T.

A false event is judged on A1 stars and on arm-B trials only, as the pre-registration says.

## Arm B: the host's own events (pre-registration 1.1 (b), (c))

For each draw the host's light curves (catalogued events masked by +-0.75 d, nothing injected) go
through the same reference rule, detrending and event search with the draw's catalogued time,
duration and depth. Its events, of any SNR, other than those within one catalogued duration of the
catalogued time, are the host's own. In the injected trial an event is the host's own if a null-run
event lies within one `dur` of it and it is not at a true transit. The event criteria are judged
with those events removed before vetting; the numbers without removal are reported beside them.
The null run is a scoring device only. It is not a veto in the rule: a real candidate has no
null run, so the R4 candidate "the same event appears with nothing injected" was not adopted.

## Probability (R5)

Every excluded period carries the false-exclusion probability of its trial's reference-SNR bin,
from the tuning set under this rule (one-TOI A2, B and Bx trials), p = (k + 0.5) / (n + 1):

| Reference SNR | k | n | p | 95% upper bound |
|---|---|---|---|---|
| unmeasured (under 5) | 0 | 507 | 0.0010 | 0.0059 |
| 5-8 | 0 | 211 | 0.0024 | 0.0141 |
| 8-12 | 0 | 202 | 0.0025 | 0.0147 |
| 12-20 | 0 | 192 | 0.0026 | 0.0155 |
| over 20 | 1 | 511 | 0.0029 | 0.0092 |

The pre-registration allowed bins of reference SNR and cadence. Only SNR is used: with one false
exclusion in 1,623 tuning trials there is nothing to split by cadence. These bins are not the bins
of the pass criteria (6-8, 8-12, 12-20, over 20), which stay as approved.

Calibration test: the predicted count is the sum of p over the scored A2 one-TOI and B trials; the
observed count of false exclusions must lie in the central 95% Poisson interval of that mean.

## Tuning-set numbers under this rule (not results)

From the exact harness (`harness.py`, then `score2.py`), full table in
`tuning/exact_harness_on_tuning_set.md`. Its "incomplete" label is because the tuning set has 76
one-TOI stars and the criterion asks for 80.

| Metric | Confirmed planets (one TOI) | Injections (round-1 draws) | Criterion |
|---|---|---|---|
| false exclusion | 0 of 76 | 0 of 515 | <= 5% |
| median ln P excluded, SNR >= 8 | 0.72 (A1) | 0.67 | >= 0.5 |
| false event | 2 of 76 | 8 of 515 (own removed); 56 of 515 (all) | <= 5% |
| wrong period | 0 of 52 | 1 of 115 (own removed); 26 of 134 (all) | <= 3% |
| repeat recovery | 61 of 76 | 120 of 188 | >= 60% |
| probability calibration | observed 0, predicted 1.37, interval 0-4 | | inside the interval |

Multi-TOI stars (not judged): 2 of 33 false exclusions.

Known weak points going into the test:
- Repeat recovery in arm B is 64% against a 60% guard.
- In arm B only the SNR > 20 bin will reach the 60 trials needed to be judged (tuning: 26, 47 and 46
  trials in the bins 6-8, 8-12 and 12-20). The pre-registration's open item on this was not changed.
- The rule was chosen on 76 + 129 stars with few failures to learn from; the test can still fail.
