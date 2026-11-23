# Round 2 test: result

**Outcome: partial.** Period exclusion meets every criterion. Automatic second-transit calls do not.
Under the pre-registration (`../prereg_round2.yaml`, 1.1) this means: exclusions can be published
with their probabilities; second-transit calls need vetting by a person.

Rule: `FROZEN_RULE.md`, frozen 2026-11-23 (`../APPROVALS.md`, `FREEZE.json`) before any test-set
light curve was loaded. Run: `run_test.py` at commit 1228a24. Score: `score2.py`, written before the
freeze. Numbers below are from `test/score_summary.md`; rows are in `test/arm_a_rows.csv` and
`test/arm_b_rows.csv`.

## Sample

- Arm A: 146 confirmed-planet stars drawn, all scored; 123 have one TOI (the primary analysis).
- Arm B: 600 injections on 150 hosts. 200 hosts were tried in the frozen order; 50 have no SPOC,
  TESS-SPOC or QLP light curve at MAST (42 from catalogue `t16`, 8 from `k2`; four checked by hand,
  all with zero products) and were replaced by the next in their bin, as pre-registered.

## Criteria

| Criterion | Confirmed planets, one-TOI stars | Injections | Limit | Verdict |
|---|---|---|---|---|
| False exclusion, overall | 0 of 123 | 0 of 600 | <= 5% | pass |
| False exclusion, SNR bins with >= 60 trials | none qualifies | 0 of 71 (12-20), 0 of 254 (>20) | <= 5% | pass |
| Median ln P excluded, SNR >= 8 | 0.65 (71 stars) | 0.67 (377 trials) | >= 0.5 | pass |
| Probability calibration | observed 0, predicted 1.63, interval 0-5 (both arms, 723 trials) | | inside | pass |
| False event | 7 of 123 (5.7%) | 11 of 600 (1.8%), own events removed | <= 5% | **miss** (arm A) |
| Wrong period | 6 of 76 (7.9%) | 6 of 151 (4.0%), own events removed | <= 3% | **miss** (both) |
| Repeat recovery | 81 of 123 (66%) | 150 of 230 (65%) | >= 60% | pass |

Reported, not judged: false exclusion on multi-TOI stars 0 of 23; arm B with the hosts' own events
kept: false event 84 of 600, wrong period 46 of 178; 85 of 600 null runs have a passing event
(51 of 150 hosts have at least one).

## What the exclusion result does and does not show

- No true period was excluded in 723 trials. The one-sided 95% upper bounds are 2.4% (arm A) and
  0.5% (arm B).
- The SNR 6-8 and 8-12 bins hold 17 and 15 arm-A stars and 46 and 52 injections: too few to be
  judged, as expected before the run. The low-SNR claim rests on the overall rate only.
- A reference the pipeline could not measure (SNR below 5) excludes nothing, by the rule.
- The arm-A sample has periods of 5 to 100 d; most candidates in the catalogue are longer.

## Review of every wrong call and false event (done after the outcome; explanations are hypotheses)

False exclusions: none to review.

**Wrong period, second transit real (6 trials).** TOI-2127.01, 3919.01, 5081.01, 1406.01, 1173.01
and injection TIC 141522677 draw 3. In each the event found is a true transit (within 0.2 durations
of the true ephemeris) but the alias table excluded the true alias. In all six the period from the
two transit times, P = dt/n, is off by enough that the predicted transit has moved by 1.2 to 5.1
durations at the far end of the data. For calls that were right, that drift has a median of 0.09
(arm A) and 0.05 (arm B) durations, and a 90th percentile of 0.55 and 0.31. So the alias table
tests a period more exactly than two transit times fix it. The single-transit exclusion has a
timing tolerance (R3); the alias table has none that grows with the number of orbits.

**Wrong period, event not a transit (6 trials).** TOI-569.01 (arm A2), and injections TIC 76804724
draws 0 and 3, 82680046 draw 1 (event 1.2 durations before a true transit), 265465927 draw 1,
271693283 draw 3. These follow from the false events below.

**False events, arm A1 (7 stars).** TOI-569.01, 6038.01, 1255.01, 1432.01, 2567.01, 1148.01,
1648.01. Each event is 11 to 68 durations from the planet's ephemeris, so it is not a leaked
transit. Four of the seven have reference SNR below 7. Causes not looked into star by star.
With six such stars the criterion would have been met (4.9%); it was missed by one star.

**False events, arm B, own events removed (11 trials on 9 hosts).** Five of them give a wrong
period; five leave the truth in the alias set; one (TIC 265465927 draw 2) ends with every alias
excluded, which is not a call.

## Limits

- Independent re-derivation (`test_independent/`, a fresh agent that did not read the harness, the
  scorer or the rows until its own were on disk): 292 of 292 arm-A and 600 of 600 arm-B trials
  agree in every compared column, and all criteria rows and the outcome agree. Its predicted
  false-exclusion count is 1.6272 against 1.6282 here (unrounded against four-decimal
  probabilities). Both use the same pipeline functions, cached light curves and trial definitions,
  so a fault shared by both would not show. Its 24 spec ambiguities are in `AMBIGUITIES.md` there.
- **Vetting by a person does not catch the six real-transit wrong periods**, because the event is
  genuine; the fault is in the alias table. Four of the five arm-A cases are called resolved at
  n = 1 (for example 14.13 d for a 7.06 d planet). So "second-transit calls need vetting" is not
  enough on its own: alias tables from two transit times should not be published as they are for
  stars where the period error drifts by more than about half a duration over the data.
- Arm B, TIC 82680046 draw 1: the injected transit was recovered 11.8 h from its true time (half a
  duration is 10.1 h), so it counts as a false event and a wrong period; its true alias is labelled
  below the minimum period from the duration, which the call rule treats as a hard limit.
- The arm-A periods drawn reach only 38.5 d. Arm B's catalogued depth and duration are exact; only
  the time has the 15 min scatter.
- The probability-calibration criterion is weak here: with 1.63 predicted, any count from 0 to 5
  passes.
- The explanations above were found after the outcome. A repair (for example a period tolerance in
  the alias table) would be a new rule and needs a new pre-registered test on stars not yet used.
- Arm B hosts carry their own signals: 14% of null runs have a passing event. The primary analysis
  removes them, as approved in 1.1; the all-events numbers are above.
- The pre-registration says the arm-A pool has 246 stars; the stated rule gives 257. The draw of
  146 and the one-TOI count match.
