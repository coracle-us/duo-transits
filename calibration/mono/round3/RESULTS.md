# Round 3 test: result

**Outcome: not defined by the pre-registration.** Every criterion on the correctness of calls was
met in both arms. The guard on repeat recovery (at least 60%) was missed in both arms. The
pre-registration (`../prereg_round3.yaml`, 1.0) names no outcome for that combination; this was
said to Jon before the freeze. It is not "pass" and, by the letter, not "calls after vetting"
either, because that outcome also requires the recovery guard.

Rule: `FROZEN_RULE3.md`, frozen 2026-11-23 (`../APPROVALS.md`, `FREEZE3.json`) before any round-3
test-set light curve was loaded. Run: `run_test3.py` at commit 97ffd3c. Score: `score3.py`, written
before the freeze. Numbers are from `test/score_summary.md`; rows in `test/arm_a_rows.csv` and
`test/arm_b_rows.csv`.

## Sample

- Arm A: 124 confirmed-planet stars in the sample; 8 have under 3 clean transits and are skipped
  (the rule of rounds 1 and 2); 116 scored, 90 with one TOI (the primary analysis).
- Arm B: 1,072 injections on 269 hosts (69 at Tmag below 10, 100 at 10-12, 100 above 12). 356
  hosts were tried in the fixed order; 87 have no SPOC, TESS-SPOC or QLP light curve.

## Criteria

| Criterion | Confirmed planets, one-TOI stars | Injections, own events removed | Limit | Verdict |
|---|---|---|---|---|
| Wrong period | 1 of 49 (2.0%) | 3 of 208 (1.4%) | <= 3% | pass |
| Wrong period, real second transit | 0 of 48 | 0 of 205 | <= 3% | pass |
| False event | 3 of 90 (3.3%) | 10 of 1072 (0.9%) | <= 5% | pass |
| Alias power against round 2's table | 0.991 (0.9757 vs 0.9846) | 0.999 (0.9838 vs 0.9848) | >= 0.8 | pass |
| False exclusion (replication of round 2) | 0 of 90 | 0 of 1072 | <= 5% | pass |
| Repeat recovery | 49 of 90 (54%) | 203 of 365 (56%) | >= 60% | **miss** |

One-sided 95% upper bounds: wrong period 9.3% (arm A) and 3.7% (arm B); with a real second
transit 6.1% and 1.5%. Arm A is small, as said before the run.

Reported, not judged: multi-TOI stars - wrong period 6 of 14, recovery 12 of 26. Arm B with the
hosts' own events kept - wrong period 54 of 242, false events 93 of 1072.

## What the result shows

- **The repair of round 2's fault held on new stars.** With a real second transit the true period
  stayed in the alias set in 48 of 48 and 205 of 205 called trials. In round 2 it was lost in 5 of
  75 and 1 of 146. The table kept 99% of round 2's excluding power.
- **The single-transit exclusion rule replicated**: no true period excluded in 1,162 further trials
  (0 of 723 in round 2).
- **Calls are right when made, but 45% of clean repeats are not found.**

## Where recovery was lost (looked at after the outcome)

| Reference SNR | Arm A recovered | Arm B recovered |
|---|---|---|
| below 6 | 2 of 34 | 4 of 117 |
| 6 to 8 | 4 of 12 | 5 of 24 |
| 8 to 12 | 9 of 10 | 17 of 33 |
| 12 to 20 | 8 of 8 | 41 of 47 |
| above 20 | 26 of 26 | 136 of 144 |

- Below a reference SNR of 6 almost nothing is recovered; above 12, 211 of 225.
- Round 2's vetting on these same trials recovers 53 of 90 and 214 of 365, both 59% - also under
  the guard. The round-3 veto cost 4 and 11 recoveries. So the miss is mostly the sample (38% of
  arm-A stars and 32% of arm-B trials with a clean repeat have reference SNR below 6), not the
  repair. This comparison was not pre-registered.
- The guard was written to stop vetting from removing real repeats. It does not separate that from
  a sample with more shallow transits.

## Review of every wrong call and false event (hypotheses, made after the outcome)

**Wrong period, arm A (1).** TOI-2459.01 (TIC 192790476, P = 19.1 d): the alias event is not a
transit of the planet; the same star gives the false event in A1.

**Wrong period, arm B (3).** TIC 339399841 draw 1, 143355733 draw 1, 139695717 draw 1. All three
alias events are not at an injected transit and were not seen in the null run; reference SNR 5.6,
8.1 and 6.5; periods 751, 744 and 158 d.

**False events, arm A1 (3 stars).** TOI-1232.01 (34 passing events on TIC 364395234, so something
periodic that is not the TOI), TOI-1443.01, TOI-2459.01.

**False events, arm B (10 trials).** Three give the wrong periods above; five leave the truth in
the alias set; two end with every alias excluded.

False exclusions: none.

## Limits

- Independent re-derivation (`test_independent/`, a fresh agent; its rows and score were hashed
  before it read any primary output): every rule column agrees in all 232 arm-A and 1,072 arm-B
  trials, and every criterion line has the same k and n. It reaches the same "not defined" outcome.
  Its only difference is a diagnostic (reference cadence, parts per million, 8 cells). It shares the
  pipeline functions, the cached light curves and the trial definitions with the primary run, and
  its script descends from the pre-freeze cross-check, so a fault common to both would not show.
- Its recovery analysis (`test_independent/recovery_analysis.md`, after the outcome) agrees with
  the one above: without the round-3 veto recovery is 53 of 90 and 214 of 365. Of the missed
  repeats, 35 (arm A) and 95 (arm B) fail round 2's vetting, mostly event SNR below 8; 4 and 11
  are lost only to the new veto; 2 and 56 are not among a sector's top three events. Trials with
  an unmeasured reference recover 0 of 26 and 3 of 103. QLP-only repeats recover 31% against 71%
  for SPOC, confounded with faint hosts.
- Points it raised: `truth_in_alias_set` counts a true alias labelled below the minimum period as
  in the set although that alias is not counted in `n_alias_ok`; four arm-B hosts have 3 draws,
  not 4 (no clean transit found for one draw); `prereg_round3.yaml` still reads "step 2 not yet
  approved" because the file is frozen by hash - the freeze is recorded in `../APPROVALS.md`.
- Arm A's periods: median 7.8 d, quartiles 6.0 and 11.5 d; only the 13 stars added above 100 d
  resemble the catalogue's periods, and some of those were among the 8 skipped.
- Few calls are "resolved" (5 in arm A, 17 in arm B); most are "narrowed", with a median of 15
  and 11 aliases left.
- Arm B hosts carry their own signals; with them kept, 22% of called trials have a wrong period.
  On a real candidate there is no null run to remove them, so a second event still needs a person.
- No unused confirmed planet with period above 5 d is left for a further round.
