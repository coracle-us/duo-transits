# Round 3 test: independent rows and score against the primary

Made after `ind_a_rows.csv`, `ind_b_rows.csv`, `ind_score.md` and `AMBIGUITIES.md` were on disk
(hashes in `SHA256_before_compare.txt`, 2026-11-23T17:59:20Z). Primary files read only then:
`../test/arm_a_rows.csv`, `../test/arm_b_rows.csv`, `../test/score_summary.md`. Nothing in the
primary files was edited, and no independent file listed in the hash file was changed afterwards.
Script: `compare.py` (tables in `compare_columns.csv`, `compare_disagreements.csv`). Numbers agree
if within 1e-6 relative; booleans, strings and empty cells must match exactly.

## Result

**Every rule column agrees in every trial: 232 of 232 arm-A rows and 1072 of 1072 arm-B rows.**
The only differing cells are 8 values of one diagnostic (the reference sector's cadence), cause
traced below. The two criteria tables are the same number for number, and both scorers say the
outcome is not defined by the pre-registration.

## Per column

| Column (independent = primary unless named) | Arm A agree / disagree | Arm B agree / disagree |
|---|---|---|
| ref_snr | 232 / 0 | 1072 / 0 |
| ref_measured | 232 / 0 | 1072 / 0 |
| depth_ref | 232 / 0 | 1072 / 0 |
| dur_ref_d | 232 / 0 | 1072 / 0 |
| n_events | 232 / 0 | 1072 / 0 |
| n_clean_repeats | 232 / 0 | 1072 / 0 |
| truth_excluded | 232 / 0 | 1072 / 0 |
| n_vetted | 232 / 0 | 1072 / 0 |
| false_event | 232 / 0 | 1072 / 0 |
| repeat_detected | 232 / 0 | 1072 / 0 |
| t2 | 232 / 0 | 1072 / 0 |
| call | 232 / 0 | 1072 / 0 |
| n_alias_ok | 232 / 0 | 1072 / 0 |
| n_alias | 232 / 0 | 1072 / 0 |
| truth_in_alias_set | 232 / 0 | 1072 / 0 |
| event_true | 232 / 0 | 1072 / 0 |
| frac_alias_excluded | 232 / 0 | 1072 / 0 |
| frac_alias_excluded_r2 | 232 / 0 | 1072 / 0 |
| the same eleven call columns with `_own` | - | 1072 / 0 each |
| n_null_events | - | 1072 / 0 |
| n_r4 = primary r2_n_vetted (events passing R4 alone) | 232 / 0 | 1072 / 0 |
| false_event_r4 = primary r2_false_event | 232 / 0 | 1072 / 0 |
| repeat_detected_r4 = primary r2_repeat_detected | 232 / 0 | 1072 / 0 |
| the same three with `_own` | - | 1072 / 0 each |
| ref_author | 232 / 0 | 1072 / 0 |
| ref_cad_min = primary ref_cadence_min (diagnostic, not in the rule) | 231 / 1 | 1065 / 7 |

Primary `status` is `ok` for all 1304 rows; the independent run also had no failed row.

Not compared, because the independent rows do not hold them: `frac_lnP_excluded`, `truth_n_epochs`,
`ref_sector`, `n_sectors`, `n_points`, `baseline_d`, `null_event`, and the primary's round-2 call
columns `r2_t2`, `r2_call`, `r2_n_alias_ok`, `r2_truth_in_alias_set` (and `_own`). Those are round-2
quantities that no round-3 criterion uses.

## Every disagreeing cell

| File | TIC | Trial | Column | Independent | Primary |
|---|---|---|---|---|---|
| A | 268604331 | A1 | ref_cad_min | 30.000635 | 30.000566 |
| B | 311271011 | draw 0 | ref_cad_min | 10.000041 | 10.000059 |
| B | 311271011 | draw 1 | ref_cad_min | 10.000041 | 10.000059 |
| B | 282209355 | draw 1 | ref_cad_min | 10.000023 | 10.000041 |
| B | 13310825 | draw 0 | ref_cad_min | 10.000599 | 10.000614 |
| B | 409141608 | draw 1 | ref_cad_min | 9.999804 | 9.999887 |
| B | 409141608 | draw 2 | ref_cad_min | 9.999804 | 9.999887 |
| B | 409141608 | draw 3 | ref_cad_min | 9.999804 | 9.999887 |

Cause (traced for TIC 268604331 A1, by recomputing both): the independent script takes the median
time step of the reference sector before detrending (30.000635 min), the primary after detrending
(30.000566 min, reproduced exactly); detrending drops points and moves the median step by a few
parts in a million. The other seven cells have the same size and sign pattern and were not
recomputed one by one. The value is a label, used by no rule and no criterion.

## Criteria tables side by side

| Criterion, arm | Independent k / n | Primary k / n | Independent verdict | Primary verdict |
|---|---|---|---|---|
| wrong period, A2 one-TOI | 1 / 49 (2.0%) | 1 / 49 | met | pass |
| wrong period, B own removed | 3 / 208 (1.4%) | 3 / 208 | met | pass |
| wrong period (true event), A2 one-TOI | 0 / 48 | 0 / 48 | met | pass |
| wrong period (true event), B own removed | 0 / 205 | 0 / 205 | met | pass |
| false event, A1 one-TOI | 3 / 90 (3.3%) | 3 / 90 | met | pass |
| false event, B own removed | 10 / 1072 (0.9%) | 10 / 1072 | met | pass |
| repeat recovery, A2 one-TOI | 49 / 90 (54.4%) | 49 / 90 | MISSED | FAIL |
| repeat recovery, B own removed | 203 / 365 (55.6%) | 203 / 365 | MISSED | FAIL |
| alias power, A2 one-TOI (n = 48) | 0.9757 vs 0.9846, ratio 0.991 | 0.9757 vs 0.9846, 0.991 | met | pass |
| alias power, B own removed (n = 205) | 0.9838 vs 0.9848, ratio 0.999 | 0.9838 vs 0.9848, 0.999 | met | pass |
| false exclusion, A2 one-TOI | 0 / 90 | 0 / 90 | met | pass |
| false exclusion, B | 0 / 1072 | 0 / 1072 | met | pass |
| reported: wrong period, A2 multi-TOI | 6 / 14 | 6 / 14 | - | - |
| reported: wrong period (true event), A2 multi-TOI | 0 / 8 | 0 / 8 | - | - |
| reported: repeat recovery, A2 multi-TOI | 12 / 26 | 12 / 26 | - | - |
| reported: wrong period, B all events | 54 / 242 | 54 / 242 | - | - |
| reported: wrong period (true event), B all events | 0 / 188 | 0 / 188 | - | - |
| reported: repeat recovery, B all events | 203 / 365 | 203 / 365 | - | - |
| reported: alias power, B all events (n = 188) | 0.9839 vs 0.9848 | 0.9839 vs 0.9848 | - | - |
| reported: false event, B all events | 93 / 1072 | 93 / 1072 | - | - |

Scored: 90 one-TOI arm-A stars and 1072 arm-B trials on 269 hosts, in both.

Outcome, independent: "NOT DEFINED by the pre-registration" - missed = repeat_recovery_rate only.
Outcome, primary: "not defined by the pre-registration: wrong period (true event) and alias power
met, repeat recovery or false exclusion missed". Same conclusion. The primary sentence names
"repeat recovery or false exclusion" as a class; in this run it is repeat recovery alone (false
exclusion is 0 of 90 and 0 of 1072).

## Differences between the two score files that are not disagreements

- The primary reports alias power for A2 multi-TOI stars alone (n = 8, 0.9241 vs 0.9756); the
  independent scorer reports A2 all stars (n = 56, 0.9751 vs 0.9846). Neither is judged.
- The primary gives a one-sided 95% bound per rate; the independent scorer does not.
- The independent scorer also prints: A1 multi-TOI stars with a passing event (8 of 26), A2 one-TOI
  trials with a false event (4 of 90), arm-B null runs with an event passing round 3's vetting
  (100 of 1072), and false exclusions on A2 multi-TOI stars (0 of 26). The primary summary has no
  line for these.

## What this agreement does and does not show

- It shows that two implementations written from the same text give the same 1304 rows. Both rest
  on the same pipeline functions (`st2core`, `st2_pipeline`), the same light-curve cache, the same
  defs files and the same `gaia_dr4/candidates_frozen.csv`. An error in any of those is common to
  both and would not show here.
- `ind_test3.py` descends from `ind3.py`, which the pre-freeze harness check had already compared
  with the primary harness on 40 stars. It is independent of the primary code, not of that earlier
  reconciliation: a reading of the text that both sides settled on then is shared now.
- The defs were taken as given. Checked: the arm-B catalogue error (reproduces for 1072 of 1072
  rows under one reading of the draw order, AMBIGUITIES 13) and the internal consistency of the
  injected depth and duration. Not checked: the selection of the arm-A stars and arm-B hosts, the
  injection draws, the arm-A truth periods and epochs.
