# Independent re-derivation of the round-3 test: where the spec was unclear, and what was chosen

Written before the primary rows and score were read (step 3). Sources used: `round3/FROZEN_RULE3.md`,
`round2/FROZEN_RULE.md`, `prereg_round3.yaml`, `prereg_round2.yaml`, `st2core.py`, `st2_pipeline.py`,
`gaia_dr4/candidates_frozen.csv`, the two test `*_defs.csv`, `round3/xcheck/independent/` (ind3.py,
its AMBIGUITIES.md and log) and `round2/test_independent/` (ind_test.py, AMBIGUITIES.md). Not read:
the round-3 harness, rule, scorer, drivers, tuning log, FREEZE3.json, samples, primary test outputs,
or the round-1/round-2 files on the exclusion list.

One accident: a first smoke test wrote to the session scratchpad, where an older `t_b.ckpt.jsonl`
left by another agent was picked up as a checkpoint and printed (four rows in the round-2
independent format: ref_snr, calls, no round-3 columns). It was not used; the smoke test was rerun
in an empty directory. Files in that scratchpad named `harness_rows.csv` and
`r3_xcheck_harness_rows.csv` were seen by name only and not opened.

## Per-trial script (`ind_test3.py`, a copy of ind3.py)

Departures of ind3.py from the frozen text, fixed in the copy:

1. **Choice of the alias event.** ind3.py chose "the highest SNR among clean events, else the
   highest SNR" (wording of a pre-freeze draft). The frozen text: "the passing event with the
   highest SNR". Under R4 every passing event is clean, so the result is the same; the dead branch
   and its `t2_boxclean` diagnostic are removed.
2. **`truth_excluded`, `n_clean_repeats`** were dropped by ind3.py; restored from the round-2
   independent script. `truth_excluded` uses round 2's rule exactly (w = 4, or 8 on a multi-TOI
   star; the reference's depth error; epochs t_cat + k P_true, k != 0, with min(T) <= epoch <= max(T)).

No other departure from the frozen text was found on reading ind3.py against it.

Ambiguities, most consequential first:

3. **`n_clean_repeats` (the denominator of repeat recovery).** "A second transit in clean data":
   true epochs other than the reference (further than dur/2 from t_cat) whose box of the reference
   `dur` in the merged detrended data has `frac >= 0.8` and is `>= 0.5 d` from a gap edge. In A2 the
   data still hold the transit, in B the injection. A trial counts with one or more such epochs.
4. **`dur` in `sec_box_snr`.** The trial's reference duration for the box length, the grid step and
   the exclusion zone around the event; the numerator is the event's fitted depth (over its own
   fitted duration). The frozen text now says this.
5. **Data for `sec_box_snr` in the reference sector.** "As the event search used it": the reference
   sector is detrended with the catalogued duration, every other sector with the reference `dur`.
6. **Support at an unconstrained epoch** (`h_k > 16`): cannot support ("neither excludes nor
   supports"), although the support test does not use the window.
7. **`frac_alias_excluded_r2`.** Same t_cat, same alias event (round 3's), same depth error,
   `h_k = w0`, no cap; denominator all aliases with P_n >= 1 d.
8. **No alias event:** call `no repeat`, `n_alias` 0, `event_true` False, fractions empty.
   **dT < 1 d:** no alias, call `all excluded`, `n_alias` 0, fractions empty (occurs: see item 20).
9. **Own events (arm B).** The null run uses the defs' catalogued values (with the 10% error);
   null events within one catalogued duration of t_cat are dropped; an injected-trial event is
   removed if it is not at a true transit and a null event lies within one (injected-trial) `dur`.
10. **A1 removal:** every epoch truth_t0 + k truth_P but the one nearest t_cat, +-1.5 catalogued
    durations, before detrending. **Arm-B masking:** +-0.75 d around the host's `t1_btjd`/`t2_btjd`
    in `gaia_dr4/candidates_frozen.csv`.
11. **Light-curve choice:** `pl.best_per_sector` (SPOC, TESS-SPOC, QLP), special sectors skipped,
    sectors under 50 points dropped. Whatever MAST lists today through the shared cache; if the
    primary run saw a different product list, rows will differ for that reason.
12. **pmin with no stellar mass or radius** (14 arm-A rows lack mstar, 2 lack rstar): NaN, nothing
    becomes `allowed_below_pmin`.
13. **Catalogue error (arm B).** Checked rather than applied. `depth_cat = inj_depth * max(0.2, 1 +
    0.1 g)` and the same for duration reproduce for all 1072 rows with numpy `default_rng([20261204,
    tic])`, two normals per trial, **only if the draws are consumed by the trials present, in row
    order**. For the four hosts with three trials (98233820, 352805517, 86433395, 301287939) a
    missing draw index does not consume random numbers; "two draws per trial in draw order" could
    also be read as indexed by draw number, which fails for 5 rows.

## Scorer (`ind_score.py`)

14. **Scored trial:** the harness finished (`truth_excluded` not empty). All 232 + 1072 did.
15. **"At least 40 called trials" for `wrong_period_rate_true_event`** ("same minimum of 40"):
    applied to that metric's own denominator (called trials whose alias event is true). The other
    reading (40 called trials of any kind) gives the same verdict here.
16. **False events in A1.** prereg: "share of A1 stars with at least one event that passes"
    (`n_vetted > 0`). FROZEN_RULE (round 2): a passing event "not at a true transit". The first is
    judged; the second is printed beside it.
17. **Repeat recovery in B:** the own-events-removed column is judged (removal never removes an
    event at the truth, so the two are equal by construction). A2: one-TOI stars.
18. **Alias power:** called = `resolved` or `narrowed`; A2 judged on one-TOI stars (the criterion
    says only "in A2"; all-star numbers beside it); B on the own-removed columns. Median of
    `frac_alias_excluded` against median of `frac_alias_excluded_r2`, same trials; met if
    r3 >= 0.8 r2.
19. **Outcome.** `incomplete` first (fewer than 60 one-TOI stars with A1 and A2 scored, or fewer
    than 600 B trials). `pass` if every criterion is met; `fail` if wrong_period_rate_true_event or
    alias_power is missed; `calls_after_vetting` if those two, repeat_recovery_rate and
    false_exclusion_rate are met and false_event_rate or wrong_period_rate is missed. **A miss of
    repeat_recovery_rate alone (or with the event criteria) has no outcome in the pre-registration**:
    the scorer says so rather than choosing one. `exclusion_regression` is a flag printed at the
    top, not an exclusive outcome.

## Things in the spec, defs or data that look wrong or worth a second look

20. **Test-set size against the pre-registration.** prereg: 124 arm-A stars (98 one-TOI), "all 124
    stars are run". Defs: 116 stars, 90 one-TOI, 26 multi (8 with period above 100 d, 7 of them
    one-TOI; prereg said 13 and 12). prereg B: "about 250 hosts"; defs: 269 hosts (69 / 100 / 100 by
    Tmag bin), 1072 trials, four hosts with three draws. Why 8 stars and 4 draws are missing is not
    in the files I may read.
21. `arm_a_defs.csv` has a column `at_grid_edge` (True for TIC 130633912 only) that no spec file
    defines.
22. Arm-B `t_cat` is up to 61 min from the true epoch (sd 14.9 min, a 4-sigma draw of the 15 min
    error); arm-A `t_cat` is offset by +3.0 min on average (sd 16 min, range -29 to +44 min) from
    truth_t0 + k truth_P, so the arm-A "catalogued time" also carries an error that the spec files
    do not describe.
23. T5 compares a depth fitted over the event's own duration (0.5 to 2 x `dur`) with the scatter of
    boxes of the reference duration: not a like-for-like SNR.
24. The cap `h_k > 16` leaves 12 steps of growth on one-TOI stars and 8 on multi-TOI stars.
25. `truth_in_alias_set` counts an `allowed_below_pmin` truth as in the set although such an alias
    is not counted in `n_alias_ok`; a trial can be "right" with its true alias outside the aliases
    the call is made from.
26. The catalogue-error seed is `[20261204, tic]`, the injection seed `[20261203, tic]`. The
    injection draws themselves (P, depth, b, t0, time error) were not re-derived: the draw order is
    not in the spec. Checked only that inj_depth = inj_k^2, inj_dur_d = the circular T14 of
    (P, a/R, b, k), truth = injected values.
27. `prereg_round3.yaml` still says "step 2 (the freeze) not yet approved"; FROZEN_RULE3.md says
    frozen and approved. `APPROVALS.md` was not opened.
28. The pre-registration's outcomes do not cover a miss of the repeat-recovery guard (item 19).
