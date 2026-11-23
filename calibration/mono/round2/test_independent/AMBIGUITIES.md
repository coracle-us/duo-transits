# Independent re-derivation of the round-2 test: where the spec was ambiguous, and what was chosen

Written before the primary rows and score were read (step 3). Sources used: `prereg_round2.yaml`
(1.1), `FROZEN_RULE.md`, `st2core.py`, `st2_pipeline.py`, `gaia_dr4/candidates_frozen.csv`,
`xcheck/independent/ind2.py`, the two test `*_defs.csv`. Also read, because FROZEN_RULE.md defers to
it ("everything not listed here is as in round 1"): `calibration/mono/prereg_mono_calibration.yaml`
(round 1's pre-registration) and `round2/FREEZE.json` (hashes and seeds only). Not read: the
harness, rule, scorer, tuning log, round-1 drivers, XCHECK.md, any primary test output.

## Per-trial harness (`ind_test.py`)

Departure of ind2.py from FROZEN_RULE.md, fixed in the copy:

1. **Depth error in the alias table.** FROZEN_RULE "Calls": "In the alias table only, `depth0_err`
   is the alias event's own depth error, not the reference's". ind2.py used the reference's error
   for `call` / `truth_in_alias_set` and kept the event's error as a side variant (`_everr`).
   ind_test.py uses the event's error for `call`, `n_alias_ok`, `truth_in_alias_set` and their
   `_own` versions; the reference-error table is kept as `*_referr*`.

Ambiguities, with the choice made (all inherited from ind2.py unless marked):

2. **Which depth error in the alias table's *exclusion* test.** The sentence in 1 could mean only
   the `supported` test or both. Chosen: both the exclusion threshold `thr` (R1-R3, including the
   R3 offsets) and the `supported` test use the event's error inside the alias table.
3. **Light-curve choice.** Round 1's pre-registration says "QLP where available, else TESS-SPOC,
   else own TESScut"; the pipeline's `best_per_sector` ranks SPOC, then TESS-SPOC, then QLP.
   Chosen: the pipeline's ranking (`pl.best_per_sector`), sectors in `pl.SPECIAL_SECTORS` skipped.
   The two documents contradict each other.
4. **Minimum sector size.** Not in the spec. Chosen: a sector with fewer than 50 points after
   loading (or after masking / removal) is dropped.
5. **Reference sector.** Chosen: the first sector (in sector order) whose time range, widened by
   0.5 d, holds t_cat (as `pl.search_one`). It is detrended with the *catalogued* duration and that
   detrended copy is the one used everywhere afterwards (as the pipeline, which caches it); all
   other sectors are detrended with the reference duration `dur`. If t_cat falls in no sector the
   reference is unmeasured and ref_snr is empty.
6. **A1 removal.** "Other transits removed" is not defined in round 2. From round 1's
   pre-registration: +-1.5 durations. Chosen: every epoch truth_t0 + k truth_P except the one
   nearest t_cat, +-1.5 catalogued durations, removed before detrending.
7. **Arm B masking.** +-0.75 d around every `t1_btjd` and `t2_btjd` of the host in
   `gaia_dr4/candidates_frozen.csv`, before injection. The injected model is batman, quadratic limb
   darkening (u1, u2 from the defs), circular, inc from inj_b and inj_ar, multiplied into the
   normalised (10-min-binned where the pipeline bins) flux, exposure-integrated at the sector
   cadence (supersampled to about 2 min).
8. **R3 offsets.** j = -w..w; j = 0 is the box at c itself and is already tested by clause 2.
9. **"Inside the data span"** for the true period's epochs: min(T) <= epoch <= max(T).
10. **ln P grid.** P_i = exp(i * step), i = 0 .. floor(ln 1000 / step); the fraction is the share of
    grid periods excluded. Time grid: min(T) + i dur/8 up to the first point at or past max(T);
    an epoch outside the grid excludes nothing. Unmeasured reference: 0.
11. **Dropping events near the reference.** "Within one `dur`": |tc - t_cat| < dur, with the
    reference (measured or catalogued) `dur`, in the reference sector only.
12. **Null run.** The null run goes through the same search (so its own reference rule and its own
    `dur`); events within one *catalogued* duration of t_cat are dropped in any sector. An injected-
    trial event is the host's own if a null event lies within one (injected-trial) `dur` of it and
    it is not at a true transit.
13. **"At a true transit".** Epochs truth_t0 + k truth_P with |epoch - t_cat| > dur/2; match within
    dur/2 (strict), `dur` = the trial's reference duration.
14. **pmin with missing stellar mass or radius** (8 arm-A rows lack mstar, 2 lack rstar):
    `core.pmin_from_duration` returns NaN and no alias is `allowed_below_pmin`. Round 1's "solar
    values" note is about injections; not applied here.
15. **dT < 1 d** (no alias at all): call `all excluded`.
16. **`sig`.** `core.duration_scale_noise(T, F, dur)` on the merged data: its cadence is the median
    over all sectors, as the function is written.

## Scorer (`ind_score.py`)

17. **Scored trial.** A row whose harness run finished (truth_excluded not empty).
18. **SNR bins.** [6,8), [8,12), [12,20), [20,inf) of `ref_snr` (the `refine_event` SNR at the
    catalogued time), whether or not the reference counts as measured. Judged per arm (A2 one-TOI
    and B separately), not pooled; a bin is judged only with >= 60 trials in that arm.
19. **ln P criterion.** Median of frac_lnP_excluded over trials with ref_snr >= 8, in A1 one-TOI
    stars and in all B trials, each >= 0.5. The pre-registration does not say whether "A1" here is
    one-TOI only; one-TOI taken because the primary analysis is one-TOI stars.
20. **Repeat recovery.** Denominator: trials with n_clean_repeats >= 1. A2 is taken as one-TOI
    stars (the criterion says only "in A2"). In B the own-events-removed version is the judged one.
21. **Wrong period.** Denominator: calls `resolved` or `narrowed`; `all excluded` is not counted.
22. **Probability calibration.** A trial's R5 bin is "unmeasured" if the reference is not measured
    (whatever its SNR), else by ref_snr: [5,8), [8,12), [12,20), [20,inf). p = (k+0.5)/(n+1) from the
    table (unrounded). Predicted = sum of p over scored A2 one-TOI and B trials, unmeasured ones
    included (as FROZEN_RULE's "Calibration test"). Interval: Poisson quantiles 0.025 and 0.975
    of the predicted mean, inclusive. Only the total is judged ("per bin and in total" in the metric,
    "the total" in the criterion).
23. **Outcome.** `incomplete` is tested first (one-TOI arm-A stars with both A1 and A2 scored < 80,
    or arm-B scored trials < 400). The pre-registration gives no outcome when false exclusion and
    ln P are met but probability calibration is missed (not pass, not partial, not fail); the
    scorer says so if it happens.
24. **Null-run rate** (reported, not judged): share of B trials whose null run has an event that
    passes the R4 cuts.
