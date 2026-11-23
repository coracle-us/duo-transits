# Round-3 independent per-trial scoring: where the spec was unclear, and what was chosen

Written before any comparison with the primary harness. Sources used: `round3/FROZEN_RULE3.md`,
`round2/FROZEN_RULE.md`, `prereg_round3.yaml`, `prereg_round2.yaml`, `st2core.py`,
`st2_pipeline.py`, `gaia_dr4/candidates_frozen.csv`, `round3/xcheck/defs.csv`, and
`round2/test_independent/ind_test.py` with its `AMBIGUITIES.md` (the starting point). Nothing else
in `round3/`, no primary harness, rule, scorer, driver or tuning log was read.

Run: 120 definitions (40 stars), 120 trials ran, 0 failed (`ind_rows.csv`, log `ind_log.txt`).

## Round 3, most consequential first

1. **Which `dur` in `sec_box_snr` (T5).** The text says boxes "of length `dur`" and defers symbols to
   round 2, where `dur` is the trial's reference duration. The event has its own fitted duration
   (0.5 to 2 x `dur`), and its `depth` is the mean over that fitted duration. Chosen: the reference
   `dur` for the box length, the grid step (`dur/2`), the `> dur` exclusion around the event, and
   the fitted event `depth` in the numerator. So numerator and denominator can be boxes of different
   length.
2. **"The detrended data of the event's own sector".** Chosen: the (t, f) the event search ran on.
   For the reference sector that is the copy detrended with the catalogued duration; every other
   sector is detrended with the reference `dur` (as the pipeline and ind_test.py). Injected trial:
   the injected light curve. Other transits in the sector (arm A2, the reference transit) are not
   masked; only boxes within `dur` of the event are left out.
3. **"Clean" in the choice of alias event.** "The highest SNR among clean events, else the highest
   SNR". R4 already requires every passing event to be clean (`frac >= 0.8`, `edge_dist_d >= 0.5`
   of the event record), so with that meaning the sentence is the same as round 2's "highest SNR".
   Chosen: that meaning. The other reading (the box of length `dur` at the event time in the merged
   data, with the merged data's gap edges) is written as the diagnostic `t2_boxclean`; it picks the
   same event in all 120 trials (and in the `_own` versions).
4. **Support at an unconstrained epoch.** "It neither excludes nor supports" is taken literally: an
   epoch with `h_k > 16` cannot support, although the support test never uses the window.
5. **Depth error in the table.** As round 2 (ind_test.py items 1-2): the alias event's own
   `depth_err` replaces `depth0_err` in `thr` (centre box and every offset box) and in the support
   test. Same in the round-2 comparison table.
6. **Exclusion clause (b) at j = 0.** The centre box is already required to be below `thr` by (a),
   so j = 0 adds nothing. Offsets are `c + j dur/8`, j = -h_k .. h_k; a box counts only if covered
   (`frac >= 0.5`).
7. **`frac_alias_excluded_r2`.** Round 2's table with the same t_cat, the same alias event (the
   round-3 one, not the event round 2's vetting would pick), the same depth error, `h_k = w0`, no
   cap. Denominator: all aliases in the table (the two tables have the same aliases).
8. **Values with no alias event (`no repeat`).** `t2` empty, `n_alias_ok` = 0, `n_alias` = 0,
   `truth_in_alias_set` = False, `event_true` = False, both fractions empty.
9. **dT < 1 d (empty table).** Call `all excluded`, `n_alias` = 0, fractions empty. Did not occur.
10. **Units of `sig_dt`.** 21.2132 min = hypot(15, 15)/1440 d, unrounded.
    `h_k = w0 + ceil(2 far_k (sig_dt/n) / (dur/8) - 1e-9)` with `dur` the reference duration.
11. **`allowed_below_pmin`.** Only an `allowed` alias changes; a `supported` alias below the minimum
    period stays `supported` and counts in `n_alias_ok`. `pmin` not finite (2 defs rows lack mstar):
    none.
12. **Own events (arm B).** As ind_test.py item 12: the null run is searched with the draw's
    catalogued values from the defs (which already carry the 10% error); its events within one
    catalogued duration of t_cat are dropped; an injected-trial event is removed if it is not at a
    true transit and a null event lies within one `dur` (the injected trial's reference duration)
    of it. Removal is before vetting, so `n_vetted_own` counts what is left.
13. **`_own` columns.** The task text says "nine call columns" and names a range of eleven
    (`n_vetted` ... `frac_alias_excluded_r2`). All eleven are written with `_own`; they are empty
    for arm A.
14. **Catalogue error (arm B).** Not applied here; `depth_cat` and `dur_cat_d` are used as given.
    The `max(0.2, 1 + 0.1 g)` rule and its seed were therefore not re-derived or checked.

## Inherited from round 2 (unchanged choices; see round2/test_independent/AMBIGUITIES.md)

Light-curve choice (`pl.best_per_sector`, special sectors skipped), minimum 50 points per sector,
the reference sector (first whose range +-0.5 d holds t_cat), A1 removal (+-1.5 catalogued
durations around every true epoch but the one nearest t_cat), arm-B masking (+-0.75 d around the
host's catalogued times) and injection (batman, exposure-integrated), dropping events within one
`dur` of t_cat in the reference sector only, "at a true transit" (strictly within `dur/2`, true
epochs further than `dur/2` from t_cat), `sig` from the merged data.

## Things in the spec that look wrong or worth a second look

- FROZEN_RULE3 "The alias event is chosen ... the highest SNR among clean events, else the highest
  SNR" has no effect under R4 (item 3). Either it is a leftover from the pipeline's
  `alias_events`, or "clean" means something else and should be defined.
- prereg_round3 `harness_check` asks for agreement on "every status in every alias table" on
  "20 tuning-set stars". The defs hold 40 stars, and the requested columns carry counts only. The
  diagnostics `n_alias_excluded`, `n_alias_supported`, `n_alias_below_pmin` and `alias_ok_n` (the n
  of every allowed or supported alias, cut at 400 characters) are written to allow a closer check.
- T5 compares a depth fitted over the event's own duration with the scatter of boxes of the
  reference duration (item 1). For an event fitted at half or twice `dur` the ratio is not a
  like-for-like SNR.
- On a multi-candidate star (`w0` = 8) the window can grow by only 8 steps before the cap of 16; on
  the others by 12. The text says the cap is "+-2 durations" and does not say this is intended.
- With an unmeasured reference nothing is excluded, so a passing event gives `narrowed` with nearly
  every alias allowed (TIC 98237607 draw 0: 1285 of 1434, fraction excluded 0). This follows
  round 2's R2 but is not said in the "call" section.
- The seed for the catalogue error is `[20261204, tic]`; the injection seed in prereg_round3 is
  `[20261203, tic]`. Probably intended (two generators); noted because they differ by one.
- FROZEN_RULE3 says "Status: proposed" and refers to `harness3.FROZEN3`; prereg_round3 still says
  step 2 is not approved. Consistent with a check before the freeze.

## What T5 and T1 did on these 120 trials (description, not a comparison)

- T5 removed at least one R4-passing event in 12 trials and turned 2 trials into `no repeat`.
- Calls: A1 15 no repeat, 3 narrowed, 2 all excluded; A2 17 narrowed, 1 resolved, 2 no repeat;
  B 22 narrowed, 58 no repeat (own events removed: 20 narrowed, 60 no repeat).
