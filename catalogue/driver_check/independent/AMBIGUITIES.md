# Independent duotransit alias table: choices made where the spec left one

Script: `ind_duo.py`. Sources read: `catalogue/prereg_catalogue.yaml` (step 4, hawthorn_duotransits),
`calibration/mono/round3/FROZEN_RULE3.md`, `calibration/mono/round2/FROZEN_RULE.md`,
`calibration/mono/monocal.py`, `st2core.py` lines 1-192 (not `alias_table`), `st2_pipeline.py`
`sector_events` / `passes`. Nothing on the forbidden list was opened.

1. **Which pipeline event gives `derr`.** The prereg says "passing event"; the task says "the pipeline
   event found within half a duration of `t2_cat`". I used *any* event in `prepare()["events"]`
   (top 3 per sector, no vetting cuts, not R4/T5), highest SNR if several. A passing-only reading
   would fall to the refit when the event near `t2_cat` fails the cuts. All 20 rows had an event.
2. **"Within half a duration"**: `|tc - t2_cat| <= dur/2` with `dur` the trial's reference duration
   (`dur_ref`), inclusive.
3. **Refit fallback** (not reached in these 20 rows): the first detrended sector (by sector number)
   with `t.min() <= t2_cat <= t.max()`, no 0.5 d margin; `core.refine_event(t, f, t2_cat, dur_ref)`;
   used if it returns a finite positive depth error, else the reference's `derr_ref`.
4. **`derr` enters every threshold.** `thr = 0.5 depth0 - 3 hypot(err, derr)` uses the second
   transit's `derr` for the box at the predicted time *and* for every offset box of the window;
   the reference's own depth error is not used anywhere in the table.
5. **Support and the cap.** An epoch with `h_k > 16` neither excludes nor supports (as written).
   Support does not require a measured reference and does not use the window.
6. **Epoch times are not rounded** to any dur/8 grid: the central box is at exactly `t1 + k P_n`,
   offset boxes at that time plus `j dur/8`, j = -h_k..h_k (j = 0 included).
7. **Condition (b)**: a window box vetoes the exclusion only if covered (`frac >= 0.5`) and
   `depth >= thr` with its own `err`/`thr`. Uncovered boxes (including empty ones, depth NaN) are
   ignored. Window boxes need not be clean.
8. **Edge distance** is measured to `core.gap_edges(T)` of all detrended data together (default
   0.5 d gap), not per sector. `sig = core.duration_scale_noise(T, F, dur)` on all data.
9. **`tmin`, `tmax`** = min and max of all detrended data `T`; k range as written, with k = 0 and
   k = n removed. `t1`, `t2` = the two catalogued times sorted.
10. **n range**: n = 1.. while `dT/n >= 1.0` (inclusive) and n <= 3000.
11. **Minimum period**: only an alias that would be `allowed` becomes `allowed_below_pmin`
    (round 3 wording); a `supported` alias below pmin stays `supported`. `pmin` uses the
    reference `dur` and `depth0`; NaN mass or radius gives no limit (3 rows).
12. **`sig_dt`** = hypot(15, 15) min = 21.213 min = 0.014731 d; `h_k = w0 + ceil(2 far_k (sig_dt/n)
    / (dur/8) - 1e-9)`; w0 = 4 (all 20 rows have `multi` False).
13. **`ref_measured`, `dur`, `depth0`** are taken from `mc.prepare` (its `out["ref_measured"]`,
    `dur_ref`, `depth_ref`), not recomputed.
14. **Clean** = `frac >= 0.8` and edge distance `>= 0.5 d` (both inclusive); the 3-sigma floor is
    `depth >= -3 err` (inclusive), upper bound `depth < thr` (strict).
