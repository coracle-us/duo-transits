# tess_wide: wide blind search of TESS FFI light curves for long-period single transits

Status 2026-11-25: **pilot on S104 done; stopped after the pilot report, as asked.** Nothing was submitted,
emailed or posted. Files that name candidate stars are in `tess_wide/data/` (a link to
`duo-transits-cache/tess_wide`, outside git) and are not public.

- `PLAN.md`: the plan, v1.0, committed (3f3da66) before any search data was downloaded. Every threshold has
  its reason there.
- `PILOT_S104.md`: the pilot report (counts only), written by `s7_report.py`.
- Scripts, in order:
  - `s1_starcut.py`: TIC rows and the star cut;
  - `s2_download.py`: QLP download, keeping compressed arrays;
  - `s3_search.py`: `search_any` (imported) on the real flux, 2 − flux and injections;
  - `s4_stages.py`: common-mode veto, one event per star, size test;
  - `s5_cnn.py`: the Salinas CNN, run in its own environment from its local clone and not copied;
  - `s6_crossmatch.py`;
  - `s7_report.py`.
- `inject.py` is the batman injection model shared by s3 and s5.

## Result in one paragraph

- **Searched:** 65,881 stars passed the cut (T ≤ 12, R★ ≤ 2 R☉, 3000–7500 K), out of S104's 1,152,719 QLP light
  curves.
- **Funnel:** 11,305 detections → 8,618 after the common-mode veto → 1,532 events on stars with exactly one event →
  1,412 after the size test → 225 pass the CNN → 203 not in any catalogue checked (22 known).
- **The control fails the rule fixed in PLAN.md:** the identical search on 2 − flux ends with 300 CNN passes
  against the real search's 225, above the 0.5× limit. Even before the CNN, the control has more single events
  than the real search (1,898 against 1,532). By that rule, **the 203 are not candidates**: in this sector the
  pipeline does not separate transits from systematics, and no threshold was changed afterwards.
- **Completeness:** at detection it levels off at 55–70 % for depths above ~3,000 ppm and falls below 10 % under
  ~650 ppm. End to end (through the CNN) it is 35–50 % above 3,000 ppm and 23 % ± 1 % over the whole injected
  range.

## What the pilot shows (read after the run; no rule was changed)

- **The CNN adds almost no discrimination.** It passes 16 % of the real single events, 17 % of the control ones
  and 67 % of injected transits. It keeps real transits, but it also measures a transit-like shape in most
  one-off systematics, which is what a re-measurement network (not a classifier) should be expected to do.
- **Systematics dominate the survivors.** Of the CNN passes:
  - 181 of 225 real and 212 of 300 control carry the PLAN 5.4 background or centroid flag;
  - against 224 of 699 injected transits (32 %).
  - Camera 1 holds the most passes, real and control.
  - With the flags applied as a cut, there would still be 88 control passes against 44 real ones. A flag cut
    alone would not rescue this sector.
- **SAP flux was the wrong trade for a blind search.** Taking SAP avoided distortion by the QLP spline, at the cost
  of scattered-light structure that both signs of the search pick up. The 0.5× control rule caught this, which is
  what it was for.

## If the search continues (proposals only; none run)

A v1.1 would need Jon's decision, a new PLAN version committed before the run, and a test on a **different
sector**, so that nothing is tuned on S104. Ideas, roughly in order of expected gain:
1. Search QLP's detrended `DET_FLUX` (or SAP after a scattered-light regression on `SAP_BKG`) rather than raw SAP.
   Measure the loss on long transits with the same injections.
2. Make the background and centroid tests cuts, with thresholds set on injections: 32 % of injected transits are
   flagged at 3σ now, so the background statistic needs a red-noise term first.
3. Replace the CNN step with a classifier trained on transits against systematics (for example inverted-flux
   events as negatives), or drop it as a vetter.
4. Per-camera control rule: camera 1 carries the most passes.

## Cost of the full search (from the pilot)

All 98 QLP sectors with bulk scripts on MAST (S74–S79 have none) come to about 5.9 million star-sectors:
- download: about 2.5 TB, about 4 days at the pilot's 18 stars/s, with 0.6 TB kept;
- search: about 2 h on 16 cores;
- the CNN takes minutes.

Free disk on 2026-11-25 was 2.7 TiB, enough for the kept arrays but not to keep the FITS files.
