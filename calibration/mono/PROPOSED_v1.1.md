# Proposed version 1.1 of the single-transit pre-registration

**Status: proposal, not approved.** Version 1.0 (`prereg_mono_calibration.yaml`) stays frozen. Written
2026-11-22, after approval of 1.0 and before any trial outcome was produced or looked at. The only
real-data run so far is the timing pilot, which discards outcomes in memory.

Building the harness showed that 1.0 leaves six things undefined or wrong. Each needs a decision
before the run, so each is a version change that needs Jon's approval.

| # | In 1.0 | Problem | Proposed for 1.1 |
|---|---|---|---|
| 1 | "the excluded set" of a single transit, from `st2core.py` | `st2core` has no such function. It builds alias tables only once a second event is found. | Define it in `monocal.py` from `core.epoch_verdicts`: a period P is excluded if any epoch t_cat + kP (k ≠ 0) is a clean epoch with verdict `excludes`. The true period is tested at its exact value. For the fraction of ln P excluded, periods are log-spaced from 1 to 1000 d with a step that moves the farthest epoch by a quarter duration, and epochs are looked up on a time grid of one eighth of a duration. |
| 2 | Pass criterion on false exclusion applies to "arms A1 and B" | In A1 every other transit is cut out, so the true period's epochs have no data and cannot be excluded. The criterion is passed by construction. | False exclusion is judged in **A2 and B**, where real or injected transits sit at the true period's epochs. A1 reports only false repeats, the wrong-period rate and the fraction of ln P excluded. |
| 3 | "wrong_period_rate among trials the pipeline calls resolved or narrowed" | "resolved" and "narrowed" are not defined. | A passing event exists and the alias table through t_cat and that event has exactly one alias `allowed` or `supported` (resolved) or more than one (narrowed). The truth is in the set if \|Δt − nP_true\| < half a duration for some n and alias n is not `excluded`. |
| 4 | Catalogued events masked by ±1.5 durations | `candidates_frozen.csv` carries no durations. | Mask ±0.75 d around each catalogued time (1.5 × a 12 h transit). |
| 5 | Limb darkening "from the TIC Teff and logg" | Needs a coefficient table that is not in the repo. | Either (a) fixed quadratic coefficients (0.35, 0.25) for every star, or (b) Claret (2017) TESS coefficients from VizieR (J/A+A/600/A30), which is a download Jon would have to allow. Recommended: (a). The exclusion rule uses a box mean, and the coefficients change that by a few percent. |
| 6 | (not stated) | The catalogued time of the single transit has an error; the harness has to choose one. | Gaussian scatter of 15 min on the catalogued time, as in the 1 November test. Stars with no TIC mass or radius use solar values and are flagged. The reference sector is also searched for repeats, except within one duration of the reference. |

Also to record in 1.1:
- `pipeline_commit`: `st2_pipeline.search_sector` was split into `search_sector` + `sector_events` so
  the flux can be changed before detrending. No threshold changed; `tests/test_synthetic.py` passes.
- Runtime (the first open item): see `results/pilot/arm_b_timing.csv` and the README. No cut to the
  sample is needed.
- Arm A needs light curves that are not on this Mac: up to 120 stars, about 10 to 40 sectors each,
  QLP or TESS-SPOC files of roughly 0.1 to 2 MB each from MAST. Arm B needs the same for the hosts
  that are not cached.

Observed while testing the harness on synthetic noise (not a result): at low reference SNR the
threshold 0.5·depth − 3σ is negative, so a strong brightening at a predicted epoch still counts as
`excludes`. Whether that matters on real data is what Q1 and Q3 measure.
