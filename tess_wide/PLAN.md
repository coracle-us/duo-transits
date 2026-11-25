# Wide blind search of TESS FFI light curves for long-period single transits

Version 1.0, 2026-11-25. **Fixed before any light curve of the search was downloaded or any candidate was looked
at.** Any later change gets a dated entry under "Changes" with its reason, and the pilot reports results under
both the old and the new rule.

Scope of this version: the plan, and a pilot on one sector (S104). No submission, email or post. **No candidate
list is made public**: the repository is public, so files that name candidate stars stay under `tess_wide/data/`
(a gitignored link to `duo-transits-cache/tess_wide`, outside git). Committed files carry counts only.

## 1. Goal

Find single transit-like events, one per star per sector, of the kind a planet on a period longer than about a
sector makes, on stars bright enough for radial-velocity follow-up. Remove events that are already known (TOIs,
CTOIs, known eclipsing binaries, the Hawthorn et al. 2024 duotransits). Measure completeness with injected
transits instead of assuming it. This version runs one sector, to measure yield, completeness, runtime and disk,
and to project the cost of all sectors.

## 2. Data

- **QLP FFI light curves** from MAST, the HLSP at doi:10.17909/t9-r086-e880, one FITS file per star per sector.
  - **Why QLP, not TESS-SPOC:** on 2026-11-25 MAST has QLP bulk scripts for S1–S104 and TESS-SPOC only to S86.
    A search of every sector needs one reduction throughout.
  - QLP limits its stars to T < 13.5. The star cut (section 3) is brighter than that, so the QLP limit never
    binds.
- **Pilot sector: S104**, the most recent QLP sector on MAST (200-s FFI cadence).
- **Flux:** `SAP_FLUX` where `QUALITY == 0`. The first 0.5 d of each orbit is dropped, as the K2 search did for
  QLP.
  - **Why SAP, not KSPSAP:** QLP's KSPSAP flux is detrended with a spline, which can absorb part of a long
    transit. The detector fits its own local baseline around each trial event (section 4), so it doesn't need
    the pre-detrended flux.
  - The cost of using SAP is more systematics. The inverted-flux control (section 4) and the common-mode veto
    (section 5) measure and remove them.
- **Kept columns:** time, `SAP_FLUX`, `QUALITY`, `SAP_BKG`, `SAP_X` and `SAP_Y`. They are kept per star as a
  compressed `.npz` and the FITS file is then deleted.
  - **Why:** a whole sector of QLP is 1,152,719 files × 682,560 B ≈ 787 GB, and keeping every FITS file for 104
    sectors is not reasonable on one laptop. The `.npz` keeps everything the search, the vetter and the
    injections need.

## 3. Star cut

The stars must be in the S104 QLP file list. The TIC 8.2 rows come from the MAST TIC TAP service, and a star is
kept if all of these hold:

| Rule | Reason |
|---|---|
| `objtype = STAR`, `disposition` empty (not SPLIT, ARTIFACT or DUPLICATE) | Removes TIC artefacts and duplicates. |
| `Tmag ≤ 12.0` | Follow-up: a single transit only becomes a planet with RVs, and warm giants at T ≈ 12 are within reach of the 1–4 m RV programmes that follow up TESS singles (Hawthorn et al. 2024 and Salinas et al. 2026 work at T ≈ 12–13). It also limits download cost, which grows steeply with depth: QLP has 1.15 million stars to T = 13.5. |
| `rad ≤ 2.0 R☉` (required) | On a larger star a 1 R_J planet gives under 0.25 % depth, and giant-star granulation dominates. Subgiants up to 2 R☉ stay in. Stars without a TIC radius are excluded, because the size test (section 5) needs one. |
| `3000 K ≤ Teff ≤ 7500 K` | The range where TIC radii are calibrated and RV follow-up is practical. Hotter stars rotate fast and have few lines; cooler ones are almost all fainter than T = 12. |
| `logg ≥ 3.5` where given | Excludes giants whose radius is missing or wrong in the TIC. |
| `contratio ≤ 1.0` where given | The target gives at least half the light in the aperture. With more dilution than 2×, the true depth, and so the size test, is uncertain by more than 2×. |

The cut uses only catalogue values, never light curves. The resulting count goes in the report.

## 4. Single-event detector

The detector is `singles/revision_search/search_any.py`, **imported, not copied**: the same code and rules as the
AAS82207 revision and the K2 blind search (AAS84033). Its rules are in its docstring. The ones that matter here:

- free depth, with a box tested at each duration on a D/4 grid;
- a local quadratic baseline from flanks of max(0.75 d, 2D) on each side;
- SNR with a red-noise factor β measured at the in-event timescale;
- **detection at SNR ≥ 8**, with coverage ≥ 0.8, n_in ≥ 3, baselines on both sides, edges ≥ 0.25 d from
  gaps, a momentum-dump test, median depth ≥ 0.5× the mean depth, and depth ≥ 3× the 1-h scatter.

**Why reuse it unchanged:** its thresholds were set and published before, with reasons (AAS82207 revision Q5). It
recovered the catalogued events and the Prša et al. 2022 positive-control EBs. Re-tuning it on this data would be
exactly the "changed after the run" problem.

Changes from `search_any`'s defaults, fixed here:

| Setting | Value | Reason |
|---|---|---|
| Duration grid | 2, 3, 5, 8, 13 h | Log-spaced by about 1.6×, so a true duration is never more than ~25 % from a grid box. 2 h is about the shortest central transit for P ≈ 27 d around an M dwarf; grazing transits shorter than that are left out. 13 h is the top of the CNN vetter's training range (1.5–15 h, section 6). The 24-h box needs 2-d flanks on each side of a 13-d orbit and is dropped. |
| Search per product | QLP only | S104 has no other pipeline light curve, so the "second reduction" re-measurement of `search_any.annotate` does not run. This is stated as a limit. |
| Momentum dumps | QLP `QUALITY` bits 32, 1 and 4 of the star's own file | The same bits as `search_any`; no TESScut cutouts are downloaded. |

**Inverted-flux control:** the identical search, with every later stage of section 5 and the vetter, runs on 2 −
flux for every star. The control counts are reported at every stage beside the real ones. **Rule fixed now:** if,
at the final stage, the control count is ≥ 0.5× the real count, the pilot reports that the search is not yet
separating signal from systematics. No threshold is then changed in this version.

## 5. From detections to single-event candidates (in order)

1. **Common-mode veto.** For each camera and CCD, event centres from all searched stars (real search) are counted
   in 0.5-d bins. A bin is flagged when the Poisson probability of its count, given the CCD's mean count per bin
   (with flagged bins left out of the mean, one iteration), is below 10⁻⁶. Events centred in a flagged bin are
   vetoed.
   - **Reason:** real transits on different stars are uncorrelated in time, while scattered light and pointing
     glitches hit many stars on a CCD at once. S104 has 16 CCDs × about 55 bins ≈ 900 bins, so a 10⁻⁶ level
     expects about 0.001 false flags per sector.
   - The control uses its own flagged bins, computed the same way.
2. **One event per star.** A star is kept only if exactly one event survives step 1.
   - Two events of consistent depth and duration (both within 0.5–2×) count as "repeating" (P less than a
     sector): out of scope here, counted and reported.
   - Any other star with ≥ 2 events is "multi-event" and dropped. Such a star is variable or affected by
     systematics, and a single real transit on it can't be told apart at this stage.
   - **Reason:** a planet with P longer than the sector makes at most one transit in it. The 0.5–2× ratio is the
     same tolerance `search_any` uses to confirm an event in a second product.
3. **Size test.** Rp = R★ √(depth (1 + contratio)) with the TIC R★ and contratio (taken as 0 if missing). The
   event is dropped as an eclipsing binary if Rp > 2 R_J.
   - **Reason:** the largest known planets and brown dwarfs are under about 2 R_J, the same rule as the K2 leads
     triage.
   - The box depth runs a little shallow, so the rule errs on the side of keeping an event.
4. **Flags, not cuts** (recorded for every surviving event, used only to describe the survivors):
   - in-event change in `SAP_BKG` above 3σ of its baseline scatter (scattered light);
   - in-event shift of `SAP_X` or `SAP_Y` above 3σ (a moment-centroid shift, with the implied source offset
     −Δc(1 − d)/d in pixels, as in the K2 leads triage).
   - **Reason they don't cut:** neither test has been validated on QLP FFI data here. A cut needs a validated
     threshold; a flag doesn't.

## 6. Neural-net vetter: the Salinas et al. 2026 CNN

- **What it is:** `helemysm/singletransiters` (arXiv:2610.24882) is a 1-D CNN that estimates depth, duration and
  mid-time with calibrated uncertainties from a 99.5-h window around a given time. It was trained on synthetic
  single transits with red noise and stellar variability: depth 300–15,000 ppm, duration 1.5–15 h, b ≤ 0.9.
- **It does not classify** planet against non-planet. It is used here as an independent re-measurement: an event
  passes if a network trained only on transits measures a transit of the same size, at the same time, at high
  significance.
- **Running it:** the existing local clone and environment in the main checkout
  (`salinas-second-transits/singletransiters`, `.venv-cnn`, Python 3.10, torch 2.14) is called from our wrapper
  through its `computeparams.param_inference_core` functions, as `salinas-second-transits/st2_pipeline.py` does.
  - The repository has no licence file, so **none of its code or weights is copied into this repository**; the
    wrapper imports it from that path.
  - Input is our QLP SAP light curve for the sector (`QUALITY == 0`), through their `preprocess_lc`.
- **Pass rule (all four must hold):**
  - CNN depth / σ_depth ≥ 5;
  - CNN depth within 0.5–2× the box depth;
  - CNN duration within 0.5–2× the box duration;
  - |CNN mid-time − box time| ≤ max(1 h, D_box / 2).
- **Reasons:**
  - 5σ asks for a depth measured significantly by a second, independent method; the injections (section 7) show
    how many real transits it loses, and the control shows how many systematics it lets through.
  - The 0.5–2× tolerances are the detector's own confirmation ratio.
  - The time window keeps the CNN on the same event (it can move its centre within its 99.5-h window).
- **Out of range:** events with box depth > 15,000 ppm or outside 1.5–15 h are outside the training range. They
  are kept, labelled "CNN out of range", and counted separately rather than passed or failed.

## 7. Injection–recovery

- **Stars:** 3,000 stars drawn at random from the searched sample (seed 20261125), one injection each, into the
  SAP flux before any processing.
- **Transits:** quadratic limb-darkened models (Mandel & Agol, `batman`) with:
  - depth log-uniform in 300–15,000 ppm (as Rp/R★, b uniform 0–0.9);
  - duration log-uniform in 2–13 h;
  - mid-time uniform over the star's good data (an injection that falls in a gap counts as not recovered,
    because that is part of real completeness).
- **Recovered at detection:** an event within max(0.1 d, D_inj / 2) of the injected time. **End to end:** the
  same event also survives sections 5.1–5.3 and the CNN pass rule.
  - The section 5.1 flagged bins are those of the real run.
  - Injections with a true Rp > 2 R_J (from the injected depth and the TIC R★) are left out of the end-to-end
    denominator, because the size test removes them by design.
- **Reported:** completeness in a 5 (depth) × 4 (duration) grid of bins at detection and end to end, with binomial
  errors.
- **Reason for 3,000:** about 150 per bin, for errors of ±4 % at 50 % completeness.

## 8. Cross-match: keep only what isn't already known

A surviving event is "known" if any of these match:

| Catalogue | Match | Source |
|---|---|---|
| TOIs (any disposition) | same TIC | ExoFOP TOI table, downloaded on the day of the run |
| CTOIs | same TIC | ExoFOP CTOI table, same day |
| TESS EB catalogue (Prša et al. 2022) | same TIC, or within 42″ | VizieR J/ApJS/258/16 |
| Gaia DR3 eclipsing binaries | within 42″ | VizieR I/358/veb (the ESA archive is down for DR4 from 2026-11-25) |
| VSX eclipsing types (E, EA, EB, EW and their subtypes) | within 42″ | VizieR B/vsx |
| Hawthorn et al. 2024 duotransits | same TIC | `data/catalog/candidates.csv` |
| Salinas et al. 2026 singles | same TIC | `salinas-second-transits/candidates.csv` |
| Confirmed planets | same TIC | NASA Exoplanet Archive `pscomppars`, same day |

- **Why 42″ (2 QLP pixels):** an eclipsing binary that close falls in or next to the QLP aperture and can make the
  dip. TIC matches are exact.
- The remaining events are the pilot's **new single-event candidates**. They are listed only under
  `tess_wide/data/` and are not made public. Eyeballing them is a step after this pilot; it doesn't feed back into
  any threshold of this version.

## 9. Measurements and report

The report (`tess_wide/PILOT_S104.md`, counts only) gives:
- stars in QLP S104, after the star cut, and searched;
- events at each stage, real and inverted control: detection, common-mode veto, one per star, size test, CNN pass
  (and out of range), cross-match;
- completeness by depth and duration, at detection and end to end;
- runtime per stage, files and bytes downloaded and kept;
- a projection of the same for all QLP sectors (S1–S104) on this MacBook: stars per sector scaled from the
  pilot, per-star costs from the pilot.

## 10. Known limits, stated before the run

- There is one reduction (QLP) per sector, so events aren't confirmed in a second product.
- There are no pixel-level (difference-image) checks in this version; the moment-centroid flag is weaker.
- The CNN re-measures and doesn't classify. Blended eclipsing binaries smaller than the size cut pass it.
- The 0.5-d orbit-start trim and the edge rule lose transits near gaps; the injections include that loss.
- The S104 TOI and CTOI tables change daily. The cross-match records the download date.

## Changes

(none yet)
