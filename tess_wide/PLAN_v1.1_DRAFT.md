# Wide single-transit search: plan v1.1 — DRAFT, NOT APPROVED, NOTHING RUN

Drafted 2026-11-25, after the S104 pilot under v1.0 (`PLAN.md`, `PILOT_S104.md`). v1.0 stays the record for S104.
This draft changes the parts of v1.0 that the pilot showed don't work. Everything not named here is
unchanged from v1.0 (data, star cut, detector rules, common-mode veto, one event per star, size test,
cross-match).

**Approval steps, as for the mono and DR4 pre-registrations:**
1. Jon approves this draft. It then becomes v1.1.0 and is frozen.
   - That approval also covers downloading S99, about 45 GB.
2. The calibration in section 3 runs on S104 injections only. Its numbers are committed before any S99 light
   curve is downloaded.
3. The S99 run, then its report. Candidate lists stay private, as in v1.0.

## 1. Why v1.1, in two lines

On S104 the inverted-flux control ended with more CNN passes than the real search (300 against 225), and about
80 % of the passes on both sides carry a background or centroid shift (32 % of injected transits do). The
pipeline was finding scattered-light and pointing structure in raw SAP flux, and the CNN, which re-measures
transits but doesn't classify them, passed it.

## 2. Development and test sectors (fixed now)

- **S104 is the development sector.** Section 3 sets numbers from its injections only. Injections carry no
  information about any real event, and no S104 candidate is looked at individually. The real and control
  counts of S104 have already been seen, so S104 can't test v1.1.
- **S99 is the test sector.** The rule: the most recent QLP sector (200-s cadence) that shares under 1 % of
  S104's cut stars and has a normal length.
  - Measured with tess-point on 3,000 S104 cut stars: S99 0.0 %, S98 0.0 % (but about 55 d long, like S97),
    S100 3.0 %, S101–S103 18–54 %.
  - S99 is searched once, under the frozen v1.1. Its result is reported whatever it is.

## 3. Changes from v1.0

### 3.1 Systematics removal before the search (new)

For each camera, CCD and orbit (ORBITID), the steps are:
- **Basis:** the 500 stars with the lowest point-to-point scatter among the cut stars on that CCD. Their SAP
  flux, divided by its median and binned to 30 min, gives the first k = 4 principal components, interpolated to
  each star's cadences.
- **Fit:** each star's SAP flux is fitted with a constant, the 4 components and its own `SAP_BKG`, by least
  squares. Cadences more than 3σ off the fit (either sign) are dropped and the fit repeated, 3 times. The
  search then runs on flux / fit.
- **Injections and the control:** injected transits go into SAP before this step, so its losses are measured.
  The control is 2 − (flux / fit).

| Choice | Reason |
|---|---|
| Common basis per CCD and orbit | Scattered light and pointing jitter hit all stars on a CCD in the same orbit together. In S104 they made most of the survivors, and the common-mode veto catches only the times when many stars cross the detection threshold together, not shared trends below it. |
| 500 quietest stars | Enough that one star's own transit or variability can't dominate a component. Quiet stars carry the instrument signal, not stellar signals. |
| k = 4 | Six free parameters per orbit (constant, 4 components, background) against ~4,000 cadences per orbit, too few to bend around a 2–13 h dip. SPOC's PDC fits up to 8 components but with priors; without priors, fewer is safer. **Fixed fallback:** if, in section 3.4, the median depth retained by S104 injections after this step is below 0.9, k = 2 is used instead. Nothing else is tried. |
| `SAP_BKG` as a regressor | 119 of 225 S104 real passes showed an in-event background change. Background leaks into SAP flux. |
| 3σ clipping, both signs, 3 passes | Keeps real dips and flares from pulling the fit. 3σ and 3 passes are the conventional values; the injections measure what they cost. |

### 3.2 Background and centroid tests become cuts (were flags)

- **Statistics:** the in-event shift of `SAP_BKG`, `SAP_X` and `SAP_Y` against the flanks, as in v1.0 5.4. Each
  is now divided by its red-noise error: the flux detector's β method applied to the ancillary series.
- **Threshold:** an event is cut if any of the three exceeds τ.
  - τ is set in section 3.4 as the 97.5th percentile of the largest of the three over S104 injections recovered
    at detection.
- **Reason:** injected transits don't change background or centroid, so their statistics are exactly the null
  distribution of these tests on real stars at random times. τ at the 97.5th percentile loses about 2.5 % of
  real on-target transits by design.
  - At v1.0's 3σ, 32 % of injected transits were flagged. The red-noise term is what was missing.
- **Limit:** a real transit on a star whose centroid moves (a blend) is cut. That is intended: the search is
  for on-target transits.

### 3.3 The CNN is a measurement, not a cut

- It still runs on every survivor and its depth, duration, time and σ are recorded, but no event is removed by
  it.
- **Reason:** on S104 it passed 16 % of real single events, 17 % of control ones and 67 % of injected transits.
  It costs a third of real transits and doesn't separate systematics. A classifier vetter (transits against
  systematics, with inverted-flux events as negatives) would be a separate v1.2, with its own validation.

### 3.4 Calibration on S104 (step 2 of the approval; injections only)

1. Re-run the 3,000 S104 injections (same seed 20261125 and draws as v1.0) through section 3.1.
   - Record the median ratio of the retained to the injected depth. That applies the k fallback rule.
2. Compute the injected events' background and centroid statistics (section 3.2) and set τ.
3. Commit τ, k and the retention numbers before any S99 download.
   - The S104 real and control searches are re-run too, for the record. They don't set any number.

### 3.5 Success criteria on S99 (fixed now)

Both must hold:
- **Control rule (unchanged from v1.0):** final-stage control count < 0.5× the real count. The final stage is
  now after the size test and the section 3.2 cuts.
- **Completeness:** end-to-end recovery ≥ 30 % for injected depths ≥ 3,137 ppm (the top two v1.0 depth bins,
  averaged over durations).
  - **Reason:** v1.0 had 35–50 % there. A fix that buys a clean control by losing more than about a third of the
    transits it can see is not worth an all-sky run.

Outcomes:
- **If both hold:** the S99 candidates are looked at by eye and go to a vetting design (pixel-level checks).
  There is still no all-sky run without a further plan.
- **If either fails:** the recommendation is to stop the wide FFI search in this form, and report the two
  pilots as a methods result. The decision is Jon's.

## 4. Same S99 measurements as the S104 report

The report gives:
- stars;
- events at each stage, real and control, with a per-camera breakdown of the final stage;
- completeness grids at detection and end to end;
- runtime and disk;
- the projection.

## 5. Limits

- S104 injections set τ and the k fallback. If S99's noise differs a lot, τ may be loose or tight there; S99's
  own injections will show it.
- A basis built from cut stars on the CCD includes some variable and eclipsing stars. Taking the 500 quietest
  limits this but doesn't remove it.
- One reduction and no pixel-level checks, as in v1.0.
