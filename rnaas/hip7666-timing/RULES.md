# HIP 7666 timing: rules fixed before any timing (2026-11-24)

This file is committed before any eclipse time or pulsation phase is measured. Results are
reported whatever they show. Any change made after timing starts is logged at the bottom with
the reason.

**Target:** HIP 7666 = TIC 308953703 = Gaia DR3 407170600540598144 (RA 24.6736°, Dec +52.5188°).

**Question:** do eclipse times and pulsation phases across all TESS sectors show a light-travel-time
effect (LTTE) from a third body? If not, what limit do they set?

## 1. Sectors and light curves

From `lightkurve.search_lightcurve` and `tess-point` (run 2026-11-24):

| Sector | Dates | Products at MAST |
|---|---|---|
| S18 | 2019 Nov | TESS-SPOC 1800 s, QLP 1800 s, TARS 1800 s |
| S58 | 2022 Oct–Nov | SPOC 120 s, TESS-SPOC 200 s, QLP 200 s, TARS 200 s |
| S85 | 2024 Oct–Nov | TESS-SPOC 200 s, QLP 200 s, TARS 200 s |
| S123 | future (tess-point) | not yet observed; not used |

**Choice per sector**, following the brief's order: SPOC 120 s first, then the TESS-SPOC FFI
product, then QLP.
- S58: SPOC 120 s PDCSAP.
- S85: TESS-SPOC 200 s PDCSAP.
- S18: TESS-SPOC 1800 s PDCSAP.
- QLP for the same sectors is a cross-check of the eclipse times only: same method, compare the
  sector-mean O−C.
- TARS is not used, being outside the brief's list.
- Lee et al. used SAP; PDCSAP is used here and SAP is checked for S58 only.

**Cuts and units**
- Quality: lightkurve's default bitmask. NaNs dropped. Times in BTJD (BJD_TDB − 2457000),
  as delivered.
- Normalise per sector by the median. Remove slow trends with a linear fit per orbit-long
  stretch to out-of-eclipse points, which leaves eclipses and pulsations untouched.

## 2. Binary and pulsation models, and their separation

- **Reference ephemeris:** Lee et al. Eq. 1, BJD 2451769.56299 + 2.37222188 E, so T0 = 4769.56299
  BTJD. All O−C against Lee are relative to this. Secondary eclipses are at phase 0.5 (circular
  orbit, e = 0.003 ± 0.002).
- **Eclipse windows:** |phase| < 0.044 around 0 and 0.5, which is Lee et al.'s out-of-eclipse
  boundary.
- **Iteration on S58** (two passes):
  1. Fold on the ephemeris and build a binary template: median in bins of 0.001 in phase, then
     smoothed with a cubic smoothing spline. The pulsations are not commensurate with P_orb, so
     they average out of the fold.
  2. Subtract the template. Fit the pulsation model to the out-of-eclipse residuals: Lee's 15
     frequencies minus f2 (= f_orb, which belongs to the binary), with frequencies fixed and an
     amplitude and phase free for each.
  3. Subtract the pulsation model from the raw light curve and rebuild the binary template.
- **Template in other sectors:**
  - The S58 template is the eclipse shape for every sector, integrated over each sector's
    exposure time (S18 1800 s, S85 200 s). Depth is free per eclipse, to allow for each
    product's own dilution.
  - In S18 and S85 the same pulsation frequencies are fitted to that sector's out-of-eclipse
    residuals, with amplitudes and phases free, before timing.
  - In S18 the pulsation model is integrated over 1800 s. f1 = 25.47 d⁻¹ is above the 24 d⁻¹
    Nyquist frequency there.
- **Known approximation:** the primary pulsates and is partly covered at primary eclipse, so a
  full-amplitude pulsation model is slightly wrong during eclipses. This is stated, not
  modelled.

## 3. Eclipse times

- **Fit:** each eclipse is fitted on its own over |phase| < 0.06. Free parameters: time shift,
  depth scale, and a linear baseline (offset and slope). The model is the template, with
  pulsations already subtracted from the data. Weighted least squares, scipy `least_squares`.
- **Coverage rules:** an eclipse is used only if at least 80% of the expected cadences in
  |phase| < 0.044 are present and both ingress and egress have data. Otherwise it is dropped and
  listed.
- **Errors:** block bootstrap of the fit residuals, with 0.05 d blocks (about one pulsation beat,
  to keep correlated residuals) and 500 resamples. The error is the standard deviation of the
  refitted times.
- **Error scaling:** in each sector and eclipse type, if the eclipse times about a linear fit have
  reduced χ² > 1, that sector's errors are scaled by √χ²_ν. The factor is reported. The cost is
  that it absorbs any signal shorter than a sector, which a dynamically stable tertiary
  (P3 ≳ 5 P_orb ≈ 12 d) could only marginally produce.
- **Primary vs secondary:** fitted separately and together.
  - Apsidal motion would move the two in opposite directions and LTTE would move them together.
  - Reported:
    - each type's own linear period;
    - the secondary − primary offset per sector.
  - The offset is compared with the inner-orbit Rømer delay predicted from Lee's K values,
    P (K_B − K_A) / (π c) ≈ 9.4 s (Kaplan 2010).
- **O−C** against Lee Eq. 1, and against our own best linear fit to all TESS eclipses.
- **QLP check:** repeat the S18/S58/S85 sector-mean O−C with QLP. A difference larger than 3σ is
  reported as a systematic.

## 4. Pulsation phases and time delays

- **Modes:** the three strongest modes other than f_orb: f1 = 25.4656, f4 = 22.4278 and
  f3 = 21.1939 d⁻¹.
- **Segments:** one orbital cycle (2.372 d) each, instead of the brief's 1–2 d.
  - One cycle averages out the primary's own inner-orbit light-time signal (a_A sin i / c ≈ 10 s)
    and any orbital-phase modulation of tidally affected modes.
  - It is still short enough to follow a tertiary with P3 ≳ 10 d.
  - Segments with less than 60% of the expected out-of-eclipse cadences are dropped.
- **Per segment:**
  - Eclipses are masked (|phase| < 0.044).
  - The binary template is subtracted, and so are the other 11 frequencies, at their sector-fit
    amplitudes and phases.
  - Then the amplitude and phase of each of f1, f3 and f4 are fitted with the frequency fixed.
  - τ = −Δφ / (2π f), with Δφ relative to the sector-wide phase.
  - Errors come from the least-squares covariance, scaled by the residual rms. The weighted
    mean τ over the three modes is also given.
- **Within-sector trends:**
  - A frequency error δf makes τ drift linearly through a sector.
  - So each sector's τ series has its linear term removed (equivalent to refitting f).
  - Only curvature within a sector is a signal.
- **Across sectors:**
  - A phase connection between sectors is tried only for S58–S85. It needs a single frequency
    alias that is preferred by Δχ² > 10 over the next alias for each of the three modes, with the
    three giving the same τ within 2σ.
  - If that fails, the report says pulsation timing does not constrain delays between sectors,
    and only eclipse timing does.
- **Validation check, not a result:** fold S58 τ on the orbital phase in 8 bins, with segments of
  0.3 d for this check only. Compare with the predicted inner-orbit delay of the primary,
  amplitude a_A sin i / c from K_A = 93.16 km/s.

## 5. Models and choice

Fitted to the eclipse O−C, with primary and secondary together and a free secondary offset:

| Model | Free parameters |
|---|---|
| M0 | linear ephemeris (T0, P) |
| M1 | quadratic (T0, P, dP/dE) |
| M2 | linear + circular LTTE, grid in P3 from 10 to 3,000 d, with free amplitude and phase at each P3 |

- The choice is by BIC = χ² + k ln N.
- For M2, k counts P3 as a parameter.
- A model is preferred only if ΔBIC < −10 relative to M0.
- An eccentric LTTE orbit is fitted only if M2 is preferred.
- The pulsation τ series is checked for the same signal, within the sectors it covers.

## 6. Limit if there is no signal

- **Injection-recovery** on the real eclipse times and errors.
  - For each P3 (log grid, 10–3,000 d) and LTTE semi-amplitude A, inject a circular LTTE at
    1,000 random phases.
  - Refit M0 and count a detection when the M0 χ² exceeds its 99.73% quantile for the dof.
  - The limit A_95(P3) is the smallest A detected in at least 95% of phases.
- **Minimum mass:** A_95 is turned into a minimum tertiary mass with the mass function
  (m3 sin i3)³ / (M12 + m3)² = 4π² (A c)³ / (G P3²), taking sin i3 = 1 and
  M12 = 1.822 + 1.246 = 3.068 M☉ (Lee et al.).
- **Caveat:** for P3 longer than the TESS baseline (S18–S85, 1,820 d), LTTE is degenerate with a
  quadratic term and the limit will be weak. That is stated.
- **Historical minima** (Escolà-Sirisi 2005; Feng 2021) are not refitted. Comparing our TESS-only
  period with Lee's 2000–2022 period is the only long-baseline test.

## Change log

Changes made after this file was committed (6010f4a). For each one, the change and why.

1. **Detrending (implementation detail; made before any timing was looked at).**
   - Change: a straight line per orbit-long stretch is replaced by a cubic least-squares spline
     with knots every 1.5 d, one per continuous chunk, fitted to the out-of-eclipse flux divided
     by the binary template. A first pass, before the template exists, uses 5 d knots.
   - Why: a line over a window that ends at a secondary eclipse would put a step in that
     eclipse. The 1.5 d spline cannot follow the 20–38 d⁻¹ pulsations, and the template carries
     the ellipsoidal and reflection shape.
2. **Ingress and egress coverage, made concrete.**
   - Change: an eclipse needs at least 50% of the expected cadences in 0.015 < |Δφ| < 0.044 on
     each side.
   - Dropped eclipses are listed in `results/dropped_eclipses_*.csv`.
3. **Common delay added (pulsations).**
   - Change: besides the per-mode phases of f1, f3 and f4, each segment also gets one delay τ,
     fitted by shifting the whole 14-frequency sector model in time (`s2_pulsation.common_tau`).
     The 0.3 d orbital-phase check uses this common τ.
   - Why: in 0.3 d segments f3 and f4 are not resolved, so a single τ is the only workable
     estimator there.
4. **Sensitivity variant "S18 from QLP" (`mixq18`), added after seeing the data.**
   - The trigger: the S18 TESS-SPOC 1800 s PDCSAP light curve has about 1% systematics over the
     first orbit, an out-of-eclipse rms of 4.5 ppt and eclipse-time scatter of about 36 s.
     QLP S18 has 0.55 ppt and about 5 s.
   - The main result stays on the main products, as written above. The variant is reported next
     to it.
5. **QLP error columns.** S58 and S85 QLP files have `DET_FLUX_ERR` instead of `KSPSAP_FLUX_ERR`.
   SAP_FLUX is still the flux used.
6. **Synthetic injection test** (`tests/test_synthetic.py`) for both clocks, including the sign of
   the pulsation delay.
