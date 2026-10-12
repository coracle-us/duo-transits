# TODO on Gaia DR4 release day (2 December 2026)

The date is from ESA's release page (https://www.cosmos.esa.int/web/gaia/release), checked
2026-10-11.

The DR4 content page (https://www.cosmos.esa.int/web/gaia/dr4) lists these epoch products:
- `rvs_epoch_parameters_single` (single-lined transits);
- `rvs_epoch_parameters_double` (**double-lined** transits);
- `rvs_epoch_spectrum`;
- `epoch_astrometry`.

That page does not state a G_RVS limit for the epoch RVs. "G_RVS ≤ 14" comes from the second
reviewer and is not verified here. Both stars below are much brighter than that.

Jon may also want to forward this to the gaia-dr4-companions project, which the reviewer
mentioned. It is not in this repo.

## 1. TIC 147035057 (Gaia DR3 2947647054662200192, G_RVS 8.28)

Our TESS ephemeris (`singles/results/vetting/147035057_ephemeris.json`):
**P = 37.31497 ± 0.00003 d, T0 = BJD_TDB 2458515.0661 ± 0.0011** (primary eclipse).

1. Fold the epoch RVs on P = 37.31497 d.
   - Fit a Keplerian with P fixed.
   - Report K, e, ω and the minimum companion mass.
   - DR3 gives an RV amplitude of 49.8 km/s peak-to-peak. That implies K ≈ 25 km/s and
     M₂ ≥ 0.64 M☉ (`singles/writeups/tic147035057/gaia_tension.json`).
2. Check the eclipse timing. The primary eclipse must fall at the right RV phase: the
   conjunction with the companion in front.
3. **SB2 test.** Does the star appear in `rvs_epoch_parameters_double`?
   - Double-lined transits, with two components about 86 km/s apart at quadrature (twins
     with e = 0), would favour the SB2-twin explanation for the low RUWE (1.12).
   - Single-lined transits with K ≈ 25 km/s would favour an F star with a K/M-dwarf companion.
4. Use `epoch_astrometry` to look for a 37.3-d photocentre wobble.
   - We expect about 0.5 mas if the companion is dark.
   - We expect about 0 for equal twins.
5. Use the result to update `singles/writeups/tic147035057/` and the RNAAS draft.

## 2. TIC 159614322 (Gaia DR3 1669046625055549312, G_RVS 10.92)

Our ephemeris (`singles/writeups/tic159614322/ephemeris.json`) gives
**P = 89.9651 ± 0.0006 d or 179.9302 ± 0.0012 d, T0 = BJD_TDB 2458703.8277 ± 0.0041**.

DR3 lists 22 RV transits, an RV error of 0.38 km/s, rv_chisq_pvalue 0.88 (not variable) and RUWE 1.51.

1. Fit the epoch RVs on both periods, with K free. Two outcomes:
   - **A planet** (K ≈ 35–40 m/s for 1 M_Jup) is far below Gaia's per-epoch precision of
     roughly 1–3 km/s at G_RVS 10.9. So the data **cannot detect a planet**.
   - **A stellar companion** at 90 or 180 d would give K of several to tens of km/s, which
     the epoch RVs can **rule in or out**. Report the upper limit on K, and the
     corresponding upper limit on M₂ sin i, for each period.
2. Use `epoch_astrometry` to check RUWE 1.51. Is there a wobble at 90 or 180 d? If not, it
   points to a wide companion or noise.
3. Update `singles/writeups/tic159614322/` and the email or RNAAS text.
