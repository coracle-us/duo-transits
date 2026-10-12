# EPIC 211924561 / TIC 238606535: independent blind re-analysis

All code in this folder was written from scratch for this analysis. Data were downloaded with lightkurve into `data/` (not committed).
Times are BTJD = BJD_TDB − 2457000 (K2 BKJD − 2167).

**Coordinates:** the RA/Dec in the brief (131.0849, +20.3497) are wrong. The K2 and TESS headers and Gaia DR3 660640844292494336 put the star at RA 132.65255, Dec +19.09137 (G = 10.70, FLAME R = 1.739 R☉, M = 1.424 M☉, so ρ = 0.382 g cm⁻³; Teff ≈ 6680 K).

## Data searched
- **K2:** C5, C16 and C18. Reductions: EVEREST, K2SFF and K2-pipeline SAP/PDC; C5 also K2SC, C18 also 1-min.
- **TESS:** S44, S45, S46 and S72. Reductions: SPOC 2-min, TESS-SPOC, QLP and TARS. There is also a 7-day TESS-SPOC product labelled "Sector 1751" (BTJD 4055.8–4063.0).
- tess-point and TESSCut list no other sectors.

## Events (joint 4-event fit, primary reductions; `fit_results.json`)
| event | reduction | Tc (BTJD) | alt. reduction Tc |
|---|---|---|---|
| K2 C5 | K2SFF | 179.89779 ± 0.00077 | EVEREST 179.89812 ± 0.00093 |
| K2 C16 | EVEREST (outlier bits ignored) | 1163.35544 ± 0.00131 | K2 SAP 1163.35496 ± 0.0011 |
| TESS S46 | SPOC 2-min PDCSAP | 2568.29833 ± 0.0030 | QLP 2568.2940 ± 0.0043 |
| TESS S72 | SPOC 2-min PDCSAP | 3270.76459 ± 0.0032 | QLP 3270.7672 ± 0.0038 |

- **Shape:** T14 = 3.00 ± 0.07 h. Single-event depths are consistent within errors (K2 2.2 ppt, TESS 2.2–2.3 ppt), so the depth is achromatic within about 10%.
- **C16 event:** EVEREST sets outlier bits on its in-transit cadences, so its QUALITY = 0 data hide it, and K2SFF has a gap at 1161.8–1166.1. The event is clear in the raw K2 SAP flux.
- **Blind prediction:** the ephemeris from C5, S46 and S72 alone predicted the C16 transit at 1163.3555; it is observed at 1163.3554.
- **Vetting:**
  - TESS events: present in SAP and PDCSAP, TESS-SPOC and QLP (weak in TARS). They are not near momentum dumps or background features.
  - Difference-image centroid offset: (0.00, −0.10) ± 0.10 px in S46 and (−0.15, −0.06) ± 0.10 px in S72. This rules out the G = 16.05 star at 25″ (1.2 px) and the G = 16.6 star at 40″. The other neighbours (G > 20) are too faint to make a 2 ppt dip.
  - K2 C5: faint pixels show no additive deficit, which rules out a rolling-band artefact.
- **Rejected candidates:** S44 2509.8 (scattered-light background step, opposite sign in QLP), S45 2539.5 (QLP only, at a gap edge), S72 3275.5 (TESS-SPOC PDC edge only), C5 EVEREST 145.17 and C16 EVEREST 1114.8 / 1116.5 (EVEREST only).

## Periods (`alias_table.json`, s14)
- **Ephemeris-consistent periods:** P = 3090.867/N with N a multiple of 22, i.e. P = 140.494/m (140 such periods with P ≥ 1 d).
- **Surviving periods:** every m ≥ 3 is excluded by at least 2 clearly flat epochs, each confirmed in 2–3 reductions. For m = 3 the excluding epochs are C18 1257.02 and S44 2521.46; for m = 4 they are C5 144.77, C16 1128.23, C18 1268.73 and S45 2533.17.
  - **P = 140.49401 ± 0.00011 d**, T0 = 179.89774 ± 0.00073 (epochs 0, 7, 17, 22)
  - **P = 70.24701 ± 0.00006 d**, same T0 (epochs 0, 14, 34, 44)
  - For P = 70.25 d, the predicted transits fall 1.5–2.4 d outside the C16, C18 and S44 data and miss the "S1751" window, so nothing tests it.
- **Secondary eclipse:** for P = 70.25 d (circular orbit), phase 0.5 was observed at 4 epochs and is flat: −0.04 ± 0.05 ppt (< 0.15 ppt at 3σ).
- **Next epoch that distinguishes the two periods** (70-d only): BTJD 4324.47. BTJD 4254.22 is shared by both.

## Shape and stellar density (s15, s15b)
- **Shape:** V-shaped and grazing. The profile χ² prefers b ≥ 0.95 (T23 = 0). Flat-bottomed fits with b ≤ 0.9 are worse by Δχ² ≈ 19–21, assuming white noise.
- **k is unconstrained upward:** k ≈ 0.055–0.16 for b = 0.95–1.10, i.e. Rp ≈ 0.9–2.7 RJ, or larger.
- **Duration vs density:** a central transit of this star would last 15.2 h at P = 140.5 d and 12.1 h at 70.2 d. The observed 3.0 h needs a grazing chord, and even then the circular-orbit density is too high:
  - **P = 140.49 d:** ρ_circ = 2.6 (2.0–4.9, 95%) g cm⁻³, 5–13× the Gaia value. Needs e ≥ 0.50 (95% lower bound; ω = 90°). A circular-orbit fit tied to the Gaia density is worse by Δχ² = 265 and pushes b to 1.27 and k to 0.31.
  - **P = 70.25 d:** ρ_circ = 1.3 (1.0–2.5) g cm⁻³, 2.7–6.4× the Gaia value. Needs e ≥ 0.31. The density-tied fit is worse by Δχ² = 69.

## Conclusion
- **Detections:** four transits on target (two in K2, two in TESS) with consistent depth, duration and shape. The planetary-period set is {140.494, 70.247} d; every shorter alias is excluded with flat data from at least 2 independent reductions. I am highly confident the true period is one of these two, and nothing in the data prefers either.
- **Caveat on the nature of the object:** the short, V-shaped transit does not fit a circular orbit around the Gaia/FLAME star. It needs a high-eccentricity orbit transiting near periastron (e ≳ 0.5 or ≳ 0.3), a grazing large companion (possibly stellar, i.e. a grazing eclipsing binary), or a host other than the primary (an unresolved bound companion; RUWE = 1.19). Treat this as a candidate, not a planet, until RV or high-resolution imaging is available.
- **Other caveats:**
  - Fixed quadratic limb darkening.
  - The white-noise χ² slightly overstates the shape preference.
  - The b = 1.15 profile point failed to optimise and was ignored.

Scripts: s01–s18 (`lcio.py` handles I/O and detrending, `tfit.py` the model). Key figures: `fig_fit_joint4_primary.png`, `fig_vet_S46_2568.30.png`, `fig_vet_S72_3270.75.png`, `fig_c16_predicted_epoch.png`, `fig_tess_diffimg_S46.png` / `fig_tess_diffimg_S72.png`, `fig_exclusions_m3_m4.png`.
