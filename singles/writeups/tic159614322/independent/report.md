# Blind re-analysis of TIC 159614322 (2026-10-11)

A separate agent ran this analysis. It was given only the two known events (S14 and S41) and the
raw data. It wrote fresh code (`lcio.py`, `model.py`, `phot.py`, `s1`–`s11`) and did not see the
main analysis or its Sector 21 finding. The agent could not write this file itself, so it was
transcribed from its final message. Numbers are in `results.json` and the intermediate JSONs in
this folder.

## Bottom line

**The period is narrowed to P = 179.93 d (n = 4) or 89.96 d (n = 8).** Existing TESS data cannot
separate the two.

This rests entirely on one partial dip in Sector 21 (BTJD ~1883.77), at an orbit end, with the
egress lost. If that dip is rejected, n = 1, 2 and 5 (719.7, 359.9 and 143.9 d) also survive.
The honest claim is "narrowed to {179.93, 89.96} d", not "resolved".

## Data and the two known events

**Sectors observed:** 14, 15, 16, 21, 22, 23, 41, 48, 49, 75 and 76.
- All have TESS-SPOC and QLP products.
- Sectors 41, 48 and 49 also have SPOC 2-min data.
- tess-point shows no further TESS coverage until Sectors 120–128.

**Joint fit of S14 (TESS-SPOC) and S41 (SPOC 2-min)**, with its own ensemble MCMC:

| Quantity | Value |
|---|---|
| Tc1 | 1703.8284 ± 0.0039 BTJD |
| Tc2 | 2423.5480 ± 0.0024 BTJD |
| dT | 719.7197 ± 0.0047 d |
| Depth | 3.00 (+0.14/−0.10) ppt |
| T14 | 9.79 ± 0.14 h |
| b | 0.30 ± 0.21 |

## Alias table

Depths are fixed-shape fits at the predicted time. A real transit would read about 3.0 ppt.

| n | P (d) | Epochs tested | Verdict | Evidence |
|---|---|---|---|---|
| 1, 2, 5 | 719.72, 359.86, 143.94 | none in data | Allowed by coverage | No third epoch falls in the data. A transit at the S21 dip would be inconsistent with these periods. |
| 3, 6, 9 | 239.91, 119.95, 79.97 | S23, S49, S76 | **Excluded** | S49 is fully covered and flat in four reductions (about 0 ± 0.15 ppt). S76 is flat. At S23, QLP and own photometry are flat. |
| 4 | 179.930 | 1883.758 ± 0.003 (S21) | **Supported** | S21 dip |
| 7 | 102.817 | S22, S48, S75 | **Excluded** | Flat, about 0.0 ± 0.15–0.2 ppt, in 2–3 products plus own photometry. |
| 8 | 89.965 | S21 only | **Supported**, but cannot be told apart from n = 4 | Its other epochs all fall in gaps. |
| 10 | 71.972 | S22, S49, S75 | **Excluded** | Flat. |
| 11 | 65.429 | S22, S48 | **Excluded** | S48 fully covered and flat in three products. |
| 12 | 59.977 | the n = 3 epochs, plus S21 | **Excluded** | Flat at the n = 3 epochs. |

The multiples of 4 from n = 16 to 48 also predict the S21 time. Each has 3–13 fully covered epochs, all flat, so all are excluded.

## The S21 dip

**Coverage.** The ingress and the floor of the transit are covered. Data end at 1883.971, while the predicted end of transit is 1883.962. Background ramps up after ~1883.90, so the egress is effectively lost.

**Depth at the predicted time:**

| Reduction | Depth (ppt) |
|---|---|
| TESS-SPOC PDCSAP | 2.26 ± 0.27 |
| TESS-SPOC SAP | 3.27 ± 0.33 |
| QLP | 3.62 ± 0.32 |
| Own TESScut photometry | 3.60 ± 0.37 |

**Timing.** An ingress-side fit gives Tc3 = 1883.767 ± 0.007 (the products span 1883.740–1883.768). That is 0.22 h, or 1.2σ, later than the n = 4/8 prediction.

**For a real transit:**
- the background is flat through ingress;
- the depth falls as the aperture grows, consistent with an on-target source;
- sliding the transit model across all 462 product-days of FFI data, the S21 dip is the only one seen in two products.

**Against:**
- Pixel-level localisation is inconclusive, because of a pointing-drift pattern (χ² 135 on 19 px).
- The in-transit centroid shift of (+1.8, +0.9) mpx is not reliable, since only a pre-transit baseline was available.
- The Gaia neighbour check was not run.
- All the products come from the same FFI pixels.

**Duration versus stellar density** (circular orbit): the fit penalty is Δχ² = 14 for n = 4, 4.7 for n = 8 and 66 for n = 1. This weakly favours the shorter periods.

## Conclusion

**Excluded:** n = 3, 6, 7, 9, 10, 11, 12, and the multiples of 4 from 16 to 48.

**Remaining:**
- n = 4: P = 179.9295 ± 0.0011 d
- n = 8: P = 89.9648 ± 0.0006 d

Linear ephemeris through all three transits: T0 = 1703.8303 ± 0.0036 BTJD (χ² = 1.5).

**Upcoming transits common to both periods:** BTJD 3683.055, 3862.984, 4042.914, 4222.843 and 4402.773. The n = 8 period adds the epochs half-way between them.

**The agent's view:** the S21 dip is probably a real, on-target transit.
