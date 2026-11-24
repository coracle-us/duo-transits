# Be star pulsations vs WISE disk phase: follow-ups (2026-11-22)

Builds on `RESULTS_FULL.md`, using the same full-run data (1,084 main-sample stars, 4,108
sectors). Both analyses are **post hoc**: they were chosen after seeing the full-run result.
- Scripts: `src/s6_freq_beta.py`, `src/s6b_inclination.py`, run with `BE_RUN=full`.
- Outputs: `freq_beta.csv`, `inclination_split.csv`, `figures/full_freq_beta.png`.
- Synthetic test: ALL CHECKS PASSED. It now also checks the log-bin merge and the
  false-positive rate of the group-difference test (0.07 for a nominal 5%).

## 1. Which frequencies grow with the disk? (`figures/full_freq_beta.png`)

**Method**
- Each saved sector spectrum (600-point log grid, 0.25–24 c/d) is merged into coarser bins.
- Each bin gets the stage-4 test: within-star β of ln(amplitude) against disk_excess,
  within-star and cross-star nulls of 10,000 draws, and a bootstrap over stars.
- Two frequency axes:
  - absolute frequency, 30 bins;
  - frequency relative to each star's dominant frequency f1 (its median over the star's
    usable sectors), 24 bins of 0.1 dex. Here the bin at f/f1 = 1 is the star's own
    pulsation group.

**Absolute frequency, main sample**

| Frequency (c/d) | β (per mag) |
|---|---|
| 0.25–1 | 0.71–0.92 |
| 1–2 | 0.59–0.76 |
| 2–5 | 0.39–0.51 |
| 5–24 | **0.31–0.37, flat** |

- Every bin has p < 0.001 under both nulls. The cross-star null 95% band is about ±0.15.
- No crowding cut: the same shape, lower overall (0.51–0.65 falling to 0.13–0.15).
- SPOC/TESS-SPOC only (329 stars): the same shape (0.57–0.84 at the low end, 0.17–0.44 at
  high frequency). Most bins are significant. The weakest are near 2.6 and 5.7 c/d
  (p_swap 0.14 and 0.20).

**Relative to the star's f1, main sample**

| f / f1 | β (per mag) |
|---|---|
| 0.07–0.7 (below the dominant frequency) | 0.73–1.16 |
| 0.9 | 0.53 |
| **1.1 (the dominant pulsation group)** | **0.35** (the minimum of the curve) |
| 1.4–14 | 0.32–0.48 |

- The SPOC-only curve has the same dip at f/f1 ≈ 1 (β 0.24, p_swap 0.08).

**Reading**
- As the disk grows, amplitude rises at *every* frequency by about +3% per 0.1 mag of W1
  excess. That is the 5–24 c/d plateau.
- On top of that, the band below the star's own pulsation frequency, and below about 2 c/d
  in absolute terms, rises about 2.5 times as fast: about +8% per 0.1 mag.
- The dominant pulsation group itself sits at the broadband level. It does not grow faster
  than the floor.
- This confirms the stage-4c ratio test frequency by frequency. What tracks the disk is
  low-frequency stochastic variability, plus a broadband rise. The coherent pulsation is not
  singled out.

**What the broadband rise could be.** It is not instrumental: the SPOC formal flux error is flat
(stage 4c). Two readings remain:
- Intrinsic stochastic variability that reaches high frequency.
- A normalisation effect: relative amplitudes scale as 1/(total optical flux).
  - A uniform positive β would need the star to get *fainter* in the optical as W1 excess
    grows. That is expected for shell (high-inclination) stars, but not for most Be stars,
    which brighten.
  - The test below does not see the inclination dependence that this reading predicts. It
    has little power, though.
- I can't separate the two readings with these data.

## 2. Does β depend on viewing angle? (`inclination_split.csv`)

**Proxy: Gaia DR3 ESP-HS `vsini_esphs`**, read-only from IRSA (`gaia_dr3_astrophys_param`,
upload join by source_id).
- It is available for 1,319 of 2,541 stars.
- Gaia `vbroad` exists for only 85 stars and is fitted with cool templates (median about
  7,000 K), so it is not used.
- The catalogue's `r_value` is the W1 brightness vs W1−W2 colour correlation, not an
  inclination proxy.

**Quality cut.** In this sample vsini_esphs depends on brightness:

| G | Median vsini (km/s) | Below 10 km/s |
|---|---|---|
| ≤ 9 | 244 | 0% |
| 11–12 | 90 | 28% |
| > 12 | 64–72 | 32–35% |

Be stars rotate fast, so this is a signal-to-noise failure. The primary split therefore uses
G ≤ 11, cut at the median vsini.

**Prediction, written into the script before running.** If the disk acts mainly through its
own optical light (dilution at low i, absorption in shell stars at high i), β should be
larger for high v sin i.

**Result: no detectable difference.**
- **Primary** (G ≤ 11, main sample): 171 high-v sin i stars (median 270 km/s) against 170
  low (median 87 km/s).

  | Metric | β high v sin i | β low v sin i | Difference | p |
  |---|---|---|---|---|
  | A1 | 0.40 | 0.60 | −0.21 ± 0.3 | 0.44 |
  | B1 | 0.89 | 1.01 | −0.12 | 0.68 |
  | red_noise | 0.83 | 1.05 | −0.22 | 0.46 |
  | flicker_rms | 0.98 | 0.90 | +0.08 | 0.80 |

  - All 7 metrics have p ≥ 0.25, and most differences have the opposite sign to the
    prediction.
  - Both groups are separately significant under both nulls.
- **No crowding cut** (214 / 214 stars): all p ≥ 0.2.
- **Sensitivity split** (all stars with vsini, median within 1-mag G bins):
  - Main sample: all p ≥ 0.47.
  - No cut: differences positive, with p 0.07–0.32 except flicker_rms at p = 0.009.
  - That is one of 28 tests, in the less reliable split, so it is not a finding.

**Power.** The 1σ uncertainty on the difference is about ±0.3 per mag. A difference as large
as the effect itself (0.4–0.6) cannot be excluded. This is a weak null: no evidence that
geometry matters, and no evidence that it doesn't.

## Bottom line for a write-up

- The robust, new-looking statement is item 1: within a Be star, TESS variability grows with
  the WISE disk excess.
  - Mostly in the low-frequency stochastic component, below the star's pulsation frequency:
    about +8% per 0.1 mag.
  - Broadband elsewhere: about +3% per 0.1 mag.
  - Not preferentially in the coherent pulsation.
- That is the figure a short note would carry.
- A literature check comes before any claim of novelty. Candidates: Labadie-Bartz et al. on
  TESS Be-star outbursts and flicker, Nazé et al., and Niwano et al. 2024 on Be X-ray
  binaries.
