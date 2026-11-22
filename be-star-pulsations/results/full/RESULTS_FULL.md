# Be star pulsations vs WISE disk phase: full run (2026-11-22)

Brief: `briefs/2026-11-04-be-star-pulsations.md` (Claude relay). Pilot: `results/pilot/RESULTS_PILOT.md`.

The full run started on Jon's "Go, relax criteria". It covers all 3,237 non-EB stars of Jia et
al. 2026 (arXiv 2611.02450). The method is the pilot's, with one fix (spline overshoot, below)
and one post-hoc stage (4c, controls). Every number comes from a file in `results/full/`
written by a script in `src/`, run with `BE_RUN=full`.

## Go / no-go (relaxed, as agreed)

| Criterion | Value | Verdict |
|---|---|---|
| Stars with ≥ 2 usable sectors whose disk_excess differs by ≥ 0.1 mag | 1,350 of 2,449 (all usable sectors); **565 of 1,084** with the crowding cut | pass (pilot: 51 / 19) |
| Noise floor: rms Δln A1 between same-star sectors with \|Δx\| < 0.03 mag | **0.54** main (3,523 pairs), 0.57 no cut; robust MAD 0.37 / 0.39 | fails the brief's 0.3. Relaxed: treated as intrinsic scatter, which the within-star fit averages over |
| Synthetic test | ALL CHECKS PASSED (`tests/test_synthetic.py`, now with an overshoot check) | pass |
| Run on the Mac | 4.5 GB cache; about 2 h wall for stages 1–2 (MAST and IRSA in parallel), 25 min for stages 3–5 | pass |

**Relaxed criteria:**
- Both the crowd-cut sample and the uncut sample are reported as primary.
- The noise floor is treated as intrinsic Be-star scatter, not a measurement limit.

## Headline

**1. Within a star, variability rises when the WISE disk is brighter.**
- Main sample (1,084 stars, 4,108 sectors), A1: β = **0.40 per mag**, 68% interval
  [0.31, 0.50]. That is about +4% amplitude per 0.1 mag of W1 excess.
- p < 0.001 under both the within-star permutation null and the cross-star disk-swap null
  (10,000 draws each).
- The low-frequency metrics rise about twice as fast:

  | Metric | β (per mag) |
  |---|---|
  | B1 (0.25–0.8 c/d) | 0.80 |
  | red_noise | 0.74 |
  | flicker_rms | 0.84 |

- The effect holds:
  - without the crowding cut;
  - after dropping mixed-pipeline sectors;
  - after dropping W1 < 8 stars;
  - in QLP stars and in SPOC/TESS-SPOC stars, except A1 in the SPOC group (below).

**2. The effect is small next to the scatter.**
- The within-star Spearman ρ is only 0.07–0.21.
- Sector-to-sector amplitude scatter at fixed disk (about 0.5 in ln) is far larger than the
  disk-driven change.

**3. It is mostly stochastic low-frequency variability, not the coherent pulsation.**
- The planned white_noise control (15–24 c/d floor) also rises with disk: β 0.33, p < 0.001.
- But the pipeline's formal flux error in SPOC/TESS-SPOC data, a purely instrumental
  quantity, is flat: β 0.09, p 0.56.
- So the 15–24 c/d floor is itself intrinsic, not a data-quality effect.
- Measured relative to that floor:
  - B1, B2, red_noise and flicker still rise: β 0.33–0.51, p ≤ 0.0002 in the main sample.
  - A1 and B3 do not: β 0.07 and 0.06, p 0.33 and 0.37.
- In SPOC/TESS-SPOC stars (329) the A1 slope itself is not significant (β 0.17, p 0.11 / 0.15).
- **Reading:** as the disk grows, the broadband red-noise or flicker component grows. The
  strongest coherent peak only rides along with the overall variability level.

**4. Spectral shape does not follow disk phase beyond time systematics.**
- The classifier's cluster–phase link is significant under the within-star permutation null
  (p 0.0005).
- It vanishes under the cross-star null (p 0.73).

## Change from the first full pass: spline overshoot (found in this run)

**The problem**
- `figures/full_examples.png` showed the FITPACK spline (s = 0.04) swinging far outside the
  measured W1 between sparse visits, typically between the last China-VO visit and a single
  IRSA extension visit.
  - It reached W1 = 30–45 mag in TIC 215343947 and TIC 348315321.
  - It oscillated by 4 mag in TIC 179637167.
  - In one star the old w1_quiet reached 5,600 mag.
- In the first pass, 97 usable sectors in 33 stars had w1_fit more than 0.1 mag outside the W1
  range of their two bracketing visits; 15 were more than 1 mag off, with a maximum of 11 mag.

**The fix** (`src/wise.py`), chosen from that figure before re-running stage 4:
- w1_quiet is clipped at the faintest real visit.
- A sector is not WISE-usable if w1_fit lies more than 0.1 mag outside its bracketing visits
  (`spline_overshoot`; `MAX_OVERSHOOT = 0.1`).
- The cross-star null applies the same rule, since it calls the same `DiskCurve.at`.
- Synthetic check: 1 of 4,000 sectors is flagged on smooth fake curves (a real bump peaking
  between visits), and an injected swing is flagged.

**Effect**
- 146 sectors dropped.
- The bad x values had *diluted* β: A1 main went from 0.24 to 0.40, and no-cut from 0.13 to 0.29.
- Phase labels changed where w1_quiet had been absurd.
- The pilot is barely affected (3 of 359 usable sectors, 1 crowd-cut), so its numbers stand.

## Sample and data (stages 0–2)

**Stage 1 draws:** all 3,237 non-EB stars (seed order 20261104).
- **2,541 pass.**
- 421 fail with fewer than 2 sectors before the last NEOWISE epoch (MJD 60519) in the four
  products.
- 275 have no TESS light curves.
- 18 runs first failed on corrupt cached FITS from an interrupted run. They passed on a resume
  pass after `retry()` was changed to delete unreadable files.

**Preferred product (stars):**

| Product | Stars |
|---|---|
| QLP | 2,020 |
| TESS-SPOC | 219 |
| SPOC 120 s | 198 |
| GSFC-ELEANOR-LITE | 104 |

**Sectors:** 10,203 in total.
- By product: QLP 8,642 / SPOC 120 s 718 / TESS-SPOC 591 / ELEANOR 252.
- 775 are flagged mixed_pipeline.

**Crowding:**
- Source: CROWDSAP for 1,309 sectors, TIC contratio for 749, the neighbour proxy for 8,145.
- Usable sectors passing the cut, by product:

  | Product | Pass / usable |
  |---|---|
  | SPOC 120 s | 530 / 689 |
  | TESS-SPOC | 364 / 561 |
  | QLP | 3,252 / 8,174 |
  | ELEANOR | 4 / 233 |

**WISE:**
- IRSA NEOWISE-R visits were fetched for all non-EB stars in 32 upload joins of 100 stars.
- 2,505 of 2,541 stars gained visits after the China-VO curve (median 1).
- Where both sources have the same visit, they agree to 0.000 mag (median).
- 9,684 of 10,203 sectors are WISE-usable. Of the rest:
  - 166 fall outside the WISE baseline;
  - 256 have a gap > 250 d;
  - 146 have spline overshoot > 0.1 mag.
- 77 stars have median W1 < 8.

**Usable sectors** (WISE-usable and metrics OK): **9,657** in 2,494 stars.
- Main sample: 4,150 usable crowd-cut sectors; 4,108 of them come from 1,084 stars with ≥ 2
  such sectors.
- Phase labels:

  | Label | All | Crowd-cut |
  |---|---|---|
  | quiescent | 3,034 | 1,430 |
  | plateau | 2,658 | 1,097 |
  | dissipating | 2,073 | 852 |
  | building | 1,892 | 771 |

- disk_excess percentiles (5/25/50/75/95) over usable sectors: 0.00 / 0.02 / 0.10 / 0.31 /
  0.87 mag.
- Median A1 is 4.3 ppt. f1 falls in 0.25–0.8 / 0.8–2 / 2–4 / 4–12 c/d for 2,678 / 3,223 /
  2,811 / 945 sectors.

## The test (stage 4): `test_results.csv`, `figures/full_pooled.png`

**Method**
- β is the change in ln(metric) per mag of disk excess, from a within-star demeaned pooled
  least-squares fit.
- The interval is a 68% bootstrap over stars.
- p-values use 10,000 draws per null.
- `p < 0.001` means at the floor (1/10,001).

**Main sample (crowd-cut): 1,084 stars, 4,108 sectors**

| metric | β (k=0) | 68% | ρ | p_perm | p_swap | β (k=0.5) | phase contrast (build − diss) |
|---|---|---|---|---|---|---|---|
| A1 | 0.40 | 0.31–0.50 | 0.10 | <0.001 | <0.001 | 0.86 | +0.17 (p <0.001 / <0.001) |
| B1 | 0.80 | 0.70–0.91 | 0.20 | <0.001 | <0.001 | 1.26 | +0.24 |
| B2 | 0.66 | 0.57–0.75 | 0.16 | <0.001 | <0.001 | 1.12 | +0.24 |
| B3 | 0.39 | 0.31–0.48 | 0.09 | <0.001 | <0.001 | 0.85 | +0.17 |
| red_noise | 0.74 | 0.65–0.84 | 0.19 | <0.001 | <0.001 | 1.20 | +0.23 |
| flicker_rms | 0.84 | 0.73–0.96 | 0.17 | <0.001 | <0.001 | 1.30 | +0.24 |
| white_noise (planned control) | 0.33 | 0.24–0.43 | 0.10 | <0.001 | <0.001 | 0.79 | +0.08 |

**No crowding cut: 2,449 stars, 9,612 sectors**

| metric | β | 68% | p_perm | p_swap |
|---|---|---|---|---|
| A1 | 0.29 | 0.23–0.34 | <0.001 | <0.001 |
| B1 | 0.59 | 0.53–0.65 | <0.001 | <0.001 |
| B2 | 0.52 | 0.46–0.57 | <0.001 | <0.001 |
| B3 | 0.26 | 0.21–0.31 | <0.001 | <0.001 |
| red_noise | 0.55 | 0.50–0.60 | <0.001 | <0.001 |
| flicker_rms | 0.62 | 0.56–0.68 | <0.001 | <0.001 |
| white_noise | 0.14 | 0.09–0.19 | <0.001 | <0.001 |

**Robustness rows** (β for A1 / B1 / red_noise / flicker_rms / white_noise, k = 0):

| Variant | Stars | A1 | B1 | red_noise | flicker_rms | white_noise | Notes |
|---|---|---|---|---|---|---|---|
| drop mixed_pipeline | 1,011 | 0.44 | 0.82 | 0.74 | 0.85 | 0.27 | all p < 0.001 |
| drop W1 < 8 stars | 1,016 | 0.49 | 0.85 | 0.77 | 0.93 | 0.30 | all p < 0.001 |
| QLP stars only | 755 | 0.52 | 0.86 | 0.78 | 0.90 | 0.29 | all p < 0.001 |
| SPOC/TESS-SPOC stars only | 329 | 0.17 (p 0.11 / 0.15) | 0.69 | 0.65 | 0.74 | 0.42 (p 0.002 / 0.014) | B1, red_noise and flicker p < 0.001; B3 0.22 (p 0.017 / 0.046) |

**Nulls and dilution**
- The cross-star null centres on 0 (swap_null_mean −0.028 to +0.015 over all 126 rows), so
  sector-era systematics shared across stars do not produce the slope.
- The dilution correction adds exactly 0.921·k to β by construction. It cannot create or remove
  the sign.

**Leverage** (`jackknife.csv`)
- Dropping any one star keeps A1 β between 0.38 and 0.46 (main) and between 0.28 and 0.31
  (no cut).
- But per-star slopes (920 stars with x range > 0.02 mag) have a median of only +0.21 for A1,
  and only 52.5% are positive.
- For B1, red_noise and flicker, 57–60% of slopes are positive, with medians of +0.53 to +0.86.
- So the pooled β is a population-level tendency carried by many stars with large disk swings.
  It is not a rule that holds star by star.

**Noise floor** (rms Δln m, main / no cut):

| Metric | Main | No cut |
|---|---|---|
| A1 | 0.54 | 0.57 |
| B1 | 0.61 | 0.65 |
| B2 | 0.54 | 0.55 |
| B3 | 0.46 | 0.48 |
| red_noise | 0.55 | 0.55 |
| flicker_rms | 0.68 | 0.74 |

These are similar to the pilot. Robust 1.4826×MAD (main): 0.32–0.54.

## Controls (stage 4c, post hoc): `control_tests.csv`, `sector_controls.csv`

Added because the planned white_noise control was not flat. Same nulls, k = 0.

**err_ppt: the median formal relative flux error per sector** (main / no cut / SPOC / QLP)

| Variant | β | p (perm / swap) |
|---|---|---|
| main | 0.25 | 0.014 / 0.061 |
| no cut | 0.13 | 0.025 / 0.087 |
| SPOC/TESS-SPOC | 0.09 | 0.56 / 0.58 |
| QLP | 0.47 | <0.001 / <0.001 |

**But QLP's KSPSAP_FLUX_ERR is not instrumental here:**
- Within a star it correlates with the light curve's own white_noise and flicker (r = 0.45 and
  0.55).
- SPOC PDCSAP_FLUX_ERR does not (r = −0.26 to −0.02).
- QLP error is consistent with being estimated from light-curve scatter.
- QLP also leaves the column empty after sector 55 (err_ppt finite for 5,794 of 8,642 QLP
  sectors).
- The clean instrumental control is therefore the SPOC one, and it is flat.

**Metric / white_noise** (β = β_metric − β_white; nulls recomputed)

| Ratio | Main | No cut | SPOC | QLP |
|---|---|---|---|---|
| A1 | 0.07 (p 0.33 / 0.53) | 0.15 (p <0.001 / 0.012) | −0.25 (p 0.06 / 0.12) | 0.23 (p 0.009 / 0.013) |
| B1 | 0.47 (<0.001) | 0.45 (<0.001) | 0.27 (p 0.04 / 0.16) | 0.57 (<0.001) |
| B2 | 0.32 (<0.001) | 0.37 (<0.001) | 0.19 (p 0.13 / 0.31) | 0.39 (<0.001) |
| B3 | 0.06 (p 0.37 / 0.61) | 0.12 (p <0.001 / 0.014) | −0.19 (p 0.11 / 0.17) | 0.18 (p 0.016 / 0.034) |
| red_noise | 0.41 (<0.001) | 0.41 (<0.001) | 0.23 (p 0.07 / 0.22) | 0.49 (<0.001) |
| flicker_rms | 0.51 (<0.001) | 0.48 (<0.001) | 0.32 (p 0.04 / 0.13) | 0.60 (<0.001) |

**Reading the ratios**
- The low-frequency (≲ 3 c/d) stochastic part grows faster than the high-frequency floor as
  the disk grows.
- The strongest coherent peak (A1) and the 2–4 c/d band (B3, where the g-mode groups sit) do
  not, in the main sample or in SPOC stars.
- In SPOC stars alone (329 stars) the ratio tests are positive in sign but no longer survive
  the cross-star null.

## Classifier (stage 5): `classifier_summary.txt`, `cluster_vs_phase.csv`, `figures/full_clusters.png`

- Method: PCA to 10 components, then a Gaussian mixture. BIC picks **K = 6**, the top of the
  3–6 search range (BIC still falling), so the true optimum may be larger.
- Cluster sizes: 1,646 / 1,603 / 1,384 / 1,506 / 1,358 / 2,160.
- The clusters split by red-noise slope, by frequency groups (near 1.9/3.7 c/d, near 2.5/5
  c/d, a hump near 0.8 c/d), and by signal-to-noise (a flat, low-SNR cluster, median 1 peak).
- 1,094 of 2,449 multi-sector stars change cluster between sectors.

**Link to disk phase**
- The cluster × phase χ² gives a within-star permutation p of 0.0003.
- The pair statistic, P(cluster changes | phase changes) − P(changes | same phase), is
  +0.028:
  - within-star permutation p = 0.0005;
  - **cross-star p = 0.73**, with a null mean of +0.033 (above the observed value).
- So other stars' disk curves at the same sector times produce the same association. It
  reflects time structure, such as QLP cadence eras, not disk phase.
- **No evidence that spectral shape follows disk phase. A CNN is not justified on this basis.**

## Pilot vs full

| | Pilot (main) | Full (main) |
|---|---|---|
| Stars / sectors | 40 / 151 | 1,084 / 4,108 |
| A1 β | 1.11 [0.71, 1.40] | 0.40 [0.31, 0.50] |
| B2 β | 1.00 | 0.66 |
| white_noise β | 0.36 (n.s.) | 0.33 (p < 0.001) |

- The pilot slope was a small-sample overestimate. The full slope lies below the pilot's 68%
  interval.
- The pilot's flat control was only a lack of power.
- The pilot's SPOC sign flip (11 stars) does not survive: SPOC stars are positive for B1,
  red_noise and flicker, and null for A1.
- The pilot classifier conclusion (no robust link) stands.

## Resources

- **Cache:** 4.5 GB in `~/astro-cache/be-star-pulsations`. TESS light curves are 4.3 GB,
  China-VO downloads 109 MB, spectra 40 MB and WISE 30 MB. None of it is committed.
- **Wall time:**

  | Stage | Time |
  |---|---|
  | 1: 6 processes | 104 min, plus a 70 s resume pass |
  | IRSA (in parallel with stage 1) | about 100 min (32 × 100 stars) |
  | 1b: crowding proxy | about 9 min |
  | 2 | about 3 min |
  | 3 | 4 min |
  | 4 | 12 min |
  | 4c | 3 min |
  | 5 and figures | a few minutes |

  An earlier ThreadPool attempt at stage 1 hung after 25 stars (astroquery downloads are not
  thread-safe) and was replaced by processes.

## Deviations and caveats

- Everything listed under deviations in the pilot write-up still applies:
  - spline detrending;
  - the restated synthetic check 1;
  - the crowding proxy;
  - no AllWISE query;
  - no stage 2.5 (the catalogue has no outburst dates).
- New in the full run:
  - the spline-overshoot guard and w1_quiet clip;
  - stage 4c (post hoc);
  - figure titles carry the run name.
- QLP dominates: 755 of the 1,084 main-sample stars. The SPOC/TESS-SPOC subset (329 stars)
  is the cleaner photometry and gives the weaker result for A1.
- Disk excess is from W1 alone. The dilution correction assumes a fixed optical/IR flux ratio
  (k), which is not measured per star.
- p-values are not corrected for the 7 metrics × 6 variants × 3 k grid. The main conclusions
  rest on p at the 10⁻⁴ floor, which survives any reasonable correction.
