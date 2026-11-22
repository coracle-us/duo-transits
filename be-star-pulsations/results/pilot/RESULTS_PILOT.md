# Be star pulsations vs WISE disk phase: pilot results (2026-11-22)

Brief: `briefs/2026-11-04-be-star-pulsations.md` (Claude relay). Catalogue: Jia et al. 2026,
arXiv 2611.02450, China-VO r101801. Every number below comes from a file in
`results/pilot/` written by a script in `src/`.

## Go / no-go

| Criterion | Value | Verdict |
|---|---|---|
| ≥ 30/100 stars with ≥ 2 usable sectors whose disk_excess differs by ≥ 0.1 mag | 51 using all usable sectors; **19** with the crowding cut that defines the main sample | pass without the cut, **fail** with it |
| Noise floor: scatter of Δln(A1) between same-star sectors with \|Δdisk_excess\| < 0.03 mag below 0.3 | rms **0.56** (122 pairs, main); 0.49 (247 pairs, no cut). Robust 1.4826×MAD: 0.35 / 0.33 | **fail** |
| Synthetic test | ALL CHECKS PASSED (`tests/test_synthetic.py`) | pass |
| Projected full run on the Mac | ~3.6 GB cache, ~11 h wall (about 6 h if stages 1 and 2 run in parallel), mostly MAST and IRSA network time | pass |

**Verdict: no-go as the brief defines it.** Two criteria fail:
- **Noise floor.** Be star amplitudes wander by about 40–55% between sectors even when the
  disk does not change. This is intrinsic variability rather than measurement noise: the
  white-noise control is flat.
- **Crowding.** The main-sample cut keeps only 40 stars with ≥ 2 sectors.

The pilot signal is nonetheless strong (below). Jon decides whether to relax the two criteria.

## Catalogue (stage 0)

China-VO r101801 needs no login. It holds five files:
- `bet_table.csv` (408 kB): 3,300 rows, all with a TIC.
- `wise_lc.zip` (4.2 MB): 3,300 per-visit W1/W2 CSVs with columns mjd, w1mpro_ep,
  w1sigmpro_ep, w2mpro_ep, w2sigmpro_ep, w1_w2, err_w1_w2.
- `bn_tess_data.zip` (20 MB): per-sector TESS CSVs for 285 stars.
- `bn_tess_lc.zip` (143 MB) and `wise_lc_fig.zip` (58 MB): not downloaded.

`bet_table.csv` columns: DR3Name, name (Gaia DR3 source_id), TIC, ra, dec, color_ptp,
w1_ptp, r_value (Pearson r), Gmag, BPmag, RPmag, P, Type, best_author.
- Type: blank for 3,163 rows, 0 for 74 (checked, not EB), 1/2/3 for 33/17/13 EBs.
- **No outburst/decay dates or durations**, so the stage 2.5 cross-check was not possible.
- No licence is stated, so `results/catalogue.parquet` is gitignored and only the loader is
  committed.

The China-VO WISE curves mostly stop in 2023 (median last visit MJD 60196). They were
extended with IRSA NEOWISE-R visits: 98 of 100 stars gained visits, median 1 extra. Where
both sources have the same visit, they agree to 0.000 mag (median).

## Sample (stage 1)

- 140 draws (seed 20261104) from the 3,237 non-EB stars gave 100 passes.
- Failures: 29 had fewer than 2 sectors before the last NEOWISE epoch (MJD 60519) in the four
  products; 11 had no TESS light curves at all.
- Preferred product: QLP 87, TESS-SPOC 7, SPOC 120 s 5, ELEANOR-LITE 1.
- 389 sectors in total; 26 are flagged mixed_pipeline.
- QLP flux column: SAP_FLUX. The QLP header lists SAP_FLUX, KSPSAP_FLUX, KSPSAP_FLUX_SML,
  KSPSAP_FLUX_LAG, QUALITY, ORBITID, …
- ELEANOR-LITE flux column: CORR_FLUX.
- **Crowding.**
  - TIC contratio is empty for nearly all QLP stars (TIC only computes it for CTL stars).
  - `src/s1b_crowding.py` adds a proxy instead: flux of TIC neighbours within 63″, weighted by
    a Gaussian PSF with σ = 16″, divided by the target flux. Rows using it are flagged
    `crowd_source = tic_proxy` (332 sectors).
  - Only 129 of 330 usable QLP sectors pass contratio ≤ 0.1, against 24 of 27 SPOC/TESS-SPOC
    sectors that pass CROWDSAP ≥ 0.9.

## Disk phase (stage 2)

- Spline: FITPACK, s = 0.04, kept.
  - s = 0.04·N/25 makes no visible difference (`figures/spline_check.png`).
  - s = 0.04 underfits low-amplitude stars, e.g. a 0.1 mag bump in TIC 407082659.
  - The spline swings across the 2011–13 WISE gap, so `w1_quiet` is taken only within 100 d of
    a real visit.
- WISE-usable sectors: 361 of 389. Of the rest, 10 fall outside the WISE baseline and 18 have
  wise_gap > 250 d.
- Phase labels, usable sectors (all / crowd-cut):

  | Label | All | Crowd-cut |
  |---|---|---|
  | plateau | 113 | 44 |
  | quiescent | 110 | 57 |
  | building | 73 | 26 |
  | dissipating | 63 | 26 |

- disk_excess percentiles (5/25/50/75/95) over usable sectors: 0.00 / 0.03 / 0.10 / 0.31 /
  0.69 mag.
- 5 stars (19 sectors) have median W1 < 8.

## Metrics (stage 3)

- 387 of 389 sectors have ≥ 400 30-min bins.
- **Usable sectors** (WISE-usable and metrics OK): 359.
- **Main sample:** usable sectors that also pass the crowding cut. That is 153 sectors, of
  which 151 come from the 40 stars with ≥ 2 such sectors.
- Median A1 is 4–5 ppt. f1 falls in 0.25–0.8 / 0.8–2 / 2–4 / 4–12 c/d for 88 / 148 / 126 / 20
  sectors.
- **Deviation from the brief:** the 5-day running median was replaced by the brief's
  alternative, a cubic spline with 5-day knots. The running median's response ripples by
  ±30% at 0.25–0.5 c/d and inflates a 1.1 c/d sinusoid by 10%. The spline's transfer is 0.88
  at 0.25 c/d, 0.97 at 0.5 and 0.99 at 1.1.

## The test (stage 4)

- β is the change in ln(metric) per mag of disk excess, from a within-star demeaned pooled
  least-squares fit.
- The interval is a 68% bootstrap over stars.
- p-values: 10,000 draws each for the within-star permutation null (p_perm) and the
  cross-star disk-swap null (p_swap).
- Full table: `test_results.csv` (7 metrics × 6 variants × k ∈ {0, 0.25, 0.5}).

Main sample (40 stars, 151 sectors):

| metric | β (k=0) | 68% | p_perm | p_swap | β (k=0.5) | p_perm | p_swap |
|---|---|---|---|---|---|---|---|
| A1 | 1.11 | 0.71–1.40 | 0.024 | 0.001 | 1.57 | 0.002 | <0.001 |
| B1 | 1.21 | 0.70–1.76 | 0.034 | 0.003 | 1.67 | 0.007 | <0.001 |
| B2 | 1.00 | 0.71–1.44 | 0.003 | 0.001 | 1.46 | <0.001 | <0.001 |
| B3 | 0.63 | 0.05–1.10 | 0.30 | 0.023 | 1.09 | 0.085 | 0.001 |
| red_noise | 0.83 | 0.44–1.27 | 0.028 | 0.008 | 1.29 | 0.002 | 0.001 |
| flicker_rms | 0.91 | 0.49–1.57 | 0.016 | 0.031 | 1.37 | 0.002 | 0.002 |
| white_noise (control) | 0.36 | −0.38–1.03 | 0.49 | 0.13 | | | |

No crowding cut (93 stars, 357 sectors):

| metric | β (k=0) | 68% interval | p_perm | p_swap |
|---|---|---|---|---|
| A1 | 0.94 | 0.64–1.19 | <0.001 | <0.001 |
| B1 | 1.00 | | <0.001 | 0.001 |
| B2 | 0.94 | | <0.001 | <0.001 |
| B3 | 0.56 | | 0.025 | 0.005 |
| red_noise | 0.77 | | 0.001 | 0.001 |
| flicker_rms | 1.04 | | <0.001 | 0.003 |
| white_noise | −0.01 | | 0.95 | 0.38 |

**Reading the table**
- Amplitudes rise with the disk. β ≈ 1 per mag means about +10% amplitude per 0.1 mag of W1
  excess.
- Dilution correction (k > 0) strengthens the effect rather than creating it, so it is not
  dilution.
- The cross-star nulls centre on 0 (swap_null_mean −0.13 to +0.17 over all rows), so sector-era
  systematics do not produce it.
- The data-quality control (white_noise) is flat at k = 0.
- The sign is opposite to the anti-correlation Niwano et al. 2024 found in Be X-ray binaries,
  which used a different metric.

**Phase contrast** (building minus dissipating, demeaned ln m) is positive for most metrics.
It is significant only without the crowding cut: A1 +0.23, p 0.001 / 0.001.

**Leverage** (`jackknife.csv`, from `src/s4b_jackknife.py`). Dropping any one star keeps β
for A1 between 0.92 and 1.31 in the main sample. The median per-star slope is +1.15, and 64%
of stars have positive slopes. The within-star null is bimodal (`figures/pilot_pooled.png`)
because two-sector stars flip sign under permutation, but no single star carries the result.

**Robustness rows that flip the sign:** only "SPOC/TESS-SPOC stars only" (11 stars, 48
sectors).
- B1 −0.81, red_noise −0.36, flicker_rms −0.56 and B3 −0.58. None of these is significant
  (p ≥ 0.2).
- A1 and B2 stay positive but are not significant.

So the positive signal comes from QLP stars (29 stars: A1 β 1.44, p 0.013 / 0.001). Whether
it holds in SPOC data needs the full run.

**Other robustness rows** (drop mixed_pipeline, drop W1 < 8) keep the sign and stay
significant for A1, B1, B2, red_noise and flicker_rms.

**Noise floor** (rms Δln m over same-star pairs with |Δx| < 0.03, main / no cut):

| Metric | Main | No cut |
|---|---|---|
| A1 | 0.56 | 0.49 |
| B1 | 0.72 | 0.63 |
| B2 | 0.56 | 0.53 |
| B3 | 0.47 | 0.46 |
| red_noise | 0.67 | 0.60 |
| flicker_rms | 0.94 | 0.82 |

## Classifier (stage 5)

- Method: PCA to 10 components, then a Gaussian mixture. BIC picks K = 4 (14,670, against
  14,676 for K = 3). Cluster sizes: 45 / 87 / 124 / 103.
- The clusters are, roughly:
  - red-noise dominated (0);
  - flat, noise-limited (1, low SNR, median 2 peaks);
  - two frequency groups near 1.8 and 3.6 c/d (2);
  - one group (3).
- They partly sort by signal-to-noise rather than by shape (`figures/pilot_clusters.png`).
- Within-star changes: 26 of 93 stars change cluster between sectors.
- Link to phase label:
  - The cluster × phase χ² gives a within-star permutation p of 0.003.
  - The direct within-star statistic, P(cluster changes | phase changes) − P(changes | same
    phase), is +0.04, with p = 0.36 (permutation) and 0.43 (cross-star).
- **No robust evidence that spectral shape follows disk phase. Not ready to justify a CNN.**

## Resources

- **Pilot disk use:** 261 MB in `~/astro-cache/be-star-pulsations` (TESS 148 MB, China-VO
  downloads 111 MB, WISE 1.2 MB, spectra 1.5 MB), plus a 500 MB venv. 3.0 TiB free.
- **Pilot wall time:**

  | Stage | Time |
  |---|---|
  | 1: MAST search and download | 14 min |
  | 2: IRSA, 4 upload-joins of 25 stars | 13 min |
  | 3: metrics | 15 s |
  | 4: 126 test runs | 3 min |
  | 5: classifier | 1 min |

- **Full run** (3,237 non-EB stars, about 2,300 passing at the pilot rate):

  | Stage | Projection |
  |---|---|
  | 1 | about 5 h and 3.4 GB |
  | 2 | about 5 h |
  | 3 | 10 min |
  | 4 | 1–1.5 h, after caching the swap table per variant (needed: it is N² DiskCurve evaluations) |
  | Total | about 3.6 GB, about 11 h, about 6 h with stages 1 and 2 in parallel |

## What did not run, or deviates from the brief

- Stage 2.5 (compare with catalogue outburst/decay dates): the catalogue has none.
- The AllWISE single-exposure table was not queried. The China-VO curves already include the
  2010 AllWISE visits, and the brief's qual_frame cut does not exist in `allwise_p3as_mep`.
- Crowding for QLP uses a TIC-neighbour proxy, because TIC contratio is empty.
- Detrending uses the spline, not the running median (see stage 3).
- Synthetic check 1 is restated. A per-sector 5% check on A1 fails about 20% of the time
  through the brief's own 2 ppt red noise. The test now checks that the pipeline's error
  against a raw sine fit is < 5% (worst 2.3%) and that the 20-sector mean is within 5% of the
  truth (0.997). 80% of single sectors fall within 5%.
- Synthetic case C: the within-star null fires in 85% of repeats; the cross-star null fires
  in 0%. This shows why the cross-star null matters.
- A white_noise control metric was added to stage 4.
