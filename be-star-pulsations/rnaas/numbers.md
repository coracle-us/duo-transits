# Be-star RNAAS draft: provenance of every number

Draft: `rnaas/ms.tex`, v2 of 2026-11-24, revised after `rnaas/REVIEW_mythos.md`. Not submitted.

Word count, from `singles/rnaas/aas_wordcount.py`:
- abstract 146 words (limit 150);
- AAS-counted text 1,264 words (limit 1,350).

That counter ran about 2% low on AAS82207R1, so the AAS count is about 1,295. Updated after the second review (`REVIEW_mythos_v2.md`).

All paths below are under `be-star-pulsations/`. Scripts run with `BE_RUN=full`. Abbreviations:
"stage-4 test" is `stats4.run_test`, and "the c5 log" is the printed output of
`src/c5_review_checks.py`.

## Sample and method

| Number | Source |
|---|---|
| 3,300 candidates; Gaia OBA selection; W1 and W1−W2 changing together; mostly unconfirmed | Jia et al. 2026 (arXiv:2611.02450); `REVIEW_mythos.md` M1 |
| 3,237 non-EB stars; 2,541 pass | `results/full/sample.csv` (`src/s1_sample.py`) |
| Products: 2,020 / 219 / 198 / 104 stars | `results/full/sample.csv` |
| 30-min bins, 5-d knots, 0.25–24 c/d | `src/metrics.py` (`prepare`, `BANDS`, `sector_metrics`) |
| Metric definitions | `src/metrics.py` |
| Flux columns: PDCSAP (SPOC), SAP_FLUX (QLP), CORR_FLUX (ELEANOR) | `src/s1_sample.py` `FLUX_COL` |
| Spline s = 0.04 | `src/wise.py` (`fit_spline`) |
| Gap ≤ 250 d; overshoot ≤ 0.1 mag; W1_quiet | `src/wise.py` (`MAX_GAP`, `MAX_OVERSHOOT`, `DiskCurve`) |
| TESS and WISE about 3 months apart | Observing geometry: TESS near opposition, WISE at 90° elongation (`REVIEW_mythos.md` S5) |
| CROWDSAP ≥ 0.9 / contamination ratio ≤ 0.1 | `src/s3_metrics.py` (`CROWD_MIN`, `CONTRATIO_MAX`) |
| 1,084 stars, 4,108 sectors, 755 QLP | `results/full/test_results.csv`, main rows |
| 10,000 draws per null; bootstrap over stars | `src/stats4.py`, `src/s4_test.py` |
| Injection tests | `tests/test_synthetic.py` |

## Results

### Main test (`results/full/test_results.csv`, main, k = 0)

| Number | Value |
|---|---|
| B1 | 0.80 [0.70, 0.91] |
| red noise | 0.74 |
| flicker | 0.84 [0.73, 0.96] |
| A1 | 0.40 |
| white noise | 0.33 |
| ρ | 0.09–0.20 |
| Noise floor (A1) | 0.54 rms |

The abstract's "about 8%" is B1: e^0.080 ≈ 1.083. The "about 3%" is the 5–24 c/d β of 0.31–0.37 in
`freq_beta.csv`.

### Robustness

| Number | Source |
|---|---|
| No crowding cut: 2,449 stars, B1 0.59 | `test_results.csv`, variant `no crowding cut` |
| Drop mixed-pipeline sectors; drop W1 < 8 | `test_results.csv` variants |
| Drop the top 100 leverage stars: B1 1.03 | `results/checks/review_checks.csv`, M5 rows |
| Linear interpolation | `results/checks/linear_interp.csv` |
| Sector fixed effects: every β changes by < 0.03 (largest +0.025, B1 and red noise) | `results/checks/review_extras.csv`, "two-way star+sector FE" |
| Gaia ESP-ELS beStar: 256 stars, B1 0.58 [0.29, 0.87], p_swap 5e-4 | `review_checks.csv`, M1 rows |

Gaia labels come from `classlabel_espels` in the cached `gaia_dr3_params.parquet`. Of the main
stars, 265 are labelled beStar and 256 have at least 2 usable sectors.

### Frequency bins

| Number | Source |
|---|---|
| Absolute: 0.71–0.92 below 1 c/d; 0.31–0.37 above 5 c/d | `results/full/freq_beta.csv` (`src/s6_freq_beta.py`), main variant |
| Relative, recentred bins: 0.74–1.17 below 0.9 f1; 0.43 [0.34, 0.53] at 0.891–1.122 f1; 0.28–0.49 above | `review_checks.csv`, "S1 relative bin main" rows |

### A1

| Number | Source |
|---|---|
| Stable f1: 510 stars, A1 0.23 [0.10, 0.36] | `review_checks.csv`, M2 rows |
| Stable f1 with f1 > 0.8 c/d: A1 0.08, p_perm 0.37, p_swap 0.29 | `review_checks.csv`, M2 rows |

### White noise

| Number | Source |
|---|---|
| Noisiest tercile: white −0.002 (p_perm 0.98, p_swap 0.87); B1 0.63 | `review_checks.csv`, M4 rows |
| Within-star correlation 0.47 | `results/checks/review_extras.csv` |
| Excess coefficient −0.033 [−0.123, 0.058] with ln B1 in the fit | `review_extras.csv` |
| SPOC formal flux error 0.09, p 0.56 | `results/full/control_tests.csv` |
| Ratios 0.41–0.51; SPOC p_swap ≥ 0.13 | `control_tests.csv` |

### Dilution

| Number | Source |
|---|---|
| QLP adds the TIC-expected flux to the difference flux | Huang et al. 2020, RNAAS 4, 204, §II item 7 and §III (arXiv HTML read 2026-11-24) |
| QLP SAP_FLUX median ≈ 1 per sector | Checked in the cached files |
| k = 0.21 (0.11–0.30), 125 stars | `results/checks/dilution_k.csv` |
| SPOC with k = 0.213: A1 0.369, B1 0.881, flicker 0.936; p_swap ≤ 0.0042 | `review_checks.csv`, M3 rows |
| QLP: 0.52, 0.86, 0.90 | `test_results.csv`, QLP variant |
| 4 shell candidates; changes < 0.01 | `results/checks/dilution_k.csv`, `shell_candidates.csv` |

### Scatter and the discussion

| Number | Source |
|---|---|
| Nearly half (519 of 1,084 = 48%) change by < 0.1 mag | c5 log; `results/checks/review_signs.csv` |
| Of the 177 stars with a range > 0.4 mag, 71% (0.712) have a positive B1 slope | `review_signs.csv` |
| 48% of stars have a negative A1 slope | `src/s7b_bexrb_slopes.py` output; `results/full/LITERATURE_CHECK.md` |
| RX J0440.9+4431 = TIC 410336237, slope −0.86 | `results/full/bexrb_slopes.csv` |

## Figure

The figure is `figures/rnaas_fig1.pdf`, copied to `rnaas/` and unchanged in v2. It is made by
`src/c4_rnaas_fig.py`.

## References and DOIs checked (2026-11-23/24)

**arXiv API**
- Jia et al. 2026 (5 authors);
- Niwano et al. 2024;
- Labadie-Bartz et al. 2025;
- Huang et al. 2020.

**Crossref**
- Niwano, MNRAS 534, 2168;
- Labadie-Bartz et al. 2022, AJ 163, 226;
- Mainzer et al. 2014, ApJ 792, 30;
- Wright et al. 2010, AJ 140, 1868;
- Huat et al. 2009, A&A 506, 95;
- Semaan et al. 2018, A&A 613, A70;
- Creevey et al. 2023, A&A 674, A26 (Apsis I).

**DataCite**
- GSFC-ELEANOR-LITE: 10.17909/j2yt-t417;
- NEOWISE-R single-exposure source table: 10.26131/IRSA144. That is the table queried in
  `src/s2_wise.py` (`neowiser_p1bs_psd`).

**UAT keyword IDs** (checked 2026-11-24 against `UAT.json` in github.com/astrothesaurus/UAT,
release 6.0.0 of 2025-12-17): Be stars 142, Circumstellar disks 235, Stellar pulsations 1625,
Time domain astronomy 2109. All four match their names and none has a change note.
