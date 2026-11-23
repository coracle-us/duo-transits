# Be star pulsations: pre-submission checks (2026-11-23)

Brief: `2026-11-23-be-star-checks.md`. Scripts: `src/c1_phase_terms.py`, `c2_linear_interp.py`,
`c3_dilution_k.py`, `c4_rnaas_fig.py` (+ `phasefit.py`, `wise.LinearDiskCurve`); run with
`BE_RUN=full`. 10,000 draws per null, seed 20261104; 68% bootstrap intervals; synthetic test passed.

## Check 1. Phase terms under the cross-star swap null (`results/checks/phase_terms.csv`)
Main sample: 1,084 stars, 4,108 sectors (766 building, 843 dissipating).

| Metric | building | p_perm / p_swap | dissipating | p_perm / p_swap |
|---|---|---|---|---|
| A1 | +0.186 | <1e-4 / <1e-4 | −0.086 | 0.004 / 0.012 |
| B1 | +0.250 | <1e-4 / <1e-4 | −0.143 | <1e-4 / 4e-4 |
| B2 | +0.227 | <1e-4 / <1e-4 | −0.156 | <1e-4 / <1e-4 |
| B3 | +0.177 | <1e-4 / <1e-4 | −0.092 | 0.002 / 0.002 |
| red_noise | +0.231 | <1e-4 / <1e-4 | −0.140 | 2e-4 / <1e-4 |
| flicker_rms | +0.231 | <1e-4 / <1e-4 | −0.155 | 3e-4 / 4e-4 |
| white_noise (control) | +0.083 | 0.023 / 0.075 | −0.054 | 0.14 / 0.14 |

- **With the spline, both terms survive the swap null** for all six variability metrics.
- QLP-only (755 stars): the same, every variability term has p_swap ≤ 0.026.
- SPOC/TESS-SPOC-only (329 stars): the dissipating term survives for B1, B2, red_noise and
  flicker_rms (p_swap 0.005–0.016). The building term is weaker (+0.08 to +0.15) and has
  p_swap 0.04–0.22.

## Check 2. Linear interpolation instead of the spline (`results/checks/linear_interp.csv`)
1,086 stars, 4,118 sectors. x_lin vs x_spline: r = 0.993. Stage 4 main test, k = 0:

| Metric | A1 | B1 | B2 | B3 | red_noise | flicker_rms | white_noise |
|---|---|---|---|---|---|---|---|
| β linear | 0.44 | 0.87 | 0.69 | 0.42 | 0.79 | 0.89 | 0.37 |
| β spline | 0.40 | 0.80 | 0.66 | 0.39 | 0.74 | 0.84 | 0.33 |

- Every linear β has p_perm = p_swap < 1e-4. Linear β is 4–11% higher, so **the main result
  does not depend on the spline.**
- **The phase terms do depend on it.** The labels agree for only 69% of sectors; 438
  spline-plateau sectors become building or dissipating.
  - Building gets stronger: +0.21 (A1) to +0.34 (B1), all p < 1e-4.
  - Dissipating shrinks to −0.08 to +0.01. Only flicker_rms keeps p_swap < 0.05
    (−0.082, p_swap 0.038).
  - The white-noise control picks up a building term too: +0.18, p < 1e-4.
- **Verdict on the phase terms:** "building is more variable" is robust in sign and is
  partly shared with the high-frequency noise floor. "Dissipating is less variable" is not
  robust to the disk model, so leave it out of the note or state it as model-dependent.
  That is 42 coefficients across 3 samples, with no multiplicity correction.

## Check 3. Dilution factor k from SAP flux (`dilution_k.csv`, `dilution_k_per_star.csv`, `shell_candidates.csv`)
- **Data:** 125 SPOC120/TESS-SPOC stars and 564 sectors. Only the star's preferred product
  was used, with at least 3 WISE-usable main sectors. Flux is the median SAP × CROWDSAP /
  FLFRCSAP (QUALITY = 0), with ZP 20.44.
- **Pooled within-star dTmag/d(disk_excess) = −0.213** (68% [−0.259, −0.167], 95% [−0.295,
  −0.105]). p_perm = p_swap < 1e-4, swap-null mean ≈ 0.
- **So k = +0.21 (95% 0.11–0.30).** The disk brightens the optical on average, which is the
  sign the dilution correction assumes. Of the stage 4 grid values, k = 0.25 is the
  closest. At that value A1 β = 0.63 and B1 β = 1.03.
- **Sign convention:** stage 4 adds 0.921·k per mag to ln(flux), so the Tmag slope is −k (the
  brief's "−0.92·k" is the ln-flux slope with its sign flipped).
- **Per-star slopes are mostly noise** (median −0.19, IQR [−0.49, +0.15]): SAP scatter at fixed
  disk state is 0.034 mag rms vs a per-star signal of ~0.025 mag. Pooled k is usable but rough.
- **Shell candidates** (slope > 0, p < 0.05): 4. TIC 55295028 is strong (34 sectors, +0.91 ± 0.07,
  p = 2e-14); 336296046, 408991617, 271614171 (3–4 sectors) are marginal. QLP stars untested.
- **Main test without the candidates** (1,080 stars, k = 0): A1 0.400 [0.303, 0.498],
  B1 0.796 [0.697, 0.898], flicker_rms 0.840 [0.730, 0.954]. All have p_perm = p_swap < 1e-4,
  so the result is unchanged.

## Check 4. Figure
`figures/rnaas_fig1.png` and `.pdf`: single column (3.4 × 5.0 in), greyscale. Panels are
B1 and flicker_rms. Each shows sectors by phase (△ building, ○ plateau/quiescent, ▽
dissipating), 10 quantile-bin means ± SE, and the pooled β with its 68% band. β and p
values come from `test_results.csv`. The y-axis is clipped at the 1st–99th percentile
× 1.1, so a few sectors fall outside.

## What went wrong or differs from the brief
- **PR #8 was already merged** (2026-11-22). Branch `be-star-pulsations` was fast-forwarded to
  `main` (it contains c4d3232). The commit went there, and a **new draft PR** was opened
  instead of pushing to #8.
- **Check 2 visit values:** the stage 2 visit tables (China-VO per-epoch W1, IRSA per-visit
  medians), not means recomputed from single exposures.
- **Check 2 sample:** the linear rule (≤ 200 d to each visit) keeps 10 more sectors, 2 more stars.
