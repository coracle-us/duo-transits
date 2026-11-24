# Blind re-derivation report (2026-11-24)

**Provenance.** A separate Claude agent wrote this report and the code in this folder. It was
given only `RULES.md`, the raw TESS FITS files, a fresh venv and arXiv:2611.17590, in a folder
outside the repo, and was told not to read the main analysis. Its harness would not let it write
this file, so the main session saved it from the agent's final message. The text below is that
message with headings tidied; the numbers are unchanged.

The working-folder paths in the original message pointed at a scratch folder. The code and
outputs were copied here. Not copied: the `out/*.npz` and `out/*.pkl` intermediates, which hold
light-curve arrays.

## Result

No LTTE signal: M0 (linear) is preferred. ΔBIC is +1.66 for M1 and +5.97 for M2. The median 95%
minimum tertiary mass is 0.199 M☉ (P3 ≥ 12 d, excluding 330–400 d).

## Eclipses used (main products, PDCSAP): 53

| Sector | Primary | Secondary | Dropped by the coverage cut |
|---|---|---|---|
| S18 TESS-SPOC 1800 s | 8 | 9 | P E=2960 (30% window, no ingress) |
| S58 SPOC 120 s | 10 | 12 | P E=3423 (69% window, 6% ingress) |
| S85 TESS-SPOC 200 s | 7 | 7 | S E=3729.5 (58% window) |

- **Errors per eclipse:** typically 4–7 s in S58 and S85, and 15–90 s in S18.
- **S18 light curve:** the TESS-SPOC PDCSAP curve has an out-of-eclipse rms of 4.5 ppt after
  removing pulsations, against 0.83–0.89 ppt in the other sectors.
- **Error scale factors** (√χ²_ν about a linear fit, applied only when > 1):

  | Sector | Primary | Secondary |
  |---|---|---|
  | S18 | 1.557 | 1.002 |
  | S58 | 1.038 | 1.286 |
  | S85 | 1.000 | 1.000 |

## Models (N = 53, primary and secondary together, free secondary offset)

| Model | χ² | dof | k | BIC | ΔBIC |
|---|---|---|---|---|---|
| M0 linear | 55.649 | 50 | 3 | 67.560 | 0 |
| M1 quadratic | 53.341 | 49 | 4 | 69.222 | +1.662 |
| M2 linear + circular LTTE | 49.706 | 47 | 6 | 73.527 | +5.967 |

- **M0:** χ²_ν = 1.11, p = 0.27.
- **M1:** dP/dE = (9.8 ± 6.5)×10⁻¹⁰ d/cycle.
- **M2:** best P3 10.23 d, A = 3.0 s. Other near-equal minima lie at 11–13 d and at 63 and 77 d.
  The grid had 20,000 P3 values, uniform in frequency. No eccentric fit was made.

## Other numbers

- **Period:**
  - P = 2.372222090 ± 0.000000068 d, and P − P_Lee = +0.018 ± 0.006 s.
  - T0 at cycle 3420 = 2882.5620367 ± 0.0000141 BTJD.
- **Secondary offset from phase 0.5:** +6.84 ± 1.95 s.
  - Per sector: S18 +10.7 ± 18.5, S58 +5.8 ± 2.5, S85 +8.0 ± 3.2 s.
  - The predicted Rømer delay is 9.36 s.
- **Primary-only P:** 2.372222072 ± 0.000000087 d (χ² 21.5 for 23 dof).
- **Secondary-only P:** 2.372222120 ± 0.000000111 d (χ² 34.0 for 26 dof).
- **Sector-mean O−C against Lee Eq. 1 (s):**

  | Sector | Primary | Secondary |
  |---|---|---|
  | S18 | +21.3 ± 12.5 | +32.0 ± 13.7 |
  | S58 | +18.7 ± 1.3 | +24.5 ± 2.1 |
  | S85 | +24.5 ± 2.0 | +32.5 ± 2.5 |

## LTTE limit

**Setup:** detection threshold χ²₀.₉₉₇₃(50) = 82.32 (observed 55.65); 1,000 phases per P3;
600 log-spaced P3 values.

**Median m3,min** (P3 ≥ 12 d, excluding 330–400 d):

| Range | Median m3,min | Median A95 |
|---|---|---|
| All | 0.199 M☉ | 28.9 s |
| 12–100 d | 0.225 M☉ | 9.7 s |
| 100–1,000 d | 0.199 M☉ | 40.5 s |
| 1,000–3,000 d | 0.178 M☉ | 108.9 s |

**Extremes:**
- Maximum excluding 330–400 d: 0.873 M☉ at 182.5 d.
- Maximum including 330–400 d: 1.96 M☉ at 362.3 d.
- Minimum: 0.124 M☉ at 129.6 d.

## Choices where RULES.md was open (summary)

1. **Template:** a weighted cubic UnivariateSpline through 0.001-phase bin medians.
2. **Eclipse centres:** the symmetry centre from minimising ∫[T(c+x) − T(c−x)]². **The fitted
   model is the symmetrised shape.**
3. **S58 iteration:** template, then pulsations, then template, with 5 d and then 1.5 d
   detrending splines.
4. **Other sectors:** the model is trend × (1 + b(T_int − 1) + Puls), with a dilution b per
   sector. Pulsation columns are multiplied by sinc(f·texp), and the eclipse model is
   sub-sampled.
5. **Pulsations:** 14 fixed frequencies plus a constant; unweighted least squares out of eclipse.
6. **Detrending:** split at gaps longer than 0.5 d.
7. **Eclipse fit:** c0 + c1·dt + s·(T(dt − shift) − 1), with `least_squares(lm)`.
8. **Bootstrap:** moving blocks of 0.05 d, 500 resamples.
9. **M2:** a sin + b cos, with k = 6.
10. **Injection–recovery:** the exact threshold amplitude per phase; A95 is the 950th of 1,000.
11. **Per-sector offset:** P fixed at the M0 value.

## Checks

- **Synthetic injections:** recovered without bias at 200 s and 1800 s. A +30 s timestamp
  shift moves every time by +30.00 s.
- **Bootstrap caveat:** on white-noise synthetic eclipses, the 0.05 d block bootstrap
  underestimates the scatter by 1.3–1.8×. On the real data, χ²_ν per sector is 0.63–1.65, so
  error scaling covers it.
- **QLP check (QLP − main):**

  | Sector | Primary | Secondary |
  |---|---|---|
  | S18 | about −14 ± 13 s | — |
  | S58 | +1.6 ± 2.0 s | +2.8 ± 3.0 s |
  | S85 | −26 ± 11 s | −46.6 ± 9.8 s (4.8σ, systematic) |

  - Cause: S85 QLP TIME minus TESS-SPOC TIME is +1.94 s until BTJD 3623.45, then −103 to
    −122 s.
  - QLP is a constant +1.94 s later in S18.
- **S18-from-QLP variant:**
  - M0 χ² 46.01 for 49 dof; ΔBIC +3.95 (M1) and +4.60 (M2).
  - P = 2.37222212(4) d; secondary offset +7.15 ± 1.74 s.
  - Median m3,min 0.169 M☉, and 0.198 / 0.167 / 0.102 M☉ by range.

## Not done

- Pulsation delays (RULES 4).
- The S58 SAP check.
- An eccentric LTTE fit.
- In-eclipse pulsation amplitude modelling.
