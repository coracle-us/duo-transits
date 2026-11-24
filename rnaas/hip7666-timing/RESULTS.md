# HIP 7666 timing: results (2026-11-24)

**Inputs:** rules in `RULES.md`, Step 0 in `STEP0_PAPER.md`, venv from `requirements.txt`.

**Reproduce:**
```
.venv/bin/python src/fetch.py
.venv/bin/python src/s1_eclipses.py main 500
.venv/bin/python src/s1_eclipses.py qlp 500
.venv/bin/python src/s1_eclipses.py sap 500
.venv/bin/python src/s3b_variant_q18.py
for k in main qlp sap mixq18; do .venv/bin/python src/s3_models.py $k; done
.venv/bin/python src/s2_pulsation.py
.venv/bin/python src/s4_figure.py
.venv/bin/python tests/test_synthetic.py
```
Light curves (`data/lc/`) and `results/prepared_*.pkl`, which holds light-curve arrays, are not committed; rerun s1 before the test. Every number below is in a file under
`results/`.

## Answer

**No third body is detected.**
- The 53 TESS eclipse times (2019–2024) fit a linear ephemeris: χ² = 56.7 for 50 dof.
- A quadratic term (ΔBIC +1.7) and a circular LTTE orbit (best P3 = 63 d, A = 4.8 s, ΔBIC +6.7)
  are both worse by BIC.
- The pulsation phases cannot serve as a clock in this star (see "Pulsation timing" below).

## Sectors and products

| Sector | Product | Cadence | Eclipses used (p + s) | Median error per eclipse (p / s) |
|---|---|---|---|---|
| S18 (2019) | TESS-SPOC PDCSAP | 1800 s | 8 + 9 | 28 / 45 s |
| S58 (2022) | SPOC PDCSAP | 120 s | 10 + 12 | 3.3 / 5.7 s |
| S85 (2024) | TESS-SPOC PDCSAP | 200 s | 7 + 7 | 4.7 / 6.9 s |

- Three eclipses were dropped for coverage (`results/dropped_eclipses_main.csv`).
- Error scale factors from scatter about a line within each sector: 1.00–1.46 (`results/models_main.txt`).
- S123 is predicted by tess-point but not yet observed.

## Eclipse timing (`results/models_main.txt`, `results/oc_main.csv`)

**Linear ephemeris (M0):**
- T0 = BTJD 3131.645322(11), i.e. BJD_TDB 2460131.645322, at cycle 3525 of Lee Eq. 1.
- P = 2.37222211(6) d.

**Comparisons**
- **Long baseline:** P(TESS) − P(Lee Eq. 1, 2000–2022) = +0.020 ± 0.016 s, so no period change
  over 22 years.
- **Primary only vs secondary only:** P = 2.37222211(7) and 2.37222212(11) d. They agree, so
  there is no apsidal motion.
- **Secondary − primary:** +12.7 ± 1.9 s overall; per sector +14.8 ± 18.4, +12.1 ± 2.4 and
  +13.2 ± 3.1 s.
  - The inner-orbit Rømer delay predicted from Lee's K values is P(K_B − K_A)/(πc) = +9.4 s.
  - The rest, 3.3 ± 1.9 s, corresponds to e cos ω ≈ (2.5 ± 1.5) × 10⁻⁵.
  - Lee et al. saw this offset (0.00019 d) and suggested pulsations as the cause; light travel
    time across the binary explains most of it.

**O−C against Lee Eq. 1 (s, weighted means):**

| Sector | Primary | Secondary |
|---|---|---|
| S18 | +18.0 ± 11.9 | +32.8 ± 14.0 |
| S58 | +15.5 ± 1.2 | +27.5 ± 2.1 |
| S85 | +22.2 ± 1.7 | +35.4 ± 2.6 |

The offset of about +15 s is the difference between Lee's least-squares ephemeris and the TESS
data; their own Fig. 2 shows the TESS points above zero as well.

**Checks**
- **Against Lee Table 1 (S58, 22 matched eclipses):** ours − theirs = −1.4 s on average, with
  5.7 s scatter (I: +0.2 ± 4.5; II: −2.8 ± 6.4).
- **S58 SAP instead of PDCSAP:** O−C +14.9 / +28.0 s against +15.5 / +27.5 s (`results/models_sap.txt`).
- **QLP** (`results/models_qlp.txt`):
  - S18 and S58 agree with the main products.
  - **S85 QLP does not.** For the same cadence numbers, the QLP `TIME` column is 103–122 s
    earlier than TESS-SPOC's from BTJD 3623.448 to the end of the sector. Before that it is
    1.9 s later, as it is in S58.
  - The QLP S85 eclipse times show the same step. It is a timestamp problem in that QLP file,
    not a signal in the star.
  - Without that sector, QLP agrees.

## LTTE limit (`results/ltte_limit_main.csv`, `results/ltte_limit_mixq18.csv`)

The limit is the semi-amplitude detected at 3σ in 95% of orbital phases (injection–recovery
on the real times). It converts to a minimum tertiary mass with M12 = 3.068 M☉ and
sin i3 = 1.

| P3 range | A95 (main) | m3 excluded above (main) | With S18 from QLP |
|---|---|---|---|
| 12–100 d | median 9 s, max 30 s | median 0.19 M☉ (max 0.40) | 0.17 M☉ |
| 100–1,000 d | median 42 s | median 0.20 M☉ | 0.15 M☉ |
| 1,000–1,820 d | median 93 s | median 0.17 M☉ | 0.11 M☉ |
| 1,820–3,000 d | median 124 s | median 0.17 M☉ | 0.08 M☉ |

- **Blind spot:** P3 ≈ 364 d. S18, S58 and S85 fall at the same season, so the limit there is
  only A < 350 s, m3 < 2.1 M☉.
- **Caveats:**
  - For P3 longer than the 1,820 d TESS baseline, LTTE trades off against the linear and
    quadratic terms.
  - Below about 12 d a tertiary would be dynamically unstable, and dynamical (non-LTTE) timing
    effects are not modelled.

## Pulsation timing (`results/pulsation_summary.txt`, `results/pulsation_delays.csv`)

**Data:** f1, f3 and f4 (and a common delay from all 14 modes), one orbit per segment.

**Errors on the common delay:**

| Sector | Error per orbit |
|---|---|
| S58 | 13 s |
| S85 | 17 s |
| S18 | about 290 s (1800 s cadence, f1 above Nyquist) |

**Excess scatter:** after removing each sector's linear term, the common delays scatter more
than their errors, χ² = 34.1 for 10 dof (S58) and 15.4 for 6 (S85). That is rms 22–24 s,
against 3–7 s per eclipse.

**Phase connection S58–S85:** not achieved under the rule. f1 and f4 pick one alias
(Δχ² = 32 and 19) but f3 does not (Δχ² = 7.5). Even a connection would leave a constant delay
between the two sectors absorbed by the frequency.

**Inner-orbit check: failed.**
- Over S58 (0.3 d segments), the delay varies with orbital phase as
  a cos 2πφ + b sin 2πφ, with a = +45 ± 21 s and b = +11 ± 14 s.
- Light travel time of the primary would give a = −10.1 s and b = 0.
- The synthetic test (`tests/test_synthetic.py`) recovers an injected a = −60 s as
  −58.8 ± 6.8 s, so the sign convention and the code are right. The phases genuinely move with
  orbital phase, in a way light travel time does not explain.
- Many of the frequencies are f_i ± n f_orb combinations (Lee Table 5), which suggests tidally
  affected modes.
- So the pulsations are reported, not used to set the limit.

## Prior work checked

- **Lee et al. 2026 (arXiv:2611.17590):** S58 only, linear ephemeris, no third body.
- **Pawar et al. 2024 (arXiv:2410.00590):** S58, truncated to four orbits; no timing.
- **Escolà-Sirisi et al. 2005 (arXiv:astro-ph/0505342):** discovery and minima.
- **Feng et al. 2021 (MNRAS 508, 529):** read only through the abstract, as the full text was
  blocked to the fetch. It gives P = 2.3722200(4) d and no third-body claim.
- **Searches:** arXiv for "HIP 7666" and two web searches found no ETV or LTTE study.
- **Not available:** SIMBAD/CDS was down on 2026-11-24, and ADS was not queried (it needs a
  token).

## Draft

`DRAFT_RNAAS.md`: 1,195 words in total including the title, caption and references, with a
148-word abstract. The count is a plain word count of the Markdown.

Excluding P3 < 12 d and the 330–400 d blind spot, the limit is below 0.3 M☉ at 94% of grid
periods. The median is 0.18 M☉ and the worst case 0.45 M☉.
