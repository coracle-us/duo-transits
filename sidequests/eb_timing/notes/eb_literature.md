# Literature check: eclipse timing variations (ETVs) of TESS eclipsing binaries

Checked 2026-10-11 with web search, arXiv abstract/HTML pages and the arXiv LaTeX sources
(`scripts/known_triples.py` pulls TIC IDs out of 27 of them). There is a systematic arXiv
API listing in `notes/arxiv_search_raw.txt`, run newest first for "eclipse timing" + TESS,
"continuous viewing zone" + eclipsing, "light travel time" + TESS + eclipsing, and
triple + "eclipse timing" + catalog. I did not have ADS API access. Coverage of papers that
are not on arXiv (e.g. some Universe, AJ or ApJS papers without preprints) is incomplete.

## Catalogues and surveys

| Work | Sample | TESS data | What was measured | ETV content |
|---|---|---|---|---|
| **Prša et al. 2022**, ApJS 258, 16 (VizieR J/ApJS/258/16) | 4,584 EBs from 2-min targets | S1–26 | ephemeris (BJD0, P), morphology, polyfit and 2-Gaussian eclipse widths, depths and phases | **None.** The VizieR table has no eclipse times and no ETV or third-body flags (column list checked with astroquery). The live Villanova site (tessebs.villanova.edu) lists sectors up to ~S88, but its only visible flag is "Ambiguous". |
| **Kostov et al. 2025**, ApJS (TESS Ten Thousand, doi 10.3847/1538-4365/ade2d8) | 10,001 FFI EBs (7,936 new), T < 15 | FFI S1–82 | ephemerides, depths, durations, centroid offsets | Free-text comment "pETVs" (potential ETV) per EB. No ETV fits. Volunteers' triple and quadruple candidates are deferred to a later paper. |
| **IJspeert et al. 2021**, A&A 652, A120; **2024**, A&A 691, A242 | ~3,000 OBA EBs (2021); 14,573 OBAF EBs (2024) | first 4 years | eccentricities, circularisation, pulsations | **No ETV analysis.** |
| **Mowlavi et al. 2023**, A&A 674, A16 (Gaia DR3) | 2.2 M EB candidates | Gaia G (sparse) | periods, two-Gaussian models | Sparse Gaia sampling, so no ETVs. Czavalinga et al. 2023 (below) use it. |
| **Czavalinga et al. 2023**, A&A 670, A75 (arXiv 2212.02903) | EBs (>1 M from several surveys) matched with Gaia DR3 NSS orbits much longer than P_EB | some TESS ETVs to validate | triple candidates from Gaia NSS | Gaia-based triple candidates, some checked with TESS ETVs (abstract). Companion paper MNRAS 526, 2830 (four triply eclipsing). |

## Systematic TESS ETV searches

| Work | Sample | Sectors | Result |
|---|---|---|---|
| **Mitnyan et al. 2024**, A&A (arXiv 2402.01486) | >3,500 EBs in and near the **northern** CVZ, own FFI photometry | S14–60 | ETVs for all, 351 non-linear, **135 triple candidates with LTTE/dynamical solutions** |
| **Borkovits, Mitnyan, Czavalinga & Rappaport 2026**, Universe 12, 107 (arXiv 2604.04004) | same NCVZ sample, re-analysed | S14–86 (200-s FFI via lightkurve) | **168 triple candidates, ≥66 certain**; 34 with longer-term variations (≥3 possible quadruples); 2 EBs stopped eclipsing; 7 show depth variations. NCVZ is next observed from S117 (May 2027). |
| **Marcadon & Prša 2024**, ApJ (arXiv 2403.07694) | ~400 Villanova EBs with \|β\| > 78° (both CVZs) → ~100 detached/semi-detached EBs with P < 27 d timed | 2-min S1–60 + eleanor FFI | **10 triple candidates**. LTTE solutions for 5, all northern CVZ: TIC 219900027, 229771234, 259006185, 356896561 (GZ Dra), 424461577. Three of these overlap Mitnyan+2024. |
| **Esmer & Daylan 2025** (arXiv 2503.20900 v2) | 152 detached Villanova EBs (morph < 0.2, P < 14 d), all-sky | 2-min to S80 (first + last sectors; up to 6 primaries per target for the allesfitter fits) | 37 with periodic signals at FAP < 0.1, 19 at FAP < 0.01; best candidate TIC 142979644 (~535 d, ~0.33 min, BD-mass). Occurrence of circumbinary BDs < 40% (2σ). |
| **Borkovits et al. 2016** (Kepler, MNRAS 455, 4136) and **Borkovits et al. 2025**, A&A (arXiv 2502.09480) | ~220 Kepler triples, revisited with TESS | Kepler + TESS 2019–2024 | 193 of the 2016 solutions confirmed or improved; 50 new |
| **Hajdu et al. 2019, 2022** (OGLE-IV bulge; arXiv 1903.00674, 2110.02391) and **Hajdu et al. 2022** CoRoT–TESS (arXiv 2209.01142) | ~80,000 OGLE EBs; CoRoT fields | OGLE; CoRoT + TESS | 992 OGLE LTTE candidates. The CoRoT–TESS paper uses TESS to extend the baseline. Not the TESS CVZ. |

## Compact triples and quadruples found with TESS (individual systems)

These papers come from the Rappaport, Kostov, Borkovits and Powell groups and the Planet
Hunters / Visual Survey Group. They include Borkovits+2020 (TIC 167692429 and 220397947, the
first southern CVZ triples, P_out = 331 d and 77 d, dynamically dominated ETVs), Borkovits+2020
(TIC 209409435), Mitnyan+2020 (TIC 278825952), Borkovits+2022 (northern triply eclipsing),
Rappaport+2023 (9 systems), Rappaport+2024 (7), Czavalinga+2023 (4), Borkovits+2025 (10
ultracompact, "more than a hundred" close triply eclipsing triples known to the group),
Kostov+2024 (TIC 290061484, P_out = 24.5 d), Kostov+2024 (101 quadruple candidates) and
Kostov+2025 (52 more), Powell+2022, 2025, 2026, Zasche+2024 (6 eccentric EBs with a third
body) and Zasche+2026 (6 doubly eclipsing quadruples), and Rappaport+2026 (TIC 366240660,
P_out = 30.5 d). Many of these systems are in the CVZs. Most were found by **eye** (the
triply eclipsing ones) and not from ETV surveys, but each one had a full ETV and photodynamical
analysis.

All TIC IDs named in these papers are in `data/known_multiples.csv` (481 unique TICs). The
harvest is over-inclusive: comparison stars mentioned in the text are included. It also
misses systems the papers list by other names.

## Where the gap is

1. **Northern CVZ: closed.** Two systematic ETV studies from the Borkovits group (S14–60,
   then S14–86) cover >3,500 NCVZ EBs. Everything we could add there is already in their list
   of 168 candidates, or would only be a marginal revision. No new NCVZ data before S117
   (May 2027).
2. **Southern CVZ: open on paper.** I found **no systematic ETV catalogue of the southern
   CVZ EBs** (none on arXiv up to 2026-10-11). Partial coverage exists:
   - Marcadon & Prša 2024 timed CVZ EBs from both poles, but only up to S60, only ~100
     detached/semi-detached systems after their cuts, and they solved only 5 (northern) systems.
   - Esmer & Daylan 2025 timed 152 detached EBs all-sky, with very few eclipses each.
   - The compact triply eclipsing triples in the SCVZ are already solved individually.

   Southern CVZ coverage now runs **S1–13, S27–39, S61–69, S87–90 and S93–98: up to ~43 sectors
   over 2,720 d (7.4 yr)**. That is longer than the NCVZ baseline the Borkovits group analysed
   (~1,950 d). Our top-200 sample (`data/sample_top200.csv`) is 100% southern (81 in the SCVZ
   proper) for exactly this reason. **No one has published a systematic ETV search of the
   southern CVZ EBs with >5 years of data.**
3. **Risk: the gap is a "not yet", not a "nobody cares".** The Borkovits/Mitnyan/Czavalinga
   group has the method, the photometry pipeline and two NCVZ papers. A southern counterpart is
   an obvious next paper for them and may already be in preparation (the 2026 paper does not
   mention one). Anything we do in the SCVZ should be treated as possibly scooped by them.

## Sources

- Prša et al. 2022: https://iopscience.iop.org/article/10.3847/1538-4365/ac324a
- Kostov et al. 2025 (TESS Ten Thousand): https://iopscience.iop.org/article/10.3847/1538-4365/ade2d8
- IJspeert et al. 2024: https://www.aanda.org/10.1051/0004-6361/202450507; 2021: VizieR J/A+A/652/A120
- Mitnyan et al. 2024: https://arxiv.org/abs/2402.01486
- Borkovits et al. 2026: https://arxiv.org/abs/2604.04004 (Universe 12, 107)
- Marcadon & Prša 2024: https://arxiv.org/abs/2403.07694
- Esmer & Daylan 2025: https://arxiv.org/abs/2503.20900
- Borkovits et al. 2025 (Kepler field revisited): https://arxiv.org/abs/2502.09480
- Hajdu et al.: https://arxiv.org/abs/1903.00674, https://arxiv.org/abs/2110.02391, https://arxiv.org/abs/2209.01142
- Czavalinga et al. 2023: https://arxiv.org/abs/2212.02903, https://arxiv.org/abs/2309.01639
- Borkovits et al. 2020: https://arxiv.org/abs/2002.06902
- Rappaport et al. 2024: https://arxiv.org/abs/2403.12041; Borkovits et al. 2025 (ten new): https://arxiv.org/abs/2510.04565
- Kostov et al. quadruples: https://arxiv.org/abs/2309.14200, https://arxiv.org/abs/2510.27470
