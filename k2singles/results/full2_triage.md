# run_full2 triage (2026-10-12)

**Inputs:**
- Results: `results/full2/summary.csv` (144 events, all status ok).
- Status check: `data/lead_status.csv` (`scripts/status_check.py`). It covers bulk TOI/CTOI, ExoFOP
  target JSON, Prša 2022 EBs, Gaia DR3 (RUWE, RV variability, NSS orbit, vari_eclipsing_binary,
  classifier), VSX and SIMBAD.

**Yield by max_support** (the number of repeats explained by the best period):

| max_support | events |
|---|---|
| 0 | 84 |
| 1 | 41 |
| 2 | 16 |
| 3 | 2 |
| 4 | 1 |

Events with any repeat: 64. Events with a K2 repeat: 26.

## Leads with support ≥ 2 (first-pass verdicts, NOT reviewed yet)

| EPIC | TIC | depth, dur | best periods (support) | status check | first-pass verdict |
|---|---|---|---|---|---|
| 211924561 | 238606535 | 2.1 ppt, 3.4 h | 70.247 / 140.495 (3) | TOI-5560.01 | positive control: recovered |
| 211411112 | 437034781 | 26 ppt, 4.4 h | **62.595 (2)**, unique | Sanders 617 = WOCS 11028 (M67), SB* | **validation**: published 62.593 d (K2 M67 Study, arXiv:2011.08365) |
| 211892898 | 184895645 | 9 ppt, 62 h | 72.506 (2), with 15.99/18.13/36.25 | RG*; Gaia DR3 SB1 **P = 72.534 d**, e = 0.48 | red-giant EB; eclipse alias matches the Gaia SB1 period. Nardiello+2016 (arXiv:1608.00459): EB, "P > 50 d". Period possibly new |
| 211390677 | 366568535 | 6.8 ppt, 40 h | 61.442 / 92.163 / 184.327 (2) | RG*; no Gaia NSS | K2 C5 + C18 events clean and identical; S45 QLP marginal. Giant host, so the companion is probably stellar |
| 211503363 | 20444116 | 4.6 ppt, 23.5 h | 61.451 / 122.903 (2) | SIMBAD Er* | both TESS repeats are QLP at SNR ~10, noisy; weak |
| 201132684 | 35021200 | 2.2 ppt, 12.8 h | 24.418 (2) | **K2-158** (confirmed b, c), CTOI | known system, drop |
| 202072917 | 308255000 | 248 ppt, 7.4 h | **34.746 (3)** | EB*; VSX P = 3.842 d | EB. VSX period is wrong: S44, S71 and S72 have no eclipse, S6/S43/S45 one each |
| 202071902 | 57051340 | 91 ppt, 5.4 h | 22.904 (2) | Em*; VSX EA/RS P = 8.697 d | EB. VSX 8.70 d is excluded by S45 and S72 having no eclipse |
| 202126877 | 86491361 | 62 ppt, 5.9 h | 21.775 / 21.952 (2) | EB*; VSX P = 3.877 d | eccentric EB. Dip pairs 7.25 d apart repeat every ~21.78 d (primary + secondary). Pipeline confused; hand-fit needed |
| 202137580 | 238479066 | 67 ppt, 51 h | 34.368 (2) | EB?; VSX EA P = 9.902 d | messy: dips at 2503.38, 2537.61, 3267.26, 3272.46. Unresolved |
| 247967714 | 268396693 | 258 ppt, 4.9 h | 61.875 / 92.813 / 185.626 (2) | nothing | EB, alias set |
| 246696483 | 436612436 | 33 ppt, 17 h | 91.122 / 182.243 (2) | Gaia RV-variable | EB likely, alias pair |
| 248811085 | 349987241 | 91 ppt, 6.4 h | 49.250 (2) | RUWE 1.93 | EB, unique |
| 246089124 | 49603732 | 139 ppt, 9.3 h | 36.981 / 73.961 (2) | nothing | EB, alias pair |
| 211489484 | 6951741 | 109 ppt, 4.4 h | 82.258 (2) | EB* | EB; K2 C18 repeat not on the 82.258 ephemeris |
| 246163364 | 38315113 | 70 ppt, 35 h | 10.69 / 20.23 / ... (2) | Gaia RS CVn | variability-dominated, unreliable |
| 231674968 | 24395505 | 11 ppt, 46 h | 6.88 (4) | RUWE 5.4, 4 Gaia sources within 5″ | variability / blend; reject |
| 212694013 | 178941241 | 7 ppt, 41 h | many (2) | nothing | QLP S91 variability; reject |
| 212026444 | 203254318 | 0.6 ppt, 25 h | many (2) | RUWE 6.5, RV-variable | 10 "repeats" = noise / variability; reject |

**Honest assessment:** apart from TOI-5560.01, no new planet-like lead. The yield is:
- one pipeline validation (Sanders 617, period matches a published orbit to 0.002 d);
- one probably-new red-giant EB period confirmed by Gaia SB1 (211892898);
- one giant-host alias set (211390677);
- about 8 deep-EB periods or alias sets, two of which correct VSX periods.

## Catalogue checks (2026-10-12, after the triage table)

**VizieR all-catalogue cone search, 3″** (`scripts/vizier_census.py` → `data/vizier_census.csv`)
for 247967714, 246696483, 248811085, 246089124, 211489484, 202137580, 211892898, 211390677.

- **No catalogued period** for 247967714, 246696483, 248811085, 246089124 or 211892898.
  - 211892898 appears only in APO-K2 red-giant tables and in the LAMOST DR4 RV-variable list
    (Tian+2020, J/ApJS/249/22: 3 epochs, ΔRV 19 km/s, log g 2.8, Teff 4810). No period.
- **211390677:** Green+2023 ellipsoidal list (J/MNRAS/522/29) gives Porb 5.998 d, but with candidate
  Score 0.0, so it is not a detection. Irrelevant.
- **211489484 = TIC 6951741.** SPOC multi-sector TCE (S1–S46), P = 106.936 d, ExoMiner label FP
  (Valizadegan+2025, J/AJ/170/287).
  - 106.936 = (S46 − S7)/10. It is excluded by the K2 C5 primary (12.31 cycles from S7).
  - 82.258 d = ΔT/13 threads K2 C5 + S7 + S46.
  - **Independent support:** the Villanova K2 EB catalogue (keplerebs.villanova.edu, overview
    ?k=211489484) lists a C5 secondary eclipse: sdepth 0.050, separation 0.1598 × its placeholder
    period of 100 d = 16.0 d after the primary.
  - Our K2 C18 half-depth dip (59 ppt, BTJD 1277.996) sits 0.195 × 82.258 = 16.0 d after a
    predicted primary. The same offset, so the system is eccentric.
- **202137580** is heavily catalogued: VSX 9.90152 d, LAMOST EB 9.9015 d (J/ApJS/280/26),
  LAMOST SB1 69.44 d (J/ApJS/282/55), ASAS-SN 63.36 d, Gaia cross-match 9.90/184.86 d.
  It is messy and already studied. **Dropped from the table.**
- **Villanova K2 EB catalogue:** 202072917, 202071902 and 202126877 are listed with placeholder
  periods (50 d, i.e. unknown; InCat False). 211489484 has a 100 d placeholder. 247967714,
  246696483, 248811085, 246089124, 211892898 and 211390677 are not in it.

## Second Fable review (`review_fable2/REPORT.md`) and follow-up (2026-10-12)

**Verdicts:**

| EPIC | verdict | detail |
|---|---|---|
| 211892898 | confirmed, with caveats | fitted mid-times 189.3319 ± 0.016, 1131.9402 ± 0.017, 1276.9610 ± 0.010; **P = 72.5086 ± 0.0013 d, unique**, K2 alone excludes every sub-multiple (the alias list was pipeline bug P1) |
| 211390677 | confirmed, with caveats | S45 dip real and on target in own TESScut (6.5 ppt, SNR 30); 61.44 d excluded by K2 C5 (edge window, P2); **92.154 or 184.309 d**; possible partial 4th eclipse at the end of S72; red-clump host, so a stellar companion |
| 211411112 | confirmed (validation) | Sandquist+2021 give 62.59484 ± 0.00003 d (their eclipses are secondaries); 62.593 is Geller+2021 RV. Quote 62.59484 |
| 202072917 | confirmed, not novel | 34.74565 d is already a SPOC S1–46 TCE (MES 419; Valizadegan+2025). VSX 3.842 d comes from an Armstrong+2016 ML classification |
| 202071902 | orbit already published | Kovalev+2024 (MNRAS 527, 521) SB2 22.859 ± 0.015 d; our eclipse P 22.90415 d (3σ from theirs, far more precise) |
| 202126877 | reviewer: 21.7753 d resolved | see the follow-up below: **flagged** |
| 246194159 | pipeline miss | five eclipses; two were dropped by the M6 rules (P3) |

**211892898, Gaia test:**
- The Gaia SB1 orbit predicts the primary-eclipse conjunctions at O−C −0.9 to −1.3 d (1.4–1.8σ).
- The secondary-conjunction hypothesis is off by 17 d.
- f(M) = 0.060 M☉. With asteroseismic M = 1.11 M☉ and R = 6.82 R☉, M2 ≈ 0.55 M☉, consistent with the depth and duration.
- Novelty: the orbital period is already known (Gaia DR3 NSS). Nardiello+2016 already reported one primary and one secondary eclipse. New: the eclipse ephemeris, identifying the eclipse with Gaia's primary conjunction, and the companion radius.

**Pipeline problems P1–P4:** fixed in k2lib v3 (see its docstring); singles code untouched.

**v3 control run** (`results/full3/`, 8 events):

| EPIC | v3 result | note |
|---|---|---|
| 211924561 | 70.25 / 140.49 (support 3) | |
| 211411112 | 62.60 (2) | |
| 211892898 | **72.51, unique** | |
| 211390677 | 61.43 / 92.15 / 184.29 | 61.43 needs the edge window, P2 |
| 231674968 | 0 repeats | noise, now clean |
| 212026444 | 1 repeat | noise, now clean |
| 246194159 | 4 repeats recovered | |

**Two-phase ephemerides** (`scripts/two_phase_ephem.py`; every predicted primary and secondary
covered by data must show a dip):
- **246194159** (TIC 301328232): **P = 73.242 d, unique.**
  - Eclipses at 747.39 and 3823.68.
  - Similar-depth eclipses (22–25%) at phase 0.868: 810.85, 3227.97, 3813.93.
  - So e ≳ 0.58.
- **211489484** (TIC 6951741): **P = 82.2555 d, unique**; secondary at phase 0.1952.
  - It agrees with Villanova's C5 secondary 16.0 d after the primary.
  - The SPOC TCE period, 106.94 d, is excluded.
- **202126877** (TIC 86491361): no survivor.
  - The K2 C0 event lies exactly on the TESS primary ephemeris (21.7753 d, cycle 124.000).
  - But it is 4.6% deep, against 17–24% in TESS, and K2 shows NO secondaries at −226.09 / −204.32, where
    TESS-era secondaries are 4–9%.
  - The only Gaia neighbour is G 14.9 at 21″, which cannot make a 20% dip, so the eclipses are on target.
  - Possible eclipse-depth variation (tilted orbit, third body?) or something odd in K2 C0.
  - **Flagged; not claimed.** The TESS-only period 21.7753 d stands.
