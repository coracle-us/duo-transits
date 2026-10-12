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
