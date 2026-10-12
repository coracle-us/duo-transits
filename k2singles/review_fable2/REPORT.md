# Adversarial review 2 of the K2-singles leads (Fable 5.1, 2026-10-12)

Scope: `results/full2_triage.md`, `results/full2/*.json`, `scripts/k2lib.py`, `scripts/run_full2.py`.
All scripts and outputs are in this folder; nothing outside it was modified. Venv python used throughout.

Files: `ecl_211892898.py/.txt/.png` (eclipse timing, alias test, Gaia conjunctions; also run for 211390677 ->
`ecl_211390677.*`), `sec_211892898.py/.txt/.png` (secondary eclipse, mass function), `tesscut_s45_366568535.py/.txt/*.png`
(own FFI photometry S44-46, S72), `s72_and_c5_366568535.txt`, `qlp_s72_366568535.txt`, `vsx_ebs.py/.txt/*.png` (three VSX
EBs), `completeness.py/.txt`, `completeness_misses.txt/*.png`, `gaia_nss.py/.csv`, `gaia_source.csv`, `vizier_cone*.txt`,
`vizier_rows.py/.txt`, `catalog_ids.txt`, `simbad_bib.txt`, `vsx_detail.txt`, `sandquist2021.txt` (text of
arXiv:2011.08365), `pipeline_verdicts_211892898.txt`.

## Verdicts, worst first

| # | Claim | Verdict |
|---|---|---|
| 1 | EPIC 211892898: eclipse period 72.506 d, aliases 15.99/18.13/36.25/72.51 survive, matches Gaia SB1 | **CONFIRMED WITH CAVEATS.** Period right (72.5086 +- 0.0013 d) and it is the **only** alias: K2 alone excludes every sub-multiple. The alias list in the triage is wrong (pipeline failure, P1). The orbital period is **already published** (Gaia DR3 NSS, 72.534 +- 0.029 d) and the eclipses (primary *and secondary*, "P > 50 d") were reported by Nardiello et al. 2016. New: the eclipse ephemeris, the proof that the eclipses sit on the Gaia primary conjunction, and the companion size. |
| 2 | EPIC 211390677: aliases 61.442 / 92.163 / 184.327 d; S45 QLP event "marginal"; companion "probably stellar" | **CONFIRMED WITH CAVEATS, stronger than claimed.** S45 dip is real and on-target (own FFI photometry, 6.5 ppt, SNR 34, difference-image centroid 4" from target). A **fourth, partial eclipse** sits at the end of S72 exactly where the ephemeris puts it (3284.2-3285.5, 5-7 ppt; also in the QLP product the pipeline used). 61.44 d is **excluded** by K2 C5 (predicted eclipse at 213.2, flat). Survivors: **92.154 or 184.309 d**. Companion is stellar by necessity (6 ppt on a 6.4-Rsun red-clump giant => R2 ~ 0.5 Rsun). Nothing published. |
| 3 | EPIC 211411112 validation: ours 62.595 vs published 62.593 (arXiv:2011.08365) | **CONFIRMED WITH CAVEATS (mis-citation).** 62.593 +- 0.003 d is Geller et al. 2021 (RV). Sandquist et al. 2021 (AJ 161, 59 = arXiv:2011.08365) give **P = 62.59484 +- 0.00003 d**, eclipse times BJD 2457178.9531 and 2458117.8789, and say the eclipses are **secondary** eclipses. Our grid times (step D/4 = 1.1 h) reproduce these to < 30 min and predict S46 to 13 min. Valid as validation; quote 62.59484. |
| 4a | EPIC 202072917: VSX 3.842 d wrong; ours 34.746 d | **CONFIRMED, not novel.** Three identical 237-ppt TESS eclipses on P = 34.74565 d from the K2 epoch; VSX period (Armstrong et al. 2016 ML catalogue) impossible (32 predicted eclipses in TESS data, none seen). The 34.7457-d period is already in the SPOC S1-46 multi-sector TCE (MES 419; listed in Valizadegan et al. 2025, J/AJ/170/287, labelled FP). |
| 4b | EPIC 202071902: VSX 8.697 d wrong; ours 22.904 d | **ALREADY PUBLISHED (orbit).** Kovalev et al. 2024 (MNRAS 527, 521, LAMOST-MRS SB2) give a Keplerian RV orbit P = 22.859 +- 0.015 d, e = 0.37, K = 49 km/s. Our eclipse period 22.90415 d is 3 sigma from theirs but far more precise (151 cycles). VSX's 8.697 d is from Armstrong et al. 2016 (rotation/RS CVn, probably). P/2 excluded (no dips at half-period epochs in S44/S45). |
| 4c | EPIC 202126877: "21.775 / 21.952 / 43.550, pipeline confused" | **WRONG as listed; resolved here.** P = **21.7753 d**, eccentric, secondary eclipse 7.24 d after primary (phase 0.332). 21.952 d excluded (4 covered epochs with no dip); 43.55 d is 2P (even/odd primaries identical per product). VSX 3.877 d (Armstrong 2016) impossible. LaCourse et al. 2015 list it as a known C0 EB without a period. No published period found. |
| 5 | Completeness (10 sampled low-support events) | **8 fine, 1 real miss, 1 false alarm on my side.** EPIC 246194159: pipeline dropped two 25-30 %-deep eclipses (K2 C12 810.91, SNR 34; S92 3813.93, SNR 63) because of the K2 beta <= 2 cut and the one-sided-baseline rule, so its period list is wrong. EPIC 240410914: my extra "dips" are scattered-light ramps in own FFI photometry; the pipeline was right to ignore them. |

## Pipeline-level problems found (new, after the first review)

- **P1. beta is identically 1 for long events, so exclusion fails on giants.** `alias_pipeline.red_noise_beta` returns 1.0 when
  the baseline has fewer than 3 x n_bin points. With BASE_W = 0.75 d per side (73 K2 cadences) and n_bin = 0.7 D / cadence,
  beta is dead for D > ~0.5 d at 30-min cadence (and D > ~17 h at 10-min FFI cadence). The window test then demands
  |d| < 3 sigma_white (0.03-0.05 ppt here) and any 0.3-1 ppt granulation/systematic wiggle makes the window "ambiguous".
  EPIC 211892898: 264 of 325 K2 windows are "ambiguous", only 50 "excluded" (`pipeline_verdicts_211892898.txt`).
  Result: aliases that the data exclude by eye (6-ppt, 2-day eclipses missing at 9 predicted epochs) are reported as
  surviving. Same failure gives the spurious 61.44-d alias of 211390677.
- **P2. Partial eclipses at sector/campaign edges are never tested** (`test_window` needs baseline on both sides). The
  S72 ingress of 211390677 is a 6-ppt event in the QLP file the pipeline loaded, and it was invisible to it.
- **P3. The M6 fixes over-reject.** K2 beta <= 2 and the "n_oot >= 0.5 BASE_W/cadence on both sides" rule discard a 34-sigma
  and a 63-sigma eclipse of 246194159 (`completeness_misses.txt`). A beta cut should not apply to detections at SNR >> 10.
- **P4. Repeat times are grid points (step D/4) but are used as measurements.** 211892898's pipeline times are up to 0.2 d
  off the fitted mid-times; the quoted 72.506 d carried an implied error of ~0.03 d, so comparing it with Gaia's 72.534 d
  was meaningless either way (see 1d).

## 1. EPIC 211892898 = TIC 184895645 (Gaia DR3 660920219027943552)

**(a) Eclipse times, re-derived** (`ecl_211892898.py`: trapezoid x quadratic baseline, residual-bootstrap errors):

| K2 | t0 (BTJD) | depth (ppt) | T14 (h) |
|---|---|---|---|
| C5 | 189.3319 +- 0.0164 | 6.92 +- 0.19 | 55 +- 1.4 |
| C16 | 1131.9402 +- 0.0168 | 5.90 +- 0.24 | 52 +- 1.6 |
| C18 | 1276.9610 +- 0.0101 | 6.20 +- 0.10 | 51 +- 0.7 |

Pipeline grid times were 189.352 / 1132.149 / 1276.810, i.e. 0.02-0.21 d off (P4). Linear ephemeris, cycles 0/13/15:
**P = 72.5086 +- 0.0013 d, T0 = 189.3316**, max |O-C| = 0.003 d. All sub-multiples P/2 ... P/14 fit the three times
equally well (they must), so the period comes from exclusions, not from timing.

**(b) Aliases.** Box test (in-eclipse +-0.35 D vs linear baseline, red-noise allowance x2) at every predicted epoch in K2
C5/C16/C18 and TESS S72 (`ecl_211892898.txt`): 36.25 d has 3 flat predicted epochs (153.08, 1168.20, 3270.95), 24.17 d
has 4, 18.13 d has 9, 14.50 d has 12, shorter ones more. **Only 72.51 d survives, from K2 alone, with no help from Gaia.**
The triage's "15.99/18.13/36.25 survive" is a P1 artefact (15.988 d is not even a sub-multiple of the fitted period; it
came from the quantised grid times).

**(c) Gaia NSS conjunction test** (`gaia_nss.csv`: P = 72.5336 +- 0.0288, Tp = BTJD 418.905 +- 0.254, e = 0.480 +- 0.019,
omega = 201.5 +- 1.8 deg, K1 = 22.75 +- 0.85 km/s; 20 000-sample MC, errors independent, correlations ignored):

| conjunction | predicted (BTJD) | observed | O-C |
|---|---|---|---|
| primary eclipsed (nu+omega = 90) | 190.24 +- 0.66 / 1133.18 +- 0.71 / 1278.25 +- 0.74 | 189.33 / 1131.94 / 1276.96 | -0.9 to -1.3 d (1.4-1.8 sigma, common offset) |
| secondary eclipsed (nu+omega = 270) | 206.3 / 1149.3 / 1294.3 | | -17 d (> 30 sigma) |

The K2 eclipses are the Gaia primary conjunction (giant eclipsed by the companion) within ~1.6 sigma; the alternative is
excluded. The common -1.2 d offset is what a 1-sigma error in Tp/omega does. Bonus: at the three predicted *secondary*
conjunctions (all inside K2 coverage) the flux is lower by 0.51 +- 0.06, 0.82 +- 0.06, 0.88 +- 0.07 ppt
(`sec_211892898.txt`; folded 3-d bins elsewhere scatter by ~0.1 ppt, worst outlier -0.26). A 0.5-0.9 ppt secondary there
is what a 0.55-Msun dwarf behind a 7-Rsun giant should give. Nardiello et al. 2016 already report "one primary and one
secondary eclipse" in C5. Treat as supportive, not as a detection claim, until detrended properly.

**(d) 72.506 vs 72.534 d:** ours 72.5086 +- 0.0013 vs Gaia 72.5336 +- 0.0288 -> 0.87 sigma. **Not significant.**

**Physical consistency.** f(M) = 0.0597 Msun. Asteroseismology: APO-K2 (Schonhut-Stasik et al. 2024, J/AJ/167/50) M = 1.11
+- 0.09 Msun, R = 6.82 +- 0.18 Rsun, RGB; Khan et al. 2023 (J/A+A/677/A21) M = 1.73, R = 8.44; Zinn et al. 2022 numax 81 uHz,
RGB. For M1 = 1.1-1.2: M2 sin i ~ 0.55 Msun, a ~ 86-89 Rsun, R1/a ~ 0.08-0.10, central-transit duration 47-58 h (observed
51-55 h), (R2/R1)^2 with R2 = 0.55 Rsun, R1 = 6.8 Rsun = 6.5 ppt (observed 6.2-6.9). A K/M dwarf transiting a red giant on an
e = 0.48, 72.5-d orbit. Companion is stellar.

**(e) Literature** (`vizier_cone*.txt`, `vizier_rows.txt`, `simbad_bib.txt`, web searches):
- SIMBAD: 0 references. Not in VSX. Not in Gaia DR3 vari_eclipsing_binary or the Gaia classifier. Not in Kruse et al. 2019
  (J/ApJS/244/11 absent from a 5" cone over all VizieR catalogues), not in Armstrong et al.
- Nardiello et al. 2016 (MNRAS 463, 1831, arXiv:1608.00459) Sect. 4.3.1: "eclipsing or transiting field object with a
  period greater than 50 d", one primary and one secondary eclipse, ~2-d duration, RV follow-up needed. (The fetched text
  quotes a primary depth "0.085 Kp mag", impossible for a 7-ppt event; read the PDF before citing.)
- Gaia DR3 NSS (2022): SB1 orbit, P = 72.534 d. A published orbital period.
- Tian et al. 2020 (ApJS 249, 22): LAMOST RV-variable candidate (3 epochs, dRV = 19 km/s), no period.
- Rebull et al. 2017 (J/AJ/154/224): Praesepe candidate, Prot = 8.51 d (it is a 1.1-kpc background giant, plx 0.88 mas).
- TESS: only S72 covers it (ecliptic latitude 0.46 deg). No primary (3234.6, 3307.1) or secondary (3252.7, 3325.3) falls in
  S72 (3260.5-3285.6), so TESS neither confirms nor contradicts. Not checked: ADS full-text, ExoFOP target page.

Bottom line: the eclipse ephemeris and the Gaia-orbit identification are new and solid; the orbital period is not.

## 2. EPIC 211390677 = TIC 366568535

- Re-timed K2 eclipses: 151.7809 +- 0.0055 and 1257.6333 +- 0.0086 (depth 6.3 / 5.9 ppt, T14 39 h). dT = 1105.852.
- **S45 event is real and on-target** (`tesscut_s45_366568535.txt`): own TESScut photometry depth 6.45 / 6.56 / 6.61 ppt for
  1.0 / 1.6 / 2.5-px apertures, white SNR 30-34 (QLP had 21); difference-image peak SNR 14.5, centroid 0.18 px = 4" from
  the target; per-pixel depths peak on the target pixels. Mid-time 2547.79 (box); the pipeline's 2548.04 is a grid point.
- **Fourth eclipse (partial) at the end of S72.** Every alias through S45 predicts an eclipse at 3285.03 +- 0.03
  (ingress 3284.21); S72 ends 3285.59. Own photometry: baseline flat to +-0.2 ppt, then a step to -5 to -7 ppt from 3284.3
  on (`tesscut_s72_366568535_lc.png`, `s72_and_c5_366568535.txt`); the QLP file shows the same (0.2-d bins: +0.2, -2.2,
  -5.4, -6.0, -6.6, -6.8, -5.4, -5.5 ppt; `qlp_s72_366568535.txt`). Caveats: last 1.3 d of a sector, scattered light rising
  (other stars in the cutout brighten by 2-4 ppt over the same window), pixel-level drop centroid 0.6 px off. The timing
  coincidence (predicted to 0.03 d over 3100 d) is hard to dismiss. It does not split the aliases.
- **Aliases of the pair (my coverage test, K2 C5/C18 + TESS S7/44/45/46/72):** without S45, ten survive (n = 1-6, 8, 10,
  12, 15). With S45: n = 6 (184.309 +- 0.002), 12 (92.154 +- 0.001), 18 (61.436), 24, 30, 36 ... The 61.44-d alias predicts
  a full eclipse at 213.22 (T1 = 212.40, T4 = 214.04; C5 runs to 214.41): the data drift up by ~1 ppt and show no 6-ppt dip,
  so it is excluded (pipeline: "ambiguous", P1). n = 24, 30, 36 are excluded by K2 (`ecl_211390677.txt`).
  **Survivors: 92.154 d or 184.309 d.** Next eclipses: 3377.2 (92-d only) and 3469.3 (both).
- **Companion is stellar.** Host is a red-clump giant (Zinn et al. 2022 EV = RC; Khan et al. 2023 R = 6.41 +- 0.25 Rsun,
  M = 1.14 Msun; LAMOST logg 2.87). A 6-ppt eclipse needs R2/R1 >= 0.077, i.e. R2 >= 0.5 Rsun; a brown dwarf or planet
  (~0.1 Rsun) would give 0.25 ppt. Gaia: no NSS solution, RUWE 1.15, rv_nb_transits 17 but no amplitude (G = 13.3).
- Literature: SIMBAD 0 refs; not in VSX; Green et al. 2023 (J/MNRAS/522/29, ellipsoidal search) list it with score 0.0 and a
  meaningless Porb = 5.998 d. Not checked: ExoFOP page, ADS full text.

## 3. EPIC 211411112 (WOCS 11028 / Sanders 617)

Sandquist et al. 2021, AJ 161, 59 (arXiv:2011.08365; `sandquist2021.txt`): K2 C5/C16/C18, single eclipses in C5 and C16,
"no eclipse was observed in campaign 18 ... one was not expected"; eclipse times BJD 2457178.9531 +- 0.00088 and
2458117.8789 +- 0.00075; "these are secondary eclipses (of the cooler, less massive star)"; joint solution P = 62.59484 +-
0.00003 d, e = 0.6227, i = 86.7-86.9 deg, masses 1.222 / 0.909 Msun. Geller et al. 2021 (AJ 161, 190, J/AJ/161/190): P =
62.593 +- 0.003 d from RVs. Our grid times: 178.9526 (catalogue), 1117.861 (26 min from Sandquist), S46 2557.5664
(predicted 2557.557: 13 min). Our 62.59517 differs from 62.59484 by 0.0003 d = 11 of their sigmas, but that is the
D/4 = 1.1-h grid over 38 cycles. No TESS detection of this eclipse is published (web search), but it is a trivial prediction
of a solved system. Keep it as the validation case; fix the citation and the number.

## 4. VSX corrections (`vsx_ebs.txt`)

All three VSX entries trace to Armstrong et al. 2016 (MNRAS 456, 2260; K2 variable catalogue II, machine-learning
classification of C0-C4), uploaded 2019 (202071902 revised by K. Bernhard, "period from 2016MNRAS.456.2260A"). Those
periods are periodogram outputs on short C0 light curves, not eclipse periods, and the VSX amplitudes (0.025 / 0.2 / 0.008
mag) are not the eclipse depths either. "VSX is wrong" is right, but the correction is against an automated catalogue.

- **202072917 (TIC 308255000):** TESS eclipses at 1489.1650, 2496.7871, 2531.5321 (237 ppt, T14 6.4 h, identical), all on
  P = 34.74565 d from the K2 epoch (O-C < 0.05 h). No secondary anywhere (phase scan max 6 +- 3 ppt; phase 0.5: 0.5 ppt).
  Odd/even identical, so 2P would need twin stars; P/2 excluded by empty epochs. SPOC multi-sector TCE (S1-46) already has
  34.7457 d (J/AJ/170/287). Gaia: Teff 13 000 K, rv_amplitude_robust 315 km/s (hot-star RVS, unreliable).
- **202071902 (TIC 57051340):** TESS eclipses 2517.1644 (S44), 3250.0969 (S71), 91 ppt, 4.6 h; P = 22.90415 d. No secondary
  (< 2 ppt at any phase). Only odd cycles observed, but P/2 epochs in S44/S45 are empty. Kovalev et al. 2024 orbit:
  22.859 +- 0.015 d, e = 0.37, SB2.
- **202126877 (TIC 86491361):** 11 TESS dips: primaries 2488.583, 2510.368, 2532.147, 3250.695, 3272.484 (205 ppt in
  TESS-SPOC S43-45, 285 ppt in QLP S71-72: product dilution), secondaries 2495.827, 2517.592, 3236.092, 3257.846, 3279.603
  (43 / 127 ppt, 8-10 h), plus one 12-h, 45-ppt feature at 2480.06 (S43 start, probably systematics). Primary period
  21.7753 d (refit from the K2 epoch, O-C < 0.5 h); secondary at +7.24 d = phase 0.332 (eccentric). "Dip pairs 7.25 d
  apart every 21.78 d" is exactly right.

## 5. Completeness spot-check (`completeness.txt`, `completeness_misses.txt`)

Sample: the 3 deepest and 3 highest-SNR max_support <= 1 events, 3 highest-SNR support-1 events, 2 random zero-K2-repeat
events (10 after de-duplication). Independent sliding-box search on the same combined light curves.

- 8/10: no missed repeat. (Catalogue depths are 1.2-1.7x the box depths throughout, a known convention issue.)
- **246194159 (TIC 301328232), real miss:** five 25-30 % eclipses (K2 747.41, 810.91; TESS 3228.00, 3813.93, 3823.68),
  pipeline kept two. 810.91 passed the window test (supported, SNR 34, cover 0.88) and was dropped by the K2 beta <= 2
  rule; 3813.93 (SNR 63) was dropped because the right-hand baseline had 134 < 162 points. The reported support-1 period
  list is wrong for this star. Shapes differ (3823.68 is narrower): eccentric EB or two eclipsing sources; not pursued.
- **240410914 (TIC 201980990):** my search flagged 3788.98 and 3802.91 in own FFI photometry; both are scattered-light
  ramps at segment starts (`completeness_240410914.png`). The pipeline was right. The K2 "event" itself is a 2-3 cadence
  spike to negative flux; the catalogue entry looks like an artefact.
- Systemic: partial eclipses at data edges are never counted (P2), so zero-repeat counts are lower limits.

## Not checked / limits

- ADS full-text and ExoFOP per-target pages were not queried (no ADS token; ExoFOP not fetched). The literature check used
  SIMBAD (query_bibobj), a 5" VizieR cone over all catalogues, VSX, Gaia DR3 tables and web search.
- Gaia NSS error propagation ignored parameter correlations (corr_vec not decoded).
- Nardiello et al. 2016 depth values come from a fetched HTML summary, not from the PDF.
- The S72 partial eclipse of 211390677 has no post-egress baseline by construction; the 5-7 ppt amplitude depends on the
  background model.
- No TESS TCE/QLP-alert lookup for TIC 86491361 or TIC 57051340.
