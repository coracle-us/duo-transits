# Blind re-derivation of K2 single-event periods (independent_table/)

This is a blind analysis by a separate agent, written from scratch without access to the project's
pipeline results, triage, reviews or catalogue checks. It covers 14 stars: 12 leads (including the
validation star Sanders 617) and 2 decoys with no pipeline repeats. Saved by the main session from
the agent's hand-back; the agent itself could not write this file.

Contents: `results.csv` (one row per star), `scripts/` (all code), `plots/`, `data/final_details.json`
(per-event, per-product fits, the epochs that exclude P/2, P/3 and P/4, and phase-0.5 checks) and
`data/allowed_windows_*.csv`. `data/` (1.4 GB) is gitignored.

## Summary

12 of the 14 stars show the K2 event again, in later K2 campaigns and/or TESS. Only 230887315 and
212152316 (the two decoys) are still single. Three stars have distinct secondary eclipses, so they
are eclipsing binaries.

| EPIC | Events found (BTJD) | Allowed P [d] | 2× possible? | Confidence |
|---|---|---|---|---|
| 248811085 | 948.3792 (K2 C14), 2573.6329 (S46), 3263.1346 (S72); ~9.0% deep, 5.7 h | **49.25012 ± 0.00002** | only if twin eclipses; depths equal to ~1%; no dip at phase 0.5 | high |
| 246696483 | 855.588 (C13), 2495.7952 (S43), 3224.7760 (S70); ~5.5% (QLP/PDCSAP), 16.6 h | **182.2452 ± 0.0003** | only if twin | medium-high: P/2 is excluded by a single S98 epoch (4044.88, flat, checked by eye) |
| 247967714 | 845.0776 (C13), 2515.7176 (S44), 3258.2245 (S71); 25%, V-shaped | **185.6267 / 92.8133 / 61.8756** | yes if twin; probably an EB | high for events; period 3-fold ambiguous |
| 246089124 | 774.4696 (C12), 3215.2285 (S70), 3806.9285 (S92); 13%; plus a **secondary** at 810.986 (0.64%, phase 0.494) | **73.96241 ± 0.00002** | no; it is an EB | high |
| 211892898 | 189.312 (C5), 1131.939 (C16), 1276.962 (C18); 0.62–0.69%, ~53 h | **72.5100 ± 0.0014** | only if twin | high |
| 211503363 | 197.5175 (C5), 2532.717 (S45), 3270.136 (S72); 0.31–0.37%, 23 h | **122.9049 ± 0.0007 or 61.4525 ± 0.0003** | if twin | medium: each TESS dip is ~10σ |
| 211390677 | 151.7865 (C5), 1257.634 (C18), 2547.824 (S45), plus a partial at the S72 end (3284.96); 0.6%, 39 h | **184.3103 ± 0.0010 or 92.1552 ± 0.0005** | if twin | high for events; period 2-fold ambiguous |
| 211411112 | 178.9534 (C5), 1117.8758 (C16), 2557.556 (S46); 2.6–3.0%, 3.7 h | **62.59479 ± 0.00005** | only if twin | high |
| 211489484 | 6 primaries (C5, C18, S7, S34, S46, S71) and 3 secondaries (208.606, 1277.968, 3252.17) | **82.25843 ± 0.00007** | no; eccentric EB, secondary at phase 0.194 | high |
| 202071902 | −208.4266 (C0), 2517.1648 (S44), 3250.0970 (S71); 9%, 4.7 h; plus an S72 ingress | **22.904130 ± 0.000002** | only if twin | high |
| 202126877 | 7 primaries (C0, S43–45, S71–72, S1751) and 5 secondaries | **21.77521** (mean; timing variations ±0.017 d) | no; EB | high |
| 202072917 | −213.3695 (C0), 1489.1650 (S6), 2496.7874 (S43), 2531.5334 (S45), 4060.341 (S1751); 24%; plus an S72-edge eclipse | **34.745606 ± 0.000007** | only if twin; consecutive events are 23.7% and 23.8% deep | high |
| 230887315 | only 660.3097 (C11a); 3.0%, V-shaped | P > 70.0 d (86% of 70–200 d allowed) | n/a | high that no other event is in the data |
| 212152316 | only 207.1933 (C5); 6.3%, 29 h, V-shaped | P > 68.5 d (29% of 68.5–200 d allowed) | n/a | high that no other event is in the data |

## Surprises

- **202126877's eclipses are getting deeper.** The primary was 6% deep in 2014 (K2), 20.5% in 2021,
  28% in 2023 and 32% in 2025. The secondary went from 4.3% to 11–12%. Eclipse times vary by
  ±0.017 d and the secondary phase drifts (0.333 → 0.328). This looks like a precessing orbit,
  probably driven by a third body.
- **246089124 has a K2 secondary eclipse** at phase 0.494, seen in PDCSAP, EVEREST and K2SFF. The
  13% event is therefore an EB, not a planet.
- **Two predicted epochs were confirmed in data the pipelines had dropped**, in cadences taken
  during scattered light:
  - 202071902: an ingress at 3273.0023, predicted 3273.0011 (TARS and own unmasked photometry);
  - 202072917: a full 19% eclipse at 3261.1888, predicted 3261.1912 (own photometry at the start
    of S72).
- **Catalogue depths can be badly off.** For 246696483, K2SFF and wide TESS apertures give 3.1–3.3%,
  but PDCSAP and QLP give ~5.5%. 202126877's catalogue depth reflects the 2014 state of a changing EB.
- **"Sector 1751"** in tess-point/TESScut is a real 7-day segment (BTJD 4055.8–4063.0) and
  contributed events.

## Method in brief

For each star:
- All K2 campaigns, found by EPIC ID and by coordinates. This added C16/C18 for several C5-field
  stars. For each campaign: PDCSAP, EVEREST and K2SFF, plus raw SAP for inspection only.
- Every TESS pipeline product on MAST, plus own aperture photometry on TESScut cutouts for every
  sector.
- A local-polynomial box scan at the catalogue duration, with a red-noise error term. Every
  detection was checked by eye, and the C5-field stars were cross-checked against each other for
  shared systematics.
- "Covered and flat" means:
  - at least 50% of the expected in-event cadences are present;
  - the expected depth is excluded at more than 3σ, and the measured depth is below half of it.
    For TESS the expected depth is taken as half the catalogue depth, to allow for dilution;
  - the product is a trusted reduction. For events longer than 0.5 d, products that detrend away
    long dips are not allowed to exclude;
  - no product shows a detection at that time.
- Candidate periods are Δt/n. A candidate is excluded if any predicted epoch lands on a
  covered-and-flat epoch. The decisive epochs for excluding P/2 were checked by eye.

## Catalogue check (done only after `results.csv` was saved; `data/catalogue_check_after_blind.txt`)

- VSX lists the three C0 stars as EBs with P = 8.697 d (202071902), 3.877 d (202126877) and
  3.842 d (202072917). Those periods contradict the 36 days of C0 data, which show a single
  eclipse each. They look like wrong automated classifications.
- SIMBAD has 211411112 as an M67 member (NGC 2682 Sanders 617) and 230887315 as HD 154148.
- No other VSX entries. ExoFOP and the literature were not checked.

## Per-star notes

- **248811085**: P/2 is excluded by K2 epochs 923.75 and 973.00 (80 ± 113 and −57 ± 26 ppm,
  against a 9% expected depth).
- **246696483**: P/3 and P/4 are excluded by S5, S32 and S98.
- **247967714**: every other epoch of the 61.88-d alias falls in data gaps; P/4 and shorter are
  excluded. The K2 SAP "event" at 821.3 is a pointing excursion, seen at the same time in 246696483.
- **246089124**: secondary at 810.986 ± 0.001, ~6400 ppm, ~8 h long. P/2 is excluded by the flat
  epoch at 811.45.
- **211892898**: P/2 is excluded at 153.06 (C5) and 1168.20 (C16). No predicted epoch falls in S72.
  A 53-h duration at a 72.5-d period implies a large (evolved) host star.
- **211503363**: the TESS dips appear in QLP SAP and own photometry at 10–12σ. QLP's detrended flux
  removes them, and they are absent in neighbouring stars. P/3 and P/4 are excluded by C5 and C18.
- **211390677**: the partial S72 event is offset −0.11 d from prediction, which is acceptable for a
  partial fit. P/3 is excluded at 213.22 and P/4 at 197.86.
- **211411112**: P/2 is excluded at 147.66, 210.25 and 1149.17. The dips at 1136.08 and 1251.89 are
  K2 roll artefacts.
- **202071902**: the period is unique because cycle numbers 119 and 151 share no common factor.
  Predicted epochs at 2540.07 (S45) and −231.33 (C0) fall in gaps.
- **202126877**: linear ephemeris T0 = 2488.5793. The K2 event sits exactly on a primary cycle
  (n = −124). No secondary in C0 (< 0.1%).
- **202072917**: P/2 is excluded at C0 −230.74 (−543 ± 62 ppm) and at S6.
- **Caveat for all "unique" periods**: they assume every event is the same kind of eclipse. Where
  no distinct secondary was found, a twice-as-long orbit with two identical eclipses cannot be
  ruled out by photometry alone.
