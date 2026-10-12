# K2 single transits re-observed by TESS

**Question.** For the single transits found in K2 (LaCourse & Jacobs 2018, RNAAS 2, 28; 164 events
in C0–C17), can TESS (2018–2026) find a repeat, or narrow the period?

**Literature check (2026-10-12).** No systematic TESS-based recovery of these periods was found.
- A 2019 paper (arXiv:1912.04287) used only overlapping K2 campaigns, for six systems.
- GOT 'EM III (Mann+2023) concerns a TESS single.
Re-check before any outreach or paper.

## Feasibility (`scripts/build_catalog.py` → `data/k2_singles.csv`)

- **Matching.** All 164 events matched EPIC → TIC within 3″. 158 stars have at least one TESS sector
  observed (ecliptic sectors 42–46, 70–72, 91–92, plus some others).
- **Sample.** About half are binary-star eclipses (flagged deep or with a secondary eclipse). 65 are
  planet-like (depth ≤ 2%, not flagged). Of those, **50** have TESS coverage and an expected
  per-transit TESS SNR ≥ 7, giving 136 sector-observations.
- **Already followed up.** None of the stars has a TOI. 3 have CTOIs and 2 are confirmed planet
  hosts; both of those are validation cases.

## Pilot plan (pre-registered)

- **Targets.** The 50 planet-like, TESS-covered, SNR ≥ 7 candidates.
- **Pipeline.** The single-transit pipeline (`singles/scripts/pilot.py`) is reused with the K2 epoch
  as reference. It runs the repeat search in all TESS data and our own TESScut photometry for
  sectors ≥ 80 that lack pipeline products. The flat-window period constraints are combined with a
  lower limit from the K2 campaign's own coverage. Difference images are made for any repeat.
- **Vetting.** As in the singles project: depth consistency, products, centroids, Gaia binarity and
  an independent re-check before any claim.

## Pilot results (2026-10-12)

- **Run.** 50 targets: 47 analysed, 3 with no TESS light curves (`results/pilot/`,
  `results_pilot_stdout.txt`). The median fraction of log P (5–1500 d) excluded is 0.38.
- **Bug fixed after the run.** The K2 period floor first used min/max over all K2 campaigns, giving
  a spurious P_min of about 1000 d for re-observed fields. `scripts/fix_pmin.py` re-applied it
  using only the campaign that contains the epoch.
- **15 stars show TESS repeats.** Triage:
  - **Systematics or variable stars:** 24395505 (5 repeats in S12), 178941241 (3 in S91 QLP,
    diff images >100″ off) and 203254318 (4 shallow, inconsistent).
  - **Off-target or weak:** 366568535 (diff 146″), 206241017 (diff 152″) and 56621708 (SNR 9.6,
    diff inconclusive).
  - **Stellar companions:** 51759568 (R* 6 R☉), 257217497 (≈2 R_J), 156599401, 369478173 (TESS
    depth 1.8× K2), 20444116 (giant host; P 122.9/61.5 d) and 425280270 (≈1.8 R_J).
    248354054 is also stellar: Gaia gives R* = 3.8 R☉, not the catalogue's 0.8, so the companion is
    ≈3.5 R_J. These give new periods for K2 single eclipses.
  - **Validation:** 35021200 = EPIC 201132684, a known K2 multi-planet host with a CTOI. The TESS
    S46 repeat is on target.
  - **Lead: EPIC 211924561 = TIC 238606535** (G 10.7 F star, Gaia R 1.74 R☉, M 1.42 M☉,
    RV-stable, RUWE 1.19).
    - Transits: K2 C5 (BTJD 179.89; K2SFF depth 1536 ± 36 ppm), TESS S46 (2568.32) and S72
      (3270.78). Both TESS events are detected in SPOC 2-min, TESS-SPOC and QLP at 1.5–2.0 ppt.
    - **P = 140.495 or 70.248 d**, with O−C of −4, +19 and −15 min. Rp ≈ 7.4 R⊕.
    - The K2 C16 epoch falls in a data gap, the C18 epoch just after the campaign ends, and the
      TESS S44 epoch (n = 34) in a gap. Nothing yet separates the two periods.
    - Concerns: T14 ≈ 3.4 h is short for the FLAME star (b ≈ 0.96–0.97, possibly grazing). The S72
      diff image is weak (SNR 3.3).
    - To do: shape fit (V vs U), independent re-analysis, ExoFOP and literature check, predicted
      transits.
