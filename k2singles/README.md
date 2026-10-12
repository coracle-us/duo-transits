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

### EPIC 211924561: blind re-analysis (`independent_211924561/`, 2026-10-12)

A separate agent was given only the K2 C5 event and wrote all its own code. (Its brief carried wrong
coordinates; the agent corrected them from the data headers. Our pipeline used the right position:
same Gaia source, 0.27″.)

**Four transits:**

| Event | Tc (BTJD) | Note |
|---|---|---|
| K2 C5 | 179.8978 ± 0.0008 | |
| K2 C16 | 1163.3554 ± 0.0013 | hidden in standard products: EVEREST flags it as outliers, K2SFF has a gap; clear in raw SAP |
| TESS S46 | 2568.2983 ± 0.0030 | |
| TESS S72 | 3270.7646 ± 0.0032 | |

- The C5, S46 and S72 ephemeris predicted C16 at 1163.3555; it was observed at 1163.3554.
- All m ≥ 3 are excluded by at least two flat epochs, each in 2–3 reductions.
  **P = 140.49401 ± 0.00011 d or 70.24701 ± 0.00006 d**, with T0 = BTJD 179.89774 ± 0.00073.
- **On target:** TESS centroid offsets are about 0.1 ± 0.1 px, ruling out the G = 16.0 and 16.6
  neighbours; all others are too faint. The K2 pixel test is clean.
- **No secondary eclipse** at phase 0.5 for P = 70.25 d: −0.04 ± 0.05 ppt.
- **V-shaped and grazing:** b ≥ 0.95, and flat-bottomed fits are worse by Δχ² ≈ 20.
  - T14 = 3.0 h, against a central-transit 15.2 h (140.5 d) or 12.1 h (70.2 d) for the Gaia star.
  - A circular orbit with the Gaia density is strongly disfavoured (Δχ² 265 and 69).
  - Explanations: e ≥ 0.5 (140.5 d) or ≥ 0.31 (70.2 d), a grazing large or stellar companion, or
    a different host.
  - R2/R* is about 0.055–0.16 or more, i.e. ≈0.9–2.7 R_J or more.
- **Status:** the period is narrowed to two values with high confidence. The companion's nature is
  uncertain and possibly stellar; RVs or high-resolution imaging would settle it.
- **Next distinguishing epoch:** BTJD 4324.47 (70.2 d only), which was 2026-10-11, just missed. The
  next shared epoch is BTJD 4394.7 (about 2026-12-20).
