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
