# Independent re-derivation of the single-transit calibration (2026-11-22)

A fresh agent re-implemented both arms from `prereg_mono_calibration.yaml` (version 1.1), the two
sample files and the pipeline code. It was told not to open `monocal.py`, the drivers, `score.py`,
the tests, `RESULTS.md`, `PROPOSED_v1.1.md`, `results/` or `calibration/README.md`, and reports that
it did not. Its code and per-trial tables are in `independent/`; its tables are
`independent/score_tables.md`. It could not write its own report file, so this file is written by
the first analyst from the agent's hand-back message.

**Same verdict: fail.** It also found three places where the first harness did not do what the
pipeline does. Those are fixed and both arms were rerun into `results_corrected/`.

## The three corrections to `monocal.py`

| # | First harness | Pipeline (and the independent run) | Effect |
|---|---|---|---|
| 1 | A reference too weak to measure got a flat 20% depth error. | `search_one` uses the catalogue's depth error when there is one. | 3 more false exclusions in arm A2 (TOI-2449, TOI-4504, TOI-1154): 4 of 109 became **7 of 109**. |
| 2 | A second event earlier than the reference gave an empty alias table, labelled "all excluded". | Times are passed in order. | Arm B calls rose from 178 to 201. |
| 3 | The alias table used the reference's depth error. | `step_aliases` uses the second event's depth error. | A few alias statuses. |

## Numbers side by side

| Metric | First run (`results/`) | Corrected (`results_corrected/`) | Independent |
|---|---|---|---|
| A2 false exclusion | 4 of 109 | 7 of 109 (6.4%) | 7 of 109 |
| A2 false exclusion, SNR 8-12 | 2 of 17 | 2 of 17 | 2 of 17 |
| A2 false exclusion, SNR above 20 | 2 of 38 | 2 of 38 | 2 of 38 |
| B false exclusion | 4 of 515 | 4 of 515 (0.8%) | 4 of 516 |
| A2 wrong period | 16 of 75 | 15 of 74 (20%) | 16 of 74 |
| B wrong period | 65 of 178 | 84 of 201 (42%) | 83 of 200 |
| A1 calls (all false by construction) | 47 | 48 | 47 |

Arm A uses the same stars, reference times and refined periods in both analyses, so it compares
star by star: after the corrections the excluded-or-not verdict agrees on 218 of 218 trials and the
call on 214 of 218. The 4 remaining differences (TIC 119584412, 243014114, 4646810, 66561343) have
not been traced. Arm B uses different random draws and is compared as rates.

**Correction to the first report.** `RESULTS.md` as first written said false exclusion "passes
overall" in A2. With the pipeline's depth error it does not: 7 of 109 is above the 5% limit.

## What the independent run adds (its classification, from numbers only, not from light curves)

- 5 of the 7 A2 false exclusions look like ephemeris problems: TOI-216, TOI-2525 and TOI-4504
  (which it calls known large timing-variation systems; not checked against the literature here),
  TOI-2449 (reference SNR 1.6, period at the edge of the refinement grid) and TOI-7510 (3-TOI star).
  The other 2 are the rule firing on a brightening at a real transit epoch.
- **The exclusion rule has no floor.** Depth below 0.5 x depth0 - 3 sigma lets a brightening or a
  bad stretch exclude a period. All 4 arm-B false exclusions are single epochs with a strongly
  negative box depth, on three hosts fainter than Tmag 12 with a 30-minute QLP reference sector.
- In arm B, 110 of 303 distinct false events recur in two or more of a host's four draws: they are
  features of the host's real light curve. The hosts were not clean.
- A false event passes the cuts on 26 of 76 single-TOI stars in A1.
- `core.alias_table` returns an empty table when t2 < t1, silently.
- The pre-registration says QLP first; `best_per_sector` prefers SPOC. Both analyses followed the code.
- The 21 arm-B hosts without light curves have no SPOC, TESS-SPOC or QLP product. Neither analysis
  made its own TESScut photometry, which the pre-registration lists as the last fallback.

## Limits

- The independent run shares the samples, the pipeline code and the pre-registration with the first
  analysis. It checks the harness and the scoring, not the pipeline or the choice of test.
- The corrected rerun was made after the outcome was known.
