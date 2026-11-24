# Driver check (prereg_catalogue.yaml, driver_check), run 2026-11-23 before the freeze

Only calibration stars were loaded. `run_catalogue.py --check` refuses any star of `catalogue/defs.csv`.

## 1. Against round 3's rows

`run_catalogue.py tables --check a2_check_defs.csv` (round 3's 116 arm-A2 trials, every passing event
taken as accepted) against `calibration/mono/round3/test/arm_a_rows.csv`. Output in `out_a2/`.

| Column | Equal (exactly) |
|---|---|
| fraction of ln P excluded | 116 of 116 |
| reference SNR, reference measured | 116 of 116 |
| number of passing events | 116 of 116 |
| call, second-event time | 116 of 116 |
| aliases, aliases allowed or supported, fraction of aliases excluded | 116 of 116 |
| true period excluded by step 1 (`period_excluded`), and its number of epochs | 116 of 116 |
| true period in the alias set (`alias_status`), called trials | 63 of 63 |

No sector failed to load.

## 2. Duotransit table against a blind re-implementation

20 one-TOI round-3 stars whose second transit was recovered; `t2_cat` is the true epoch nearest that
event (`duo_check_defs.csv`). A fresh agent wrote `independent/ind_duo.py` from the pre-registration
and the two frozen rule files, without reading `catrun.py`, `rule3.py` or either harness
(`independent/AMBIGUITIES.md`). Driver output in `out_duo/`.

- Alias statuses: **24,966 of 24,966 equal**; periods equal to the last digit.
- The true period: not excluded by step 1 in 20 of 20; in the duotransit table supported in 19, allowed in 1.

Limits:
- In all 20 rows the second transit's depth error came from a pipeline event at `t2_cat`. The two
  fallbacks (refit at `t2_cat`; the reference's error) were not exercised by either implementation.
- No row is on a multi-candidate star, so the wider window (w0 = 8) was not exercised here. It was in
  round 3's cross-check.
- The "third event" columns have no independent check; they are covered only by a synthetic test.
- Both implementations share the loading, detrending, reference and event-search functions.

## 3. The gate

`run_catalogue.py search` without `--check`, before the freeze: "The catalogue run is not frozen:
catalogue/APPROVALS.md has no step-2 approval. Nothing was loaded."
