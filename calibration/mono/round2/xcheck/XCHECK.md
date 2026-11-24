# Harness cross-check before the freeze (prereg_round2.yaml, harness_check)

2026-11-22. A fresh agent re-implemented the round-2 harness from `../FROZEN_RULE.md`, the
pre-registrations and the pipeline (`st2core.py`, `st2_pipeline.py`). It did not read
`rule2.py`, `harness.py`, `dump.py`, `tune.py`, `monocal.py`, the runners or the tests, and
`harness_rows.csv` was moved out of the folder while it worked.

- Trials: `defs.csv`, 60 trials on 20 tuning-set stars (10 confirmed-planet hosts, A1 and A2; 10
  injection hosts, 4 draws each), drawn with numpy seed 20261125 (`sample_a.csv`, `sample_b.csv`).
- Its code and rows: `independent/ind2.py`, `independent/ind_rows.csv`.
- Ours: `harness_rows.csv` (`harness_rows_referr.csv` is the run before the fix below).

## Result

Agreement on 60 of 60 trials in every compared column: reference measured, reference SNR, depth
and duration, true period excluded, fraction of ln P excluded, number of vetted events, false
event, repeat detected, clean repeats, alias event time, call, number of live aliases, truth in
the alias set, and for arm B the null-run event flag and all the own-events-removed columns.

Limits of this check: none of the 60 trials has a false exclusion, 15 have a call, and the two
codes share the pipeline's functions, so it tests the harness around the pipeline, not the pipeline.

## What the agent found in the spec

1. The spec did not say which depth error the alias table uses. Our harness used the
   reference's; the pipeline (and round 1's third correction) uses the second event's. Fixed in
   `harness.py` and stated in the spec. On these 60 trials both readings give the same calls.
2. Sixteen smaller ambiguities (strict or inclusive limits, which duration is used where, the
   null run's own reference). Its literal readings matched ours in every case. The ones that
   could matter are now written into the spec.
3. Differences from the pre-registration that the spec now states: probabilities are binned by
   SNR only; the null run is a scoring device, not a veto.
