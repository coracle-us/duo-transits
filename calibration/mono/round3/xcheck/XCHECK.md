# Round 3: harness cross-check before the freeze

A fresh agent re-implemented the round-3 scoring from `FROZEN_RULE3.md`, round 2's `FROZEN_RULE.md`,
the two pre-registrations and the pipeline, without reading `round3/` code, the round-2 harness or
the tests (`independent/ind3.py`, started from the round-2 blind implementation). It did not see
the harness rows; they were kept outside the repository until it had finished.

Trials: `defs.csv`, 120 trials on 40 already-seen stars - 20 drawn at random (seed 20261206) and 20
drawn among stars with a call (seed 20261207), half arm A and half arm B (arm B with the 10%
catalogue error already applied).

Result: **120 of 120 trials agree in every compared column** - ref_snr, ref_measured, n_vetted,
false_event, repeat_detected, t2, call, n_alias_ok, n_alias, truth_in_alias_set, event_true,
frac_alias_excluded, frac_alias_excluded_r2, and the same eleven call columns with own events
removed on the 80 arm-B trials. Numbers agree to 1e-9.

Limits:
- Alias tables were compared through their counts and excluded fractions, not status by status.
- Among the 120 trials there are 42 narrowed, 1 resolved and 2 all-excluded calls.
- The catalogue error of arm B (`harness3.add_noise`) was given to the agent already applied, so
  it is not cross-checked.
- Both implementations call the same pipeline functions on the same cached light curves.

The agent's ambiguities are in `independent/AMBIGUITIES.md`. Four led to wording changes in
`FROZEN_RULE3.md` (no change of behaviour): which duration `sec_box_snr` uses, which sector data,
that the alias event is simply the highest-SNR passing event, and the note on multi-candidate stars
and unmeasured references.
