# Round 2 tuning log

Pre-registration: `../prereg_round2.yaml` version 1.0 (step 1 approved by Jon, 2026-11-22).
**Status: tuning done on the tuning set only. No rule is frozen. No test-set star has been loaded.**
Every number below is from stars already seen in round 1, so none of it is a result.

## Tuning set

| Set | Stars | Trials | Source |
|---|---|---|---|
| A1, A2 | 109 confirmed-planet hosts (76 with one TOI, 33 with more) | 218 | round 1's `arm_a_sample.csv`, same reference transits |
| B | 129 frozen hosts (21 of the 150 have no light curve) | 515 | round 1's 4 draws per host, same seed, same planets |
| Bx | the same 129 hosts | 1032 | 8 new draws per host, depth log-uniform 500-5000 ppm, seed 20261124 |
| null | the same 129 hosts | one per draw | the star searched with nothing injected, with that draw's depth and duration |

The 83 stars of the 2026-11-01 test were not used: they ran through a different harness.

Code: `dump.py` stores each trial's box depths and events once (`calibration/mono/work/round2_tuning/`,
not in git); `rule2.py` holds the rule with its tunable parts; `tune.py` scores variants
(`tuning/*.csv`). Check: with the round-1 setting the stored summaries reproduce round 1's arm A
verdicts on 218 of 218 trials and its calls on 213 of 218 (the 5 differences come from looking alias
epochs up on the time grid, which `tune.py` does only to rank variants).

## Exclusion rule (R1, R2, R3): 72 variants, `tuning/exclusion.csv`

Floor 3, 2 or none; minimum threshold none, 0, 1, 2 sigma; unmeasured reference may exclude or not;
timing tolerance 0, 0.5, 1 duration.

| Variant | A2 one-TOI | A2 multi-TOI | B | B + Bx | median ln P excluded, A1 / B (SNR >= 8) |
|---|---|---|---|---|---|
| round 1 | 2 of 76 | 5 of 33 | 4 of 515 | 13 of 1547 | 0.729 / 0.676 |
| floor 3 sigma only | 2 of 76 | 5 of 33 | 1 of 515 | 4 of 1547 | 0.729 / 0.674 |
| unmeasured reference does not exclude, only | 0 of 76 | 4 of 33 | 1 of 515 | 3 of 1547 | 0.729 / 0.676 |
| **E: floor 3 sigma, unmeasured reference does not exclude, tolerance 0.5 duration (1 duration on multi-candidate stars)** | 0 of 76 | 2 of 33 | 0 of 515 | 1 of 1547 | 0.710 / 0.666 |

- On one-TOI stars the round-1 rule already met 5% overall; round 1 failed on the multi-TOI stars
  and on small SNR bins.
- All but two of the round-1 false exclusions on one-planet stars had an unmeasured reference or a
  brightening at the epoch. E removes them.
- The one false exclusion left under E (TIC 230130556, Bx) is a harness effect: the reference
  transit fell on a dip of the host's own (reference 9969 ppm, the injected planet's other
  transits 2744 ppm).
- Cost of E: a candidate whose reference the pipeline cannot measure (SNR below 5) gets no
  exclusions at all. The ln P criterion is judged at SNR 8 or more, so it does not see this cost.
- The minimum-threshold part (R1 sensitivity) added nothing once the floor was in; left off.
- Multi-TOI stars: 2 of 33 remain under E (TOI-216, TOI-7510). Reported, not judged.

## Event vetting (R4): 1080 variants, `tuning/vetting.csv`, and six full rules, `tuning/final.csv`

| Rule | false event, A1 one-TOI | false event, B null runs (hosts) | repeat recovery A2 / B | wrong period A2 one-TOI | wrong period B |
|---|---|---|---|---|---|
| round 1 | 26 of 76 | 142 of 515 (67 of 129) | 63 of 76 / 131 of 188 | 8 of 54 | 84 of 201 |
| E + best of the approved vetoes (depth and duration ratio within 1.5) | 11 of 76 | 51 of 515 (32) | 60 of 76 / 115 of 188 | 3 of 54 | 29 of 142 |
| E + event must be clean | 8 of 76 | 83 of 515 (48) | 63 of 76 / 127 of 188 | 4 of 57 | 49 of 162 |
| E + clean + depth ratio within 1.5 | 6 of 76 | 53 of 515 (34) | 62 of 76 / 121 of 188 | 2 of 56 | 32 of 144 |
| E + clean + depth ratio 1.5 + SNR 8 + anti-ratio 2 | 1 of 76 | 39 of 515 (26) | 60 of 76 / 113 of 188 | 1 of 53 | 22 of 127 |

**Said before the freeze, as the pre-registration requires: tuning cannot reach the event criteria
as approved.**

1. *The approved vetoes are not enough.* No combination of them gets below 11 of 76 (14%) false
   events on confirmed-planet hosts while keeping repeat recovery at 60%. The target is 5%.
2. *One veto that is not on the approved list does most of the work:* requiring the second event
   itself to be clean (at least 80% covered and 0.5 d from a gap or edge, the pipeline's own
   definition). Half of the false events are near an edge or gap; real repeats almost never are.
   Using it needs a new version of the pre-registration.
3. *Arm B cannot meet the event criteria under any vetting, because of its hosts.* The hosts are
   stars that carry real single- and duo-transit candidates. Under the two strictest rules, 31 of
   32 and 21 of 22 wrong arm-B calls used an event that is also found on that host with nothing
   injected: the host's own signal (a real repeat of its candidate, an eclipse, or variability).
   The null-run false-event rate stays at 8% for the same reason. On confirmed-planet hosts the
   same rules give 2 of 56 and 1 of 53 wrong calls.
4. *The recovery guard is tight.* The strictest rule recovers 113 of 188 in arm B (60.1%); the
   guard is 60%. "Clean + depth ratio 1.5" recovers 121 of 188 (64%).
5. *A gap in the approved text:* "each host is first run once with nothing injected" does not say
   which depth and duration the search uses. Here it is run once per draw, with that draw's values.

## Probability (R5), under E

One-TOI A2 + B + Bx trials, by reference SNR. p = (k + 0.5) / (n + 1).

| Reference SNR | false exclusions | trials | p | 95% upper bound |
|---|---|---|---|---|
| under 5 or unmeasured | 0 | 507 | 0.0010 | 0.0059 |
| 5-8 | 0 | 211 | 0.0024 | 0.0141 |
| 8-12 | 0 | 202 | 0.0025 | 0.0147 |
| 12-20 | 0 | 192 | 0.0026 | 0.0155 |
| over 20 | 1 | 511 | 0.0029 | 0.0092 |

Counts are too small to split by cadence as well. With rates this low the test set (about 700
trials) is predicted to hold about 2 false exclusions.

## Proposed for Jon's decision (a version 1.1 of the pre-registration; not applied)

- (a) Add "the event must be clean" to the vetoes R4 may use.
- (b) Judge arm B's event criteria (false-event rate, wrong-period rate) with each host's own
  null-run events removed, and report the all-trials numbers beside them. Arm A stays as approved.
- (c) State that the null run is made once per draw, with the draw's depth and duration.

## Not done yet

- The frozen harness (exact alias tables under the new rule) and its cross-check by a fresh agent
  on 20 tuning stars.
- The freeze itself (step 2), which needs Jon's approval.
