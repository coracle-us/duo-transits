# Single-transit calibration: result (run 2026-11-22, pre-registration version 1.1)

**Outcome under the pre-registered criteria: fail.** The full table is `results/score_summary.md`;
every rate is in `results/score_rates.csv`, from `score.py` (committed before the run finished, bb965eb).

The two halves of the result point in opposite directions.

> **Corrected 2026-11-22 after the independent check (`INDEPENDENT_CHECK.md`).** The tables below are
> the first run. With three harness errors fixed (`results_corrected/`): A2 false exclusion is
> **7 of 109 (6.4%), which fails the overall limit too**; A2 wrong period 15 of 74; B wrong period
> 84 of 201; B false exclusion unchanged at 4 of 515. The independent run gives 7 of 109, 16 of 74,
> 83 of 200 and 4 of 516. "No independent re-derivation has been run" under Limits no longer holds.

## 1. Period exclusion from one transit: reliable overall, but two bins miss

| Arm | True period wrongly excluded | 95% upper bound | Criterion (≤ 5%) |
|---|---|---|---|
| A2, confirmed planets, transits visible | 4 of 109 (3.7%) | 8.2% | overall pass |
| A2, SNR 8-12 | 2 of 17 (11.8%) | 32.6% | **fail** |
| A2, SNR above 20 | 2 of 38 (5.3%) | 15.7% | **fail** |
| B, injections before detrending | 4 of 515 (0.8%) | 1.8% | pass, every bin |
| B, trials where the true period had a clean epoch | 4 of 189 (2.1%) | 4.8% | (not a criterion) |

The excluded set is also useful: the median trial excludes 61% (arm B) to 67% (arm A) of ln P between
1 and 1000 d, falling to 20% below reference SNR 6.

## 2. Automatic "second transit" calls: not reliable

| Arm | Calls (resolved or narrowed) whose alias set lacks the true period | Criterion (≤ 3%) |
|---|---|---|
| A1, other transits removed | 47 of 47 (every call is false by construction) | **fail** |
| A2 | 16 of 75 (21%) | **fail** |
| B | 65 of 178 (37%) | **fail** |

The useful A1 number is how often a false event passes the cuts at all: 51 of 109 stars (47%).
Real repeats are found in 81% of A2 stars and 70% of arm-B trials that have one in clean data.

## Diagnostics made after seeing the outcome (post-hoc, not pre-registered)

- **All 4 A2 false exclusions are on stars with more than one TOI** (TOI-7510, 2525, 216, 261). On
  the 76 single-TOI stars there are 0 of 76. Two of those systems are, from memory, known for large
  transit-timing variations (TOI-216, TOI-2525); that has not been checked against the literature.
  A linear ephemeris is the calibration's "truth", so a planet that does not follow one looks like
  a false exclusion.
- **Wrong periods come from the event, not from the exclusion rule.** In arm B, 64 of the 65 wrong
  calls have a second event that is not on the true ephemeris at all; only 1 is an alias the rule
  excluded. In 49 of them no true repeat was detected and something else passed the cuts.
- **Arm B's hosts are not clean.** They are candidate hosts, and only the catalogued transit times
  were masked. A host with a real periodic signal, an eclipsing binary, or a second dip keeps those
  in its light curve, and they pass as repeats. This inflates arm B's wrong-period rate by an
  unknown amount and is a flaw in the arm's design.
- **That does not explain arm A.** On single-TOI confirmed-planet hosts with the planet's transits
  removed, 26 of 76 (34%) still give a passing event. Transit wings left by timing variations, other
  planets without a TOI, and systematics are all possible. Not investigated.

## What this means for the catalogue

- A period excluded from a single transit by the half-depth rule can be published, with the measured
  false-exclusion rate attached. Multi-planet and timing-variation systems need a wider tolerance
  or a flag.
- A second transit found by the box search alone must not be turned into an alias set without
  vetting. The Salinas work already vetted each event by hand; the calibration shows that step is
  necessary, and that the automatic cuts (SNR ≥ 7, depth and duration within a factor of 2) are
  not enough on their own.

## Limits

- No independent re-derivation has been run. CLAUDE.md section 5 requires one before any of this is
  claimed.
- 21 of the 150 arm-B hosts have no light curve on MAST (515 trials on 129 hosts). 11 of the 120
  arm-A stars have fewer than 3 clean transits (109 scored).
- The `implementation_notes` in the pre-registration were not individually approved.
- Arm B's draws use a circular orbit and one planet per trial.
- The pass criteria were not changed after the run, and nothing was re-thresholded.
