# Calibration (phase 3 of `docs/STRATEGY_living_alias_catalogue.md`)

| Folder | What | State |
|---|---|---|
| `published_lists/` | Do published alias lists keep the periods that were later solved? | Done for Hawthorn et al. 2024 (exploratory) |
| `mono/` | Injection and confirmed-planet calibration for single transits | **Run 2026-11-22 under pre-registration 1.1. Outcome: fail, confirmed by an independent re-derivation.** See `mono/RESULTS.md` and `mono/INDEPENDENT_CHECK.md`. |

## Published lists: result (2026-11-22)

`published_lists/hawthorn_solved.py` compares Tables A2 and A3 of Hawthorn et al. (2024) with the 26
candidates whose periods have since been solved (`data/catalog/status.csv`). It was not pre-registered:
the solved periods had been seen before the script was written.

| Question | Answer |
|---|---|
| Is the solved period at or above the paper's shortest allowed period? | 26 of 26 |
| Is it certainly in the published list (list has no gaps)? | 22 of 26. The other 4 (TIC 77437543, 118339710, 188620407, 404518509) have gaps, and the full online list is needed. |
| Is it inside the quoted errors of the paper's most probable period? | 17 of 26 (65%). An 18th, TIC 349091983, misses by 0.0001 d of rounding. |
| Is it the shortest allowed period? | 2 of 26. Median rank is 5th shortest of a median 31 aliases. |

**What it means.** The weakness found on 2026-11-01, where a depth-based epoch rule threw out the true
period for 46 of 83 confirmed planets, was in this project's own rule. Hawthorn's lists, which exclude
on coverage, do not show it. So the epoch-rule result is not a finding about published catalogues. It
stays what it was: a reason our pipeline needed the half-depth rule.

**Period priors.** With weights proportional to P^-a over each candidate's allowed aliases:

| a | Sum of ln p(true alias) | Mean PIT | KS p against uniform |
|---|---|---|---|
| 8/3 | −77.31 | 0.57 | 0.19 |
| 5/3 | −75.55 | 0.50 | 0.77 |
| 2/3 | −77.55 | 0.39 | 0.10 |
| flat | −83.39 | 0.29 | 0.003 |

A flat prior is rejected. The data prefer a = 5/3 but cannot separate it from 8/3 or 2/3 (1.8 to 2.0 in
log-likelihood, 26 objects). The solved sample is biased toward short periods, because those repeat in
TESS data and get solved first, so the true exponent may be shallower than this suggests. The value 8/3
for MonoTools is from memory; check Osborn et al. 2022 before quoting it.

**Not done.** The same test on the 12 published periods used in arXiv 2611.03019, and on the full online
alias lists for the 4 gapped candidates.

## Single-transit calibration: state (2026-11-22)

| File | What |
|---|---|
| `mono/prereg_mono_calibration.yaml`, `mono/APPROVALS.md` | Version 1.0, approved by Jon and frozen (hash checked by `mono/tests/test_prereg_mono.py`). |
| `mono/PROPOSED_v1.1.md` | Six gaps in 1.0 found while building the harness, with a proposed fix for each. **Needs Jon's approval before any trial is scored.** |
| `mono/monocal.py` | One trial: detrend and search as `st2_pipeline` does, then score against the truth. Defines the single-transit excluded set. |
| `mono/run_arm_b.py` | Arm B driver: injects a batman planet into the raw sectors before detrending. `--timing-only` discards outcomes. |
| `mono/tests/` | `pytest calibration/mono/tests`: synthetic light curves only. |
| `mono/results/pilot/arm_b_timing.csv` | The timing pilot. Seconds only. |

The pipeline change is in `salinas-second-transits/st2_pipeline.py`: `search_sector` now calls a new
`sector_events`, so the harness can change the flux between loading and detrending. No threshold moved.

**Timing pilot** (allowed by the pre-registration before the sample is fixed). 10 frozen hosts
whose raw light curves were already cached, 40 injected trials, outcomes discarded in memory:
10 s to load a star and 3.2 s per trial on average (longest 11 s, for stars with
about 40 sectors). The full run (600 arm-B trials and up to 240 arm-A trials) is therefore about an
hour of computing on one core. The pilot did not time downloads: all ten stars were cached. Most
of the 150 arm-B hosts and all arm-A stars are not, so MAST downloads will dominate. No cut to the
sample is needed.

**Not done.** The arm-A driver (confirmed planets as fake singles), the draw of the 150 hosts (needs
a TESS magnitude for all 1301 frozen hosts), and any scored trial.

## Round 2 (drafted 2026-11-22, not approved, nothing run)

`mono/prereg_round2.yaml`. Round 1's diagnostics were made after the outcome, so the fixes they suggest
are hypotheses. Round 2 tunes a repaired rule on stars already seen (a floor on the exclusion rule, a
timing tolerance, automatic vetting of second events, a false-exclusion probability per excluded
period), freezes it with Jon's approval, and tests it once on about 146 confirmed planets and 150
hosts that no earlier run has touched. Two harnesses must agree on 20 tuning stars before the freeze.
