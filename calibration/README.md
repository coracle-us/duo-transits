# Calibration (phase 3 of `docs/STRATEGY_living_alias_catalogue.md`)

| Folder | What | State |
|---|---|---|
| `published_lists/` | Do published alias lists keep the periods that were later solved? | Done for Hawthorn et al. 2024 (exploratory) |
| `mono/` | Injection and confirmed-planet calibration for single transits | Pre-registration drafted; **not approved, nothing run** |

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
