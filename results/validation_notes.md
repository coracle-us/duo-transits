# Alias pipeline validation (plan step 4), 2026-10-11

Pipeline: `scripts/alias_pipeline.py`. Validation: `scripts/validate_pipeline.py` →
`results/validation_known.csv` (A), `results/validation_paper.csv` (B). Per-epoch and
per-alias tables in `results/aliases/{all,maxs39}/` are gitignored; regenerate with
`python scripts/validate_pipeline.py` (~1 min).

## A. Known periods, all TESS data (26 candidates)

- **26/26 PASS**: the true alias is never excluded.
- **19/26 uniquely recovered** from TESS alone (only the true alias survives).
- The 7 non-unique cases keep the true alias among the survivors. In these, the extra transits
  fall on the shared Tc1 + m·dT grid or in gaps, or are too weak to count as detections:
  - TOI-2589: 17 survivors
  - TOI-6695: 26 survivors (TTVs)
  - TOI-2490: 10 survivors
  - TOI-4310: 13 survivors
  - TIC 439491923: 10 survivors
  - TIC 242241304: 29 survivors
  - TOI-1861: 12 survivors (resolved only with ground-based data)
- **By tier, all pass:**
  - 9 peer-reviewed
  - 1 from the catalog paper
  - 5 from the thesis
  - 8 ExoFOP community claims (independently confirmed where unique: 81089255, 188620407, 265465927, 289840544, 304339227, 457649900)
  - 3 TOI/QLP periods

## B. Paper-era data (S ≤ 39) vs Hawthorn+24 Table A3 N_alias (85 candidates)

- **Alias counts against N_alias:** equal for 40, within +1 for 65, within +3 for 76. Only one candidate has fewer than the paper.
- **The one with fewer:** HD 21520 (TIC 404518509) has 1 survivor vs the paper's 18. The SPOC 2-min
  data show a third transit in S3 (BTJD 1406.14) that fixes n = 28 (25.13 d), which matches Nies+2024.
  The paper used FFI light curves only.
- **Large excesses, where we are more conservative:**
  - shallow signals: TIC 20904104 (0.9 ppt), 332697924 (1.2 ppt), 145006304 (2.4 ppt)
  - grazing transits (b ≳ 0.93): 107113345, 211409161, 215402824, 317923092

  In these cases a missing transit would not be clearly visible, so our depth-based test declines
  to exclude. MonoTools excludes on coverage alone.
- **Our light curves cannot see TOI-709.02 (TIC 396720998).** Neither known transit is recovered,
  giving 48 survivors vs 12. Look into this before working on the candidate.
- **Known transits:** 166/170 are recovered as "supported".
- **Paper table inconsistency:** in Table A3, P_min = dT/(N_alias − 1) for most rows. Our survivor sets run from n = 1 to N_alias, so N_alias is the count and P_min is off by one.

## Decisions made during validation (post hoc; disclose)

Thresholds were set before the run except for two fixes, each prompted by a failure:

1. **Detected-transit cut (prompted by HD 21520).** A dip can exclude *other* aliases only if:
   - SNR ≥ 8
   - depth between 0.5 and 2 times the expected depth
   - coverage ≥ 80%
   - baseline on both sides

   A gap-edge systematic (SNR 7.3, one-sided baseline) had wrongly killed the true alias.
2. **Clean-window rule for exclusion (prompted by TIC 265465927 and TOI-6695).** No exclusion
   from windows with red-noise factor β > 3 or a significant brightening. Both false exclusions
   came from such windows: an edge ramp with β ≈ 9.5, and a β ≈ 8.5 window.

Part B results did not change after fix 2. Both fixes make the pipeline more conservative,
never less. Thresholds are now frozen. Any further change should be re-validated on all of
A and B.

## Caveats

- TTVs are not modelled; for TOI-6695 the true alias survives only as "allowed".
- Single-product light curves are used (SPOC 2-min > TESS-SPOC > QLP). The vetting step (6)
  still needs a second, independent reduction for any new detection.
