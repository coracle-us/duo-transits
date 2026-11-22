# Independent period re-derivation: comparison with round 3, 2026-11-22

The CLAUDE.md rules (§5) require every claimed period to be re-derived independently before outreach.
This folder holds that re-derivation for the four periodic Salinas stars.

## How the blind run was set up

- It ran in a separate Claude Code session, in a new folder outside this repository
  (`~/Documents/GitHub/coracle/salinas-blind-periods`).
- It was given only the TIC numbers and the paper's single-transit mid-time, depth and duration.
- It was told not to read anything in this repository or its history, the Claude relay folder, ExoFOP,
  or any period estimate. It reports that it didn't, and came across no period.
- It downloaded its own data from MAST and used its own method:
  - a batman shape fit;
  - Δχ²-profile timing errors;
  - its own TESScut photometry, re-barycentred to the target;
  - an alias scan on raw flux.
- The data (591 MB) stayed in that folder. Its `RESULTS.md`, `scripts/` and `out/` are copied here
  unchanged.

## Periods

Round 3 values are from `../oddeven_summary.csv`; blind values are from `RESULTS.md`. The difference is
in units of the two errors added in quadrature.

| TIC | Round 3 P (d) | Blind P (d) | ΔP | Round 3 T0 (BTJD) | Blind T0 (BTJD) | ΔT0 | Transits (round 3 / blind) |
|---|---|---|---|---|---|---|---|
| 275647095 | 17.064019 ± 0.000032 | 17.064042 ± 0.000028 | 0.5σ | 2315.5125 ± 0.0018 | 2315.5125 ± 0.0020 | 0.0σ | 7 / 10 |
| 241550976 | 21.196941 ± 0.000038 | 21.196873 ± 0.000018 | 1.6σ | 2339.4442 ± 0.0019 | 2339.4451 ± 0.0011 | 0.4σ | 5 / 7 |
| 107113048 | 49.701517 ± 0.000168 | 49.701483 ± 0.000129 | 0.2σ | 2312.9429 ± 0.0033 | 2312.9433 ± 0.0030 | 0.1σ | 4 / 4 |
| 230133616 | 23.752973 ± 0.000022 | 23.752987 ± 0.000023 | 0.4σ | 2396.1280 ± 0.0007 | 2396.1275 ± 0.0007 | 0.6σ | 31 / 37 |

All four periods agree. The blind run excludes every shorter alias (P/2 and below, down to 1 d) at 15σ
or more, and every longer one through transits at epochs with no common factor. Its weakest exclusion
is P/2 for 107113048, which rests on two clean flat epochs (S10 1592.27 and S100 4077.35).

## Differences worth knowing

- **241550976 (the 1.6σ gap):** the blind run times two transits round 3 dropped.
  - Its TESScut light curve covers S100 E83 fully; ours was 64% covered and 0.03 d from an edge.
  - It also times S102 E85, which ours had 0.42 d from an edge.
  - The differences come from the two runs' quality masks and photometry.
  - Its fit has χ²ᵣ 1.02 against our 3.2, so its period is the better-constrained one.
- **275647095:** the blind run also times three partial transits that round 3's clean-transit rule
  excludes: S10 E−43, S37 E1 (half-covered at the sector end) and S99 E103.
  - Using only full transits, it gets P = 17.064035 ± 0.000027, 0.4σ from round 3.
  - The partial S37 E1 transit may be why the paper listed the star as a single.
- **230133616:** the blind run finds 37 transits, against our 31. It adds QLP times where TESS-SPOC's
  cadence mask drops the transit, and partials.

## Concerns both runs reached on their own

| | Round 3 | Blind |
|---|---|---|
| 275647095 odd/even | 1.3σ | 1.9σ |
| 107113048 odd/even | 1.0σ | 2.2σ (2 vs 2 transits; the odd ones are both QLP) |
| 241550976, 230133616 odd/even | 0.7σ, 0.6σ | 0.55σ, 0.2σ |
| Secondary at phase 0.5 | none (3σ limits 44–705 ppm) | none |
| 275647095 transit too long for the TIC star | ρ/ρ_TIC = 0.22 | 0.22 ρ☉ vs TIC 0.63 ρ☉ |
| 230133616 density | 1.65 × TIC | ~3 × TIC (b poorly constrained) |
| Depth scatter between sectors | above errors; TESScut diluted | above errors; source-dependent dilution |

The blind run also asks why 230133616 was listed as a single transit at all, with 37 transits in
S16–S86.

## Later follow-up on 275647095 (blind run only)

After this comparison was written, the blind run checked the long transit of 275647095 with Gaia DR3 and its
own difference images (scripts `s11`–`s14`, section "Follow-up" in `RESULTS.md`):
- R★ = 1.14 ± 0.06 R☉ from the parallax, as in the TIC. A bigger star does not explain the duration.
- The dip comes from within about 2″ of the target in 5 of 7 transits. Every Gaia neighbour, including
  the G = 14.3 star at 19.6″ behind the TIC contamination ratio of 0.21, is ruled out.
- What is left: an eccentric orbit (e ≳ 0.34), or an eclipsing system blended within ~2″ that Gaia does
  not resolve. Round 3's pixel check also found the signal on target.
