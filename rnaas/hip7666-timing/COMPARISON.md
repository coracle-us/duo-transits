# Blind re-derivation vs main analysis (2026-11-24)

**Setup.** A separate agent got only `RULES.md`, the six raw TESS FITS files, an empty venv and
the public Lee et al. paper. It worked in a separate folder and was told not to read this repo.
It wrote its code from scratch; the code and outputs are now in `independent/`, its report in
`independent/REPORT.md`.

**Re-run.** The main limit was re-run on a 600-point P3 grid after the comparison (change log
item 7). The main numbers below come from that run.

## Comparison

| Quantity | Main | Blind | Difference | Beyond 1σ? |
|---|---|---|---|---|
| Eclipses used (S18 / S58 / S85, p + s) | 8+9 / 10+12 / 7+7 = 53 | same 53, same 3 dropped | — | no |
| P (M0) | 2.37222211(6) d | 2.372222090(68) d | 2.0e-8 d (0.2σ) | no |
| P − P(Lee Eq. 1) | +0.020 ± 0.016 s¹ | +0.018 ± 0.006 s¹ | 0.002 s | no |
| M0 χ² / dof | 56.7 / 50 | 55.6 / 50 | 1.0² | no |
| BIC M0 | 68.6 | 67.6 | 1.0 | no |
| ΔBIC M1 (quadratic) | +1.7 | +1.7 | 0.0 | no |
| ΔBIC M2 (LTTE) | +6.7 (best P3 63 d) | +6.0 (best P3 10.2 d; 63 d nearly as good) | 0.7 | no, same choice (M0) |
| P primary only / secondary only | 2.37222211(7) / 2.37222212(11) d | 2.372222072(87) / 2.372222120(111) d | 0.3σ / 0.0σ | no |
| Median m3,min (P3 ≥ 12 d, excluding 330–400 d) | 0.191 M☉ | 0.199 M☉ | 4% | n/a³ |
| … 12–100 / 100–1,000 / 1,000–3,000 d | 0.210 / 0.193 / 0.169 | 0.225 / 0.199 / 0.178 | 3–7% | n/a³ |
| Worst m3,min outside 330–400 d | 0.75 M☉ at 182.5 d | 0.87 M☉ at 182.5 d | same period | n/a³ |
| S18-from-QLP variant, median m3,min | 0.152 M☉ | 0.169 M☉ | 10% | n/a³ |
| **Secondary offset from phase 0.5** | **+12.7 ± 1.9 s** | **+6.8 ± 2.0 s** | **5.9 s (2.2σ)** | **yes** |
| QLP S85 timestamp step | −103 to −122 s after BTJD 3623.448 | the same, found on its own | — | no |

¹ The main error includes Lee's period error; the blind error does not.
² The χ² difference is far below its own scatter, √(2·50) = 10.
³ The limit is a threshold, not a measurement with an error. The differences come from small
differences in the error-scaling factors (for example S18 primaries 1.46 vs 1.56) and in the grid.

## The one disagreement: the secondary-eclipse offset

**The difference is systematic.** Eclipse by eclipse, main minus blind is −2.5 s for every
primary (scatter 0.1–0.3 s in S58/S85) and +2.7 to +3.0 s for every secondary (scatter
0.2–0.8 s). Noise would not do that.

**Cause:** the template.
- The main analysis fits the S58 template as it is and measures times from that template's
  symmetry axis.
- The blind analysis fitted a symmetrized template.
- The S58 eclipse profiles have an odd (asymmetric) part of 0.27 ppt (primary) and 0.20 ppt
  (secondary), within |phase| < 0.044.

**Test.** Refitting the S58 eclipses in the main code with a symmetrized template moves the
primaries by +2.35 ± 0.09 s and the secondaries by −2.72 ± 0.22 s, the same as the main−blind
differences. So the offset depends on the method at about ±3 s per eclipse type.

**What doesn't change:**
- The period, χ², BIC and limit: a constant shift per eclipse type is absorbed by the free
  secondary-offset parameter.
- Both offsets agree with the 9.4 s Rømer delay, at 1.7σ (main) and 1.3σ (blind).

**Effect on the draft.** It now gives the offset as 7–13 s and drops the e cos ω estimate that
relied on the main value alone.

## Corrections the comparison prompted in the main analysis

The main limit had been computed on a 120-point P3 grid, which missed the narrow blind spots at
half and a third of a year. On a 600-point grid:
- the median is 0.19 M☉ (was 0.18);
- the worst case outside 330–400 d is 0.75 M☉ at 182.5 d (was 0.45 M☉);
- near 121 d it is 0.50 M☉.

The draft has been updated to match.
