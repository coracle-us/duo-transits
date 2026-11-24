# No Third Body in Five Years of TESS Eclipse Timing of the δ Scuti Binary HIP 7666

Jon McBee, Independent researcher, USA

## Abstract

HIP 7666 (TIC 308953703) is a detached 2.37 d eclipsing binary with a δ Scuti primary. I timed 53 eclipses in Sectors 18, 58 and 85 (2019–2024) against an
empirical template, with the pulsations subtracted first. A linear ephemeris fits
(χ² = 56.7 for 50 dof). A quadratic term and a light-travel-time orbit are both
disfavoured by the Bayesian information criterion. The TESS period agrees with the 2000–2022
ephemeris to 0.020 ± 0.016 s. Injection–recovery excludes light-travel-time semi-amplitudes
above about 9 s for tertiary periods of 12–100 d and 80–210 s for 1,000–3,000 d. Away from
blind spots at one year and its harmonics, this excludes tertiaries above about 0.2 M☉. A blind
re-analysis reproduces these results. Secondary eclipses arrive 7–13 s late, consistent with the
9.4 s Rømer delay across the binary. Pulsation
phases vary with orbital phase beyond light travel time, so they cannot serve as a clock here.

## 1. Introduction

Eclipsing binaries with a δ Scuti component carry two clocks: the eclipses and the pulsations.
A third body moves both together through the light-travel-time effect (LTTE; Irwin 1952;
Murphy et al. 2014).

HIP 7666 was found to be an eclipsing δ Scuti system by Escolà-Sirisi et al. (2005). It was
studied by Feng et al. (2021) and Pawar et al. (2024), and solved to about 2% in mass by Lee et
al. (2026), who give M_A + M_B = 3.068 M☉ and a circular orbit. Their timing analysis used TESS
Sector 58 with historical minima and a linear ephemeris. None of these papers tests for a third
body. TESS has observed the star in three sectors; this note uses all three.

## 2. Method

**Light curves**
- S18: TESS-SPOC 1800 s (Caldwell et al. 2020).
- S58: SPOC 120 s (Jenkins et al. 2016).
- S85: TESS-SPOC 200 s.
- All three use PDCSAP flux. QLP (Huang et al. 2020) is a check.

**Models**
- An empirical binary template comes from the S58 light curve, folded on the Lee et al.
  ephemeris and fitted with a cubic spline.
- Fourteen of Lee et al.'s fifteen frequencies (all but f_orb) are fitted to the out-of-eclipse
  residuals of each sector. The result is subtracted before the eclipses are timed.

**Eclipse times**
- Each eclipse is fitted with a time shift, a depth scale and a linear baseline, using the
  template averaged over the exposure time.
- Errors come from a block bootstrap of the residuals.
- Within each sector, errors are scaled up where the scatter about a line requires it, by
  factors of 1.0–1.5.

**Pulsation delays**
- The delays come from the phases of the three strongest modes, and from a common shift of the
  whole pulsation model, in one-orbit segments.

**Model choice and limit**
- The rules were fixed and committed before any timing.
- The models are linear, quadratic and linear plus a circular LTTE orbit (period grid
  10–3,000 d). They are compared by BIC.
- The limit comes from injecting circular LTTE signals at 1,000 random phases into the real
  times. It is the amplitude that makes the linear fit fail at 3σ in 95% of cases.
- A synthetic test recovers injected signals in both clocks, with the correct sign.
- A separate analysis, given only the rules, the raw light curves and Lee et al. (2026), and
  written from scratch without sight of this one, re-timed the eclipses and recomputed the
  models and the limit. Its results are compared in Section 3.

## 3. Results

**Eclipse times**
- Median errors per eclipse are 3–7 s in S58 and S85 and 28–45 s in S18.
- The linear ephemeris is T0 = BJD_TDB 2460131.645322(11), P = 2.37222211(6) d.
- The alternatives are worse by BIC: a quadratic by +1.7, and the best LTTE orbit
  (P3 = 63 d, 4.8 s) by +6.7 (Figure 1).
- The primary and secondary periods agree (2.37222211(7) and 2.37222212(11) d), so there is no
  apsidal motion.
- Our S58 times match Lee et al.'s Table 1 to −1.4 s on average.
- The blind re-analysis gives P = 2.372222090(68) d (0.2σ from ours), M0 χ² = 55.6 for 50 dof,
  and ΔBIC +1.7 (quadratic) and +6.0 (LTTE): the same choice of the linear model.

**Secondary-eclipse delay.** Secondary eclipses come +12.7 ± 1.9 s after phase 0.5 with our
template, and +6.8 ± 2.0 s in the blind re-analysis. The difference (2.2σ) comes from method,
not noise:
- the S58 eclipse profiles have a small odd part (0.2–0.3 ppt);
- the blind analysis fitted a symmetrized template, which moves every primary by +2.5 s and
  every secondary by −2.7 s;
- doing the same in our code reproduces the shift.

Both values agree with the Rømer delay across the inner orbit, P(K_B − K_A)/(πc) = +9.4 s
(Kaplan 2010), within 1.7σ. Lee et al. suggested pulsations as the cause.

**Limit.** For circular LTTE orbits, the 95% upper limits on the semi-amplitude are:

| P3 | Limit | Tertiary mass excluded above (sin i3 = 1) |
|---|---|---|
| 12–100 d | about 9 s (median) | about 0.21 M☉ |
| 100–1,000 d | about 40 s | about 0.19 M☉ |
| 1,000–3,000 d | 80–210 s | about 0.17 M☉ |

The overall median is 0.19 M☉ (blind: 0.20 M☉), and 94% of grid periods are below 0.3 M☉.

- With QLP for S18, whose light curve is cleaner than TESS-SPOC's, the limit at 2,000–3,000 d
  falls to about 0.08 M☉.
- The three sectors fall in the same season, so the limit is weak near one year and its
  harmonics:
  - near 364 d, only m3 < 1.8 M☉ is excluded;
  - near 182 d, m3 < 0.75 M☉ (blind: 0.87 M☉);
  - near 121 d, m3 < 0.5 M☉.
- Periods below about 12 d would be dynamically unstable (Mardling & Aarseth 2001).

**Pulsation delays**
- Per-orbit errors on the common delay are 13 s in S58 and 17 s in S85.
- They scatter more than these errors (rms 22–24 s), several times the eclipse scatter.
- In 0.3 d segments the delays vary with orbital phase with a cosine amplitude of +45 ± 21 s.
  Light travel time of the primary would give −10 s.
- Many of the frequencies are f_i ± n f_orb combinations, which points to tidally affected
  modes. The delays are therefore shown but not used for the limit.

**A QLP data issue.** For the same cadences, the QLP Sector 85 timestamps for this star are 103–122 s
earlier than TESS-SPOC's after BTJD 3623.448. This shifts the QLP eclipse times by the same
amount. QLP Sectors 18 and 58 agree with the SPOC products.

## 4. Conclusion

Five years of TESS eclipse times show no third body around HIP 7666 above about 0.2 M☉
(median 0.19 M☉) for periods of two weeks to eight years, apart from gaps near one year and its
harmonics. An independent blind re-analysis confirms the period, the model choice and the limit. Further TESS sectors and
the historical minima would extend the baseline.

## Figure 1 caption

Top: eclipse O−C against the linear ephemeris (blue; filled primary, open secondary) and
pulsation delays per orbit with each sector's linear trend removed (red), for Sectors 18, 58
and 85. Bottom: 95% upper limit on the LTTE semi-amplitude against tertiary period (solid:
main products; dashed: S18 from QLP), labelled with the minimum excluded tertiary mass in M☉.
Plotting code was written with Claude (Anthropic).

## Acknowledgments

This paper includes data collected by the TESS mission, obtained from MAST at STScI.
- SPOC 120 s: doi:10.17909/t9-nmc8-f686.
- TESS-SPOC: doi:10.17909/t9-wpz1-8s54.
- QLP: doi:10.17909/t9-r086-e880.

Software: lightkurve, astropy, numpy, scipy, matplotlib, tess-point.

AI statement: I used Anthropic's Claude, an AI assistant, under my supervision to write and run
the analysis code (including injection tests) and to draft this note. I checked the numbers
against the saved outputs and take full responsibility for the content.

## References

- Caldwell, D. A., et al. 2020, RNAAS, 4, 201
- Escolà-Sirisi, E., et al. 2005, A&A, 434, 1063
- Feng, G., et al. 2021, MNRAS, 508, 529
- Huang, C. X., et al. 2020, RNAAS, 4, 204
- Irwin, J. B. 1952, ApJ, 116, 211
- Jenkins, J. M., et al. 2016, Proc. SPIE, 9913, 99133E
- Kaplan, D. L. 2010, ApJL, 717, L108
- Lee, J. W., et al. 2026, arXiv:2611.17590
- Mardling, R. A., & Aarseth, S. J. 2001, MNRAS, 321, 398
- Murphy, S. J., et al. 2014, MNRAS, 441, 2515
- Pawar, T., et al. 2024, arXiv:2410.00590
- Ricker, G. R., et al. 2015, JATIS, 1, 014003
