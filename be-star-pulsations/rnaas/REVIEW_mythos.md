# Review of the Be-star RNAAS draft (`rnaas/ms.tex`, 2026-11-23)

Reviewer: Claude (Mythos), acting as a sceptical RNAAS moderator and Be-star specialist.
Written 2026-11-23. Nothing was edited, committed, sent or posted. Two new files only: this
review and `rnaas/review_mythos_checks.py` (my check script, read-only on `results/`).

## Verdict: submit after fixes

The low-frequency result is real as far as I can test it. It reproduces from the per-sector
table, survives sector fixed effects and dropping the 100 highest-leverage stars, and I found
no prior paper that did this test. But the draft should not go in as written, for five
reasons. (1) It calls the sample "Be stars"; the source catalogue calls them photometric
candidates. (2) The abstract leads with the A1 result, which is the weakest one and is not a
pulsation measurement. (3) The dilution paragraph applies a correction to the whole sample
that is only valid for the SPOC third of it. (4) The white-noise paragraph omits that its
ratio test fails in the SPOC subset, and misses the evidence that explains the control.
(5) "About half go the other way" reads measurement noise as stellar behaviour. All five are
text changes plus small reruns of cached analyses. None needs new data.

**How to read my numbers.** Sections marked *my check* come from `review_mythos_checks.py`
run on `results/full/sector_metrics.csv`. It reproduces the main β values exactly (A1 0.40,
B1 0.80, flicker 0.84, white 0.33; 1,084 stars, 4,108 sectors). Its uncertainties are a
500-draw bootstrap over stars, with **no permutation or swap nulls**. Treat them as
indications. Every one that goes into the note must be rerun through `stats4.run_test`.

---

## MUST FIX

### M1. The sample is Be-star *candidates*, not Be stars

> "the NEOWISE W1 light curves of 1,084 Be stars from the catalogue of \citet{Jia2026}"
> "\citet{Jia2026} identified Be stars whose disks changed state"

Jia et al. select from 4.2 million Gaia "candidate OBA stars" by W1 variability and a
W1 vs W1−W2 correlation, then visual inspection. They call the result "BeT candidates" and
note "less direct spectroscopic validation" than their 2024 catalogue. A summary of their
§4.1 (I could not read that section verbatim) gives 654 of 1,772 SIMBAD-classified sources
with a Be-type classification and 112 (6%) classed as YSOs, LPVs, AGB or WR stars.

A Be specialist will ask how many of the 1,084 are Herbig Ae/Be stars or other
variable-disk objects, where "more IR excess, more optical flicker" has a different cause
(accretion, dust).

**Fix (text).** Use "Be-star candidates" in the title-adjacent sentences, abstract, Method
and caption. Add one sentence to Method:

> "The catalogue is photometric: stars are Gaia OBA candidates with correlated W1 and
> W1$-$W2 changes, and most lack spectroscopic confirmation."

**Fix (analysis, strongly advised).** Rerun the main test on the confirmed subset:
cross-match the 1,084 stars against Jian et al. 2024 (spectroscopically confirmed OBA
parents, VizieR J/A+A/682/A59) and BeSS, then `stats4.run_test` on that subset for B1 and
flicker. If SIMBAD/VizieR are still down, the Jia catalogue's own overlap flag with Jian et
al. 2024 (557 stars) will do. One clause reports it: "The confirmed subset (N stars) gives
β = … for B1."

### M2. The abstract leads with the weakest result, and says "power" for amplitude

> "The amplitude of the strongest peak rises by 4\% per 0.1\,mag of excess, and the
> stochastic power below 1\,c\,d$^{-1}$ rises by about 8\%. Both trends pass a null test
> that swaps disk histories between stars ($p<10^{-4}$)."

Three problems.

- **A1 is not confirmed by the body.** The Results say A1 is flat relative to the white
  noise (0.07, p 0.33) and not significant in SPOC (0.17, p 0.11). A reader who goes from
  abstract to Results sees an overstated abstract.
- **A1 is not the pulsation.** *My check:* only 470 of 1,084 stars keep the same strongest
  frequency (within 10%) across their sectors. In those stars A1 β = +0.18 ± 0.15
  (+0.09 ± 0.12 if f1 > 0.8 c/d). In the 614 stars where the strongest peak moves, A1
  β = +0.48 ± 0.12. The chance that the strongest peak lies below 0.8 c/d rises with the
  excess (+0.135 ± 0.043 per mag). So most of the A1 slope is the "strongest peak" becoming
  a low-frequency stochastic peak. The per-star sign test agrees: 51.6% of A1 slopes are
  positive (binomial p = 0.32), against 57% for B1 (p = 2×10⁻⁶).
- **"Power" is wrong.** B1 is an rms amplitude. Power would rise by 16%.

**Fix.** Lead with B1. Replacement abstract (about 146 words; the two numbers in the
second-to-last sentence are from my check and need the pipeline rerun in S1 below):

> Be stars build and lose their decretion disks over months to years, and NEOWISE tracked
> the disk emission of thousands of candidates. I compare TESS light curves with the NEOWISE
> W1 light curves of 1,084 Be-star candidates from the catalogue of \citet{Jia2026}, sector
> by sector and within each star. When the W1 disk excess is larger, the star is more
> variable in the optical. The variability amplitude below 1\,c\,d$^{-1}$ rises by about 8\%
> per 0.1\,mag of excess, and passes a null test that swaps disk histories between stars
> ($p<10^{-4}$). At higher frequencies the rise is about 3\%, and I find no separate response
> of the dominant periodic signal. Stars whose excess changed most carry the trend: of 177
> with changes above 0.4\,mag, 71\% show it. Stars with small changes are uninformative, so
> this is a population tendency, not a rule for each star.

In Results, define A1 honestly where it first appears: "the amplitude $A_1$ of the strongest
peak, which is a coherent pulsation in some stars and a stochastic peak in others".

### M3. The dilution correction does not apply to QLP, which is 70% of the sample

> "Light from the disk dilutes relative amplitudes. A star whose $T$ magnitude changes by
> $-k$ per mag of W1 excess has every $\beta$ raised by $0.921k$ once corrected. … which
> raises $\beta$ to about 0.6 for $A_1$ and 1.0 for $B_1$. The uncorrected values above are
> therefore conservative."

QLP does not normalise by the star's measured flux. It measures difference-image flux and
adds the flux expected from the TIC magnitude (Huang et al. 2020, both RNAAS notes; "adding
the expected flux from each source based on the TESS magnitude estimates from the TIC").
A QLP relative amplitude is ΔF / F_TIC. Extra disk light does not dilute it. SPOC PDCSAP is
normalised by the measured flux, so there it does.

Consequences:

- The corrected values "0.6 and 1.0" are wrong for the pooled sample. The correction belongs
  to the 329 SPOC/TESS-SPOC stars only.
- Applied there, it removes most of the QLP–SPOC gap. SPOC β + 0.921k (k = 0.21, so +0.19)
  gives A1 0.36, B1 0.88, flicker 0.93, against QLP 0.52, 0.86, 0.90. *My check* of the raw
  gap: QLP − SPOC = +0.35 ± 0.22 (A1), +0.18 ± 0.22 (B1), +0.15 ± 0.26 (flicker). The
  predicted gap is +0.19. This is a better account of "SPOC does not confirm A1" than the
  draft's "cleaner photometry" framing, and it answers packet question 2.
- "Conservative" also assumes the variable light is photospheric. If the low-frequency
  signal comes from the inner disk (Labadie-Bartz et al. 2021 suggest part of the red noise
  may), there is nothing to correct.

**Fix (replacement paragraph, about the same length):**

> The two pipelines normalise differently. SPOC fluxes are relative to the star's measured
> brightness, so added disk light dilutes them; QLP adds the flux expected from the TIC
> magnitude to a difference flux \citep{Huang2020}, so it does not. From SAP photometry of
> 125 SPOC and TESS-SPOC stars the $T$ magnitude brightens by $k=0.21$\,mag (95\% range
> 0.11--0.30) per mag of W1 excess, which lowers SPOC $\beta$ by $0.921k=0.19$ if the
> variable light is photospheric. Corrected, the SPOC stars give 0.36 for $A_1$, 0.88 for
> $B_1$ and 0.93 for flicker, close to the QLP stars (0.52, 0.86, 0.90). Four stars, one of
> them (TIC~55295028) significantly, fade in the optical as W1 brightens, as shell stars
> should; removing them changes $\beta$ by less than 0.01.

Then change the second caveat (see M4) so it no longer says SPOC fails to confirm A1 without
this context. **Please verify the QLP statement against the cached files before using it:**
confirm `flux_col` for QLP sectors is the TIC-normalised `SAP_FLUX`, and state the SPOC
flux column (PDCSAP or SAP) in Method.

### M4. The white-noise paragraph omits the SPOC failure and misreads the control

> "The formal flux error of the SPOC and TESS-SPOC products does not ($\beta=0.09$,
> $p=0.56$), so this is not an instrumental noise floor. Relative to the white noise, $B_1$,
> red noise and flicker still rise ($\beta=0.41$--0.51, $p<10^{-3}$), but $A_1$ does not."

- **Omission.** In the SPOC subset the ratio tests do not survive the swap null
  (`control_tests.csv`: B1/white 0.27, p_swap 0.16; red 0.23, p 0.22; flicker 0.32, p 0.13).
  The draft says the SPOC subset "confirms the low-frequency trends". It confirms the raw
  slopes. It does not confirm them relative to the control. That has to be said.
- **The inference is too weak.** A flat formal error in 329 SPOC stars does not show the
  15–24 c/d rise is intrinsic in a sample that is 70% QLP.
- **There is direct evidence, and it is better.** *My check*, splitting stars by their own
  median white-noise level:

  | Tercile (median 15–24 c/d amplitude) | β white | β A1 | β B1 |
  |---|---|---|---|
  | quietest (0.026 ppt) | +0.61 ± 0.19 | +0.42 ± 0.24 | +0.96 ± 0.23 |
  | middle (0.062 ppt) | +0.49 ± 0.19 | +0.34 ± 0.13 | +0.87 ± 0.16 |
  | noisiest (0.122 ppt) | **0.00 ± 0.14** | +0.43 ± 0.13 | +0.63 ± 0.15 |

  The control is flat exactly where photon noise dominates the band, and rises where the
  star is bright enough for signal to show there. Also, once ln B1 is in the regression the
  white-noise slope on excess is −0.03: the 15–24 c/d rise is fully accounted for by its
  covariance with the low-frequency signal (within-star correlation 0.46). So the
  "control" is contaminated by the tail of the stochastic signal in bright stars. It is not
  a noise floor, and dividing by it over-corrects. In the noisiest third, where the control
  works, B1 and A1 still rise.

**Fix (analysis).** In `s4c_controls.py` add a variant that splits main-sample stars into
terciles of star-median `white_noise` and runs `stats4.run_test` (both nulls) on
white_noise, A1 and B1 in each. Add the two-regressor fit (ln white on excess + ln B1).

**Fix (replacement paragraph, if the rerun agrees):**

> The white-noise control also rises ($\beta=0.33$), but only in bright stars. In the third
> of stars with the highest 15--24\,c\,d$^{-1}$ level, where photon noise dominates that
> band, the control is flat ($\beta=0.00\pm0.14$) while $B_1$ still rises ($0.63\pm0.15$).
> In quieter stars the control tracks $B_1$, and its slope vanishes once $B_1$ is included
> in the fit. The formal flux error of the SPOC products is also flat ($\beta=0.09$,
> $p=0.56$). The band is therefore the tail of the stochastic signal, not a noise floor.
> Relative to it, $B_1$, red noise and flicker still rise in the full sample
> ($\beta=0.41$--0.51, $p<10^{-3}$) but not significantly in the SPOC stars alone
> ($0.23$--0.32, $p_{\rm swap}\geq0.13$).

### M5. "About half go the other way" reads noise as behaviour

> "Within individual stars the trend is weak and about half go the other way"
> "48\% of stars have a negative $A_1$ slope of their own"

Half the sample cannot show a slope of either sign. *My check:* the median star's excess
changes by 0.11 mag between its sectors, and 519 of 1,084 stars change by less than 0.1 mag.
Against a 0.54 rms scatter in ln A1 their slopes are coin flips. The sign fraction climbs
with the size of the disk change:

| W1 excess range | Stars | A1 positive | B1 positive | Flicker positive |
|---|---|---|---|---|
| < 0.05 mag | 374 | 47% | 49% | 50% |
| 0.05–0.1 | 145 | 45% | 55% | 45% |
| 0.1–0.2 | 190 | 53% | 58% | 57% |
| 0.2–0.4 | 198 | 52% | 60% | 61% |
| > 0.4 mag | 177 | 65% | 71% | 72% |

Related: 56 stars (5%) carry half of the regression weight (Σx²) and 275 carry 90%. That is
not a fault. *My check:* dropping the 10, 50 and 100 highest-weight stars leaves A1 at
0.48, 0.45, 0.43 and B1 at 0.88, 0.94, 1.03. But the note should say who carries the result.

**Fix.** Replace the first caveat with:

> "First, the effect is small against the scatter: at fixed disk excess, $\ln A_1$ varies by
> 0.54 rms between sectors of the same star, and the within-star rank correlations are
> 0.09--0.20. Half the stars change by less than 0.1\,mag between sectors and cannot show a
> trend. Of the 177 that change by more than 0.4\,mag, 71\% have a positive $B_1$ slope and
> 65\% a positive $A_1$ slope."

**Fix (analysis).** Add this table to `s4b_jackknife.py` (per-star slope sign by x-range
bin) and write it to `results/full/`. Also write the drop-top-N-leverage rows, which are a
stronger robustness statement than drop-one-star.

The Niwano comparison then needs the matching base rate. RX J0440.9+4431 has a 0.32 mag
range; in the 0.2–0.4 bin 48% of A1 slopes are negative, so the 0.16 probability stands.
But see S4.

---

## SHOULD FIX

### S1. The f/f1 bins have an edge exactly on the dominant frequency

> "$\beta$ is 0.73--1.16 below $0.7f_1$ and 0.35 at the dominant group itself, the same as
> at high frequency."

`freq_beta.csv` has relative bins 0.794–1.000 (β 0.53) and 1.000–1.259 (β 0.35). The
star's own peak sits on the boundary and is split between them. Quoting 0.35 picks the
lower half. It also drives the "Štefl frequency" remark in the literature notes.

**Fix (analysis).** In `s6_freq_beta.py` shift the relative grid by half a bin so one bin
is centred on f/f1 = 1 (0.891–1.122), rerun, and quote that bin. Given M2, also run it
only for the 470 stable-f1 stars, where f1 means something. If the centred value is about
0.45, the sentence "the same as at high frequency" must go.

> "The coherent pulsation therefore follows the broadband level, while the low-frequency
> stochastic signal rises about 2.5 times faster."

Replace with: "I find no separate response at the dominant frequency: it rises with the
broadband level, while the signal below 1 c/d rises about 2.5 times faster." The data do
not isolate a coherent mode.

### S2. Say that outbursts themselves are in B1 and flicker

A 5-d-knot spline leaves most of a days-long outburst brightening in the 0.25–0.8 c/d band
and in the 1-d-smoothed residuals. "Stochastic" (abstract, Results, Discussion) is an
interpretation. The measured thing is low-frequency amplitude, which includes the outburst
light curve. One clause in the Discussion covers it: "$B_1$ and flicker include the
days-long brightening of an outburst itself, so they count ejection events as well as
stochastic variability." This fits your "more ejection episodes per sector" reading.

### S3. Add sector fixed effects as a plain-language partner to the swap null

The swap-null sentence is clear enough for a moderator (packet question 5). A second,
simpler statement helps a reader who does not trust a bespoke null. *My check:* adding a
fixed effect for each TESS sector to the per-star fit changes nothing (A1 0.41, B1 0.83,
flicker 0.85, white 0.33). Restricting to changes within one star *and* one FFI cadence
(721 stars) gives B1 0.58 ± 0.23, flicker 0.48 ± 0.26, A1 0.22 ± 0.30: lower power, same
sign. Add: "Adding a fixed effect for each TESS sector leaves every $\beta$ unchanged."
Rerun in `s4_test.py` as a variant.

### S4. Niwano paragraph: two wording errors

> "With 48\% of single stars going that way, four of five is consistent with the positive
> population mean."

- "Single stars": nothing establishes that the sample is single. Say "sample stars".
- Niwano et al. measured coherent peaks after removing red noise, against J band. Your A1
  is not that (M2). Say the comparison is loose: "Their statistic is not my $A_1$, but
  with 48% of sample stars having a negative $A_1$ slope, four of five is unremarkable."
- The X-ray-binary cross-match used a list "from memory" with SIMBAD down. Either redo it
  or keep the claim to the one Niwano star, as the draft does. Do not add a count.

### S5. State the sector–visit geometry and the direction of the bias

TESS observes near the anti-Sun direction and WISE at 90° elongation. So for most stars
every sector falls about three months from the nearest WISE visit on both sides. The disk
excess at a sector is always interpolated, never measured. That is errors-in-variables and
biases β low. One clause in Method: "TESS and WISE view a star about three months apart,
so the excess at each sector is interpolated; this biases $\beta$ towards zero."

### S6. Small errors of fact or wording

- "not $A_1$ ($0.17$, $p\geq0.11$)": p_perm is 0.1055. Write "$p=0.11$ and 0.15".
- "stars brighter than $W1=8$": Jia et al. already require mean W1 of 7–12.5. Fine, but
  say the catalogue excludes W1 < 7 if you keep the saturation remark.
- "Of their 3,237 stars not flagged as eclipsing binaries": the catalogue has 3,300 and
  63 EBs. Fine. Say "3,300 candidates" once.
- "All analyses were checked on synthetic data first": a moderator cannot use this. Either
  cut it (saves 9 words) or say what was checked: "Injection tests recover $\beta$ and give
  calibrated nulls."
- Method does not say which flux column each product uses. It matters (M3).
- Two `%% TODO` DOI lines are still in the acknowledgments.
- Confirm the repo is public and tag the commit; the note promises "per-sector tables" at
  the URL, and `results/full/sector_metrics.csv` must actually be committed.

### S7. Disclose the relaxed go/no-go threshold, or drop the appearance of pre-planning

`RESULTS_FULL.md` records that the preset noise-floor criterion (0.3) failed (0.54) and was
relaxed. The note reports 0.54 but nothing implies a failed gate. That is acceptable for an
RNAAS. But "the planned control" language should not appear unless the failed gate does
too. The current draft says "a control", which is fine. Keep it that way.

---

## OPTIONAL

- **Title.** "Within a Be Star…" promises a per-star rule that the abstract withdraws.
  Suggest: "Be-star Candidates Are More Variable in TESS When Their WISE Disk Excess Is
  Larger".
- **Figure (packet question 8).** Keep the scatter: it shows the honesty about scatter.
  But panels (a) and (b) say the same thing (B1 and flicker are highly correlated). Replace
  (b) with β against absolute frequency for the main and SPOC samples. One two-panel figure
  then carries both results. Only do this after S1.
- **Phase sentence (packet question 9).** Keep, shortened: "Sectors on a brightening W1
  curve are more variable than the slope predicts, but the phase label depends on the
  interpolation, so I claim no phase effect."
- **Post hoc frequency analysis (packet question 6).** Fair to include. Move the caveat
  next to the result ("In a post hoc analysis, repeating the fit in 30 frequency bins…")
  so it is not three paragraphs away.
- **Shell stars (packet question 10).** The sign is right: high-inclination stars fade in
  the optical while the IR still brightens. TIC 55295028 at +0.91 mag per mag is a strong
  case. Give its common name if it has one; a specialist will look it up.
- **Physics otherwise.** W1 excess tracking inner-disk density within a star is fine. The
  phase thresholds are arbitrary but no claim rests on them. Štefl frequencies sit about
  10% below the main group, so the 0.9 f1 remark must wait for the recentred bins (S1).
- **Cite two more.** Huat et al. 2009 and Semaan et al. 2018 are the classic single-star
  and small-sample precedents (below). References do not count against the word limit.

## Word budget

The draft has 1,113 counted words, so 237 spare. M1 +20, M3 about +25, M4 about +45,
M5 about +20, S2 +25, S3 +12, S5 +22. Roughly +170, which fits. The abstract replacement
is about 146 words; recount with `aas_wordcount.py`.

---

## Answers to the packet's 10 questions

1. **Headline strength.** The claim is right for B1 and flicker. The title overstates
   ("Within a Be Star") and "about half go the other way" understates. See M5.
2. **QLP vs SPOC.** Not handled correctly yet. The gap is 1.6σ for A1 and is the size the
   pipelines' different normalisation predicts (M3). Lead with B1 (M2).
3. **Dilution.** The sign argument is right for SPOC. It does not apply to QLP.
   "Conservative" is the wrong word for the pooled sample (M3).
4. **White noise.** Not convincing as written; the SPOC ratio failure is omitted. The
   tercile split is the convincing argument (M4). I see no sign of scattered light or
   epoch effects: sector fixed effects change nothing (S3).
5. **Swap null.** Clear enough. Add the fixed-effects sentence (S3).
6. **Frequency dependence.** Fair with the caveat moved up. Fix the bin edge first (S1).
7. **Novelty.** I found no within-star, multi-year, population TESS × IR test. See below.
8. **Figure.** Keep the scatter, swap panel (b) for β vs frequency (Optional).
9. **Phase terms.** Keep one shortened sentence.
10. **Physics.** The candidates issue (M1) and "A1 is the coherent pulsation" (M2) are the
    two physics problems. Shell-star sign and Štefl placement are fine.

---

## Number spot-checks

All against `results/` files; my script also reproduces the main β values from
`sector_metrics.csv`.

| Text | File value | Match |
|---|---|---|
| 3,237 non-EB; 2,541 pass; 2,020 / 219 / 198 / 104 | `sample.csv`: 3,237 rows, 2,541 pass, same split | yes |
| 1,084 stars, 4,108 sectors; 755 QLP; 329 SPOC | `test_results.csv` main, QLP, SPOC rows | yes |
| A1 0.40 [0.31, 0.50] | 0.4043 [0.3114, 0.5037] | yes |
| B1 0.80 [0.70, 0.91]; B2 0.66; B3 0.39; red 0.74; flicker 0.84 [0.73, 0.96] | 0.8028 [0.7017, 0.9074]; 0.6578; 0.3905; 0.7402; 0.8440 [0.7293, 0.9631] | yes |
| white noise 0.33 | 0.3329 | yes |
| swap-null means −0.03 to +0.02 | −0.028 to +0.015 over 126 rows | yes |
| no cut: 2,449 stars, A1 0.29, B1 0.59 | 2,449; 0.2902; 0.5873 | yes |
| linear interpolation raises β by 4–11% | ratios 1.042–1.105 | yes |
| < 1 c/d: 0.71–0.92 | bins up to 0.98 c/d: 0.714–0.915 | yes |
| > 5 c/d: 0.31–0.37 | bins from 5.24 c/d: 0.306–0.370 (the 4.50–5.24 bin is 0.40) | yes |
| 0.73–1.16 below 0.7 f1 | bins to 0.63: 0.727–1.155 | yes |
| 0.35 "at the dominant group itself" | bin 1.00–1.26: 0.354; bin 0.79–1.00: 0.528 | number yes, description no (S1) |
| SPOC formal error 0.09, p 0.56 | 0.0941, p_perm 0.561 (p_swap 0.578) | yes |
| ratio to white: 0.41–0.51, p < 10⁻³; A1 0.07, p ≥ 0.33 | 0.407, 0.470, 0.511; 0.0714, p 0.327 / 0.531 | yes |
| k = 0.21 (0.11–0.30), 125 stars | 0.2131 (0.105–0.295), 125 | yes |
| "about 0.6 for A1 and 1.0 for B1" | 0.404 + 0.196 = 0.60; 0.803 + 0.196 = 1.00 | arithmetic yes; claim no (M3) |
| shells removed: change < 0.01 | 0.400, 0.796, 0.840 | yes |
| 0.54 rms; ρ 0.09–0.20 | 0.5439; 0.092–0.196 | yes |
| 48% negative A1 slopes | 0.484 (my recount: 525 of 1,084) | yes |
| SPOC: B1 0.69, flicker 0.74, p < 10⁻³; A1 0.17, p ≥ 0.11 | 0.685, 0.740, p 10⁻⁴; 0.173, p 0.1055 / 0.149 | yes, except "≥ 0.11" (S6) |
| RX J0440.9+4431 negative slope | −0.861, 3 sectors | yes |
| Figure legend 766 / 2,499 / 843 | 766 building, 843 dissipating, sum 4,108 | yes |
| "stochastic power … rises by about 8%" | B1 is rms amplitude; e^0.080 = 1.083 | number yes, word no (M2) |

Could not check: the 149 / 1,113 word counts (did not rerun the counter), the UAT keyword
IDs, the dataset DOIs, and the bibliography entries (I did not re-verify volumes and pages).

---

## Literature

**Novelty.** I ran 11 web searches (not ADS) and read the Jia et al. paper in part. I found
no paper that tests TESS variability against a multi-year IR or Hα disk record within
stars across a sample. Jia et al. use TESS only to find eclipsing binaries; their
conclusions reportedly list pulsation and outburst statistics as future work for the
community, with no specific plan. The catalogue is three weeks old and accepted in MNRAS,
so the obvious follow-up is open to anyone. That argues for submitting soon.

**Overlaps the draft should know about** (none pre-empts it):

- Jia et al. 2026, the source: https://arxiv.org/abs/2611.02450. "BeT candidates"; see M1.
- Jian et al. 2024, A&A 682, A59, the spectroscopic predecessor (916 stars):
  https://arxiv.org/abs/2311.08395. Use it for the confirmed-subset check.
- Huat et al. 2009, A&A 506, 95 (CoRoT, HD 49330): the pulsation spectrum changes through
  an outburst. https://www.aanda.org/articles/aa/abs/2009/40/aa11928-09/aa11928-09.html
- Semaan et al. 2018, A&A 613, A70 (15 CoRoT Be stars): in five, power-spectrum changes
  follow the outburst. Closest earlier "sample" result. Found via search; I did not open
  the paper. https://doi.org/10.1051/0004-6361/201629243 (check the DOI before citing.)
- Labadie-Bartz et al. 2018, AJ 155, 53 (KELT; outbursts in 28% of the sample, more in
  early types): https://arxiv.org/abs/1708.02594
- Labadie-Bartz et al. 2021, MNRAS (γ Cas frequency groups; brightenings usually come with
  larger group amplitudes; inner-disk inhomogeneities may add red noise):
  https://arxiv.org/abs/2012.06454. Relevant to M3 and S2.
- Nazé, Rauw & Pigulski 2020, MNRAS 498, 3171 (15 γ Cas stars; amplitudes change over time
  even without outbursts): https://arxiv.org/abs/2008.08334
- Baade et al. 2016, A&A 588, A56 (BRITE; Štefl frequencies as circumstellar):
  https://arxiv.org/abs/1602.01744
- Richardson et al., "The birth of Be star disks II" (λ Pav; 698 spectra with TESS):
  https://arxiv.org/abs/2512.08214. Newer than the 2025 paper the draft cites.
- HD 212044, A&A 2024 (brightness anticorrelates with Hα EW; combined frequency amplitude
  rises in Hα outbursts): https://www.aanda.org/articles/aa/full_html/2024/11/aa52215-24/aa52215-24.html
- "The optical photometric variability of Herbig Ae/Be stars from TESS", ApJ 2025:
  https://arxiv.org/abs/2510.08633. No clear link between light-curve morphology and
  near-IR colours. Relevant because Herbig stars are a likely contaminant (M1).
- QLP flux convention, for M3: Huang et al. 2020, RNAAS 4, 204 and 206;
  https://arxiv.org/abs/2011.06459

Not done: an ADS full-text search, and a check of the Labadie-Bartz northern-ecliptic
follow-up (no such paper surfaced).

---

## First-impression notes (written after reading only `ms.tex` and the figure)

Kept as written, before I opened the packet, results or code.

**Believable?** The low-frequency result (B1, red noise, flicker) is believable and
physically expected: Be stars are known to be photometrically noisier when they are feeding
a disk. The A1 result is not established by the text itself.

1. The abstract leads with A1 (+4% per 0.1 mag, p<1e-4) but the body concedes three things
   that undo it: the white-noise control rises by nearly the same β (0.33 vs 0.40), A1
   relative to white noise is flat (0.07, p≥0.33), and the cleaner SPOC subset does not
   confirm it (0.17, p≥0.11). A moderator reading abstract then Results will see the
   abstract as overstated. A1 looks like "the whole spectrum floor goes up", not a
   pulsation result.
2. The white-noise control failing is the biggest technical worry. The rebuttal (formal
   SPOC flux error flat) is on 329 SPOC stars while the β=0.33 is on a 70% QLP sample. Is
   white-noise β also 0.33 in the SPOC-only subset? If it is ~0 there, the rise is a QLP
   artefact; if it persists, it is probably spectral leakage / red-noise tail from the
   low-frequency signal into 15–24 c/d (which abuts the 24 c/d Nyquist of 30-min bins), and
   "relative to white noise" over-corrects. Either way the text needs to say which.
3. 5-d spline detrend + B1 band at 1.25–4 d + "flicker" = 1-d-smoothed residuals: these
   metrics will catch the photometric outburst itself (the brightening/fading event of disk
   feeding), not only "stochastic variability". So the finding may reduce to "sectors with
   larger W1 excess more often contain an outburst" — still a result, but "the star is more
   variable" and "stochastic" over-interpret; the variable light may be circumstellar.
4. "48% of stars have a negative A1 slope": 52% positive of 1,084 is not distinguishable
   from a coin flip (sign test p≈0.2). So the pooled slope is carried by a minority of
   high-leverage stars with large W1 range. Need: slope/sign fraction as a function of each
   star's W1 excess range; a leverage-robust estimate. Drop-one-star is not enough; drop
   top-N-leverage or a per-star-median is.
5. Dilution "correction" assumes the variability is photospheric and the extra light is
   disk continuum. If the low-frequency variability is itself circumstellar, no correction
   applies. "Conservative" is too strong. Also k from SAP across sectors: SAP
   sector-to-sector zero points are not reliable (different camera/CCD/aperture) — how was
   k measured?
6. Sample: Jia et al. select Be stars *whose disks changed state* — a variability-selected
   sample. Text should say the result is for disk-variable Be stars. Are they confirmed Be
   stars or candidates (LAMOST? contaminants: Herbig, B[e], YSOs)? arXiv:2611.02450 is
   weeks old — risk the same group is doing the TESS cross-match.
7. Title "Within a Be Star…" promises a per-star rule the abstract's last sentence
   withdraws.
8. Errors in the regressor (W1 spline interpolated across ≤250 d gaps) attenuate β; fine,
   conservative, but worth one clause. Also the time-scale mismatch: 27-d sector vs 6-month
   WISE cadence.
9. Epoch confound (QLP processing, FFI cadence 30→10→3.3 min, changing gap structure) — the
   swap null is a good answer. Is the swap done within pipeline / within similar sector
   sets? Otherwise stars with non-overlapping epochs get extrapolated splines.
10. Discussion: Niwano "four of five is consistent" — fine (P≈0.16) but BeXRBs are not
    "single stars". "more ejection episodes per sector" is speculation, labelled as
    plausible — OK.
11. Missing obvious prior art from memory: Huat et al. 2009 (CoRoT HD 49330: pulsation
    spectrum changes through an outburst), Baade et al. 2016–2018 (BRITE), Semaan et al.
    2018 (CoRoT Be), Nazé et al. 2020 (TESS Be red noise), Bernhard et al. 2018 /
    Labadie-Bartz 2017–18 (KELT + outbursts), Rivinius et al. 2013 review, and the earlier
    WISE Be variability catalogue (Jian et al. 2024?). Novelty claim "extends to 1,084
    stars and years" needs these checked.
12. Housekeeping: two TODOs for DOIs still in the file; "All analyses were checked on
    synthetic data first" is unverifiable as written; repo URL must be public and tagged;
    personal gmail as contact; figure is .png in tree but \plotone wants .pdf.

**New?** Probably yes as a population-scale, within-star, TESS×NEOWISE test. Not new as
physics.

**What a moderator pushes on:** abstract vs body mismatch on A1; the white-noise control;
variability-selected sample; one figure OK; length.

**How the first impressions held up.** Points 1, 3, 5, 6, 7 became M2, S2, M3, M1 and the
title note. Point 2 resolved in the draft's favour (M4): the control is contaminated by
signal, and the note can say so with evidence. Point 4 was half right: the sign test is
indeed null for A1, but the result is robust to removing the highest-leverage stars (M5).
Point 9 is fine (S3). In point 12 the figure PDF does exist in `rnaas/`.
