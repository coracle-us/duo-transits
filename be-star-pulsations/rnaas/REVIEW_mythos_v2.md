# Second review of the Be-star RNAAS draft (`rnaas/ms.tex` v2, 2026-11-24)

Reviewer: Claude (Mythos). Requested by the "BE star pulsations analysis" session on Jon's
behalf. Nothing was edited, committed, sent or posted; this file is the only new one.

## Verdict: ready to submit after the small text fixes in "Must fix"

No analysis is outstanding. All five must-fix items from the first review are resolved, and
the reruns went through the pipeline with both nulls. I checked 30 numbers in v2 against
`review_checks.csv`, `review_extras.csv`, `review_signs.csv` and the earlier result files:
all match. What remains is four wording errors (10 minutes of edits) and the pre-submission
housekeeping listed at the end. None of it changes a result.

---

## Answers to the five questions

### 1. Is Gaia ESP-ELS an acceptable substitute for the confirmed subset? Yes, as worded.

- **What it is.** ESP-ELS flags a star when the Hα pseudo-equivalent width in the
  low-resolution RP spectrum shows emission at 3σ, then a random forest assigns a class
  (Be, Herbig Ae/Be, T Tauri, active M dwarf, WC, WN, planetary nebula). The Be and Herbig
  classes also use the GSP-Phot stellar parameters.
- **What it supports.** "256 stars with Hα emission in their Gaia spectra" (abstract) and
  "Gaia DR3 classifies as Be stars from Hα emission" (Results) are both accurate. It is a
  subset with detected emission, which is what a sceptic of a photometric catalogue wants.
- **What it does not support.** It is not classical spectroscopic confirmation, and the Be
  versus Herbig split rests on photometric stellar parameters. The draft does not claim
  either, so no change is needed.
- **Quality of the labels** (my check, from the cached `gaia_dr3_params.parquet`): all 256
  have a Be-class probability above 0.5; the median is 0.96. None of the 1,084 main-sample
  stars is labelled Herbig, and one is labelled WN. Restricting to probability > 0.8
  (227 stars) gives B1 0.79 ± 0.32 (bootstrap only, no nulls). The result does not depend
  on the weaker labels.
- **Citation.** Creevey et al. 2023, A&A 674, A26 is "Gaia DR3: Apsis I, methods and
  content overview" (arXiv:2206.05864), and it describes ESP-ELS. Correct.
- **Not blocking:** the Jian et al. 2024 cross-match is still worth doing when CDS is back,
  for your own records.

### 2. Is the dilution paragraph now correct? Yes.

The QLP statement is right (difference flux plus TIC-expected flux; Huang et al. 2020), the
"if the variable light is photospheric" condition is stated, and the numbers match:
SPOC corrected 0.369 / 0.881 / 0.936 with p_swap ≤ 0.0042, QLP 0.52 / 0.86 / 0.90. "Close
to" is fair for A1 (0.37 against 0.52) given the ±0.2 uncertainty on the gap.

### 3. Abstract and white-noise paragraph: accurate, with two wording slips (M2, M3 below)

- Abstract (146 words): every claim traces to a file. "The dominant periodic signal shows
  no separate response" is supported: the bin centred on f1 gives 0.43 [0.34, 0.53] and
  the bins above it 0.28–0.49.
- White-noise paragraph: the numbers are right (noisiest third −0.002, p 0.98 / 0.87;
  B1 0.63; slope with ln B1 in the fit −0.03 [−0.12, 0.06]). Two phrases need fixing.

### 4. New errors introduced by the revision: four small ones, listed under Must fix

### 5. Ready? Yes, after Must fix and the housekeeping list

---

## MUST FIX (wording only)

### M1. "Half the stars" is 48%

> "Half the stars change by less than 0.1\,mag between sectors and cannot show a trend"

519 of 1,084 is 47.9%, and the median range is 0.108 mag, so fewer than half are below
0.1. Write "Nearly half the stars (48\%) change by less than 0.1\,mag…".

### M2. "−0.00"

> "it is flat ($-0.00$, $p\geq0.87$)"

Write "$0.00$". A signed zero reads as a typo.

### M3. "Elsewhere it tracks B1" misdescribes the correlation

> "Elsewhere it tracks $B_1$ (within-star correlation 0.47)"

The 0.47 in `review_extras.csv` is over all 1,084 stars, not the other two thirds. Write
"Across the sample it tracks $B_1$ (within-star correlation 0.47)".

### M4. B2 and B3 are defined but never reported

> "the rms in the bands 0.25--0.8, 0.8--2 and 2--4\,c\,d$^{-1}$ ($B_1$--$B_3$)"

v2 dropped the B2 and B3 results but kept their definitions. Either define only B1
("the rms at 0.25--0.8\,c\,d$^{-1}$ ($B_1$)", which saves about 8 words) or restore
"0.66 for $B_2$, 0.39 for $B_3$" in the first Results sentence. I would cut: the narrow
frequency bins already cover those bands.

---

## SHOULD FIX

### S1. Keep a wider margin on length

The counter gives 1,283 and you estimate about 1,310 by the AAS count, against 1,350.
AAS82207R1 was bounced for length, so 40 words is thin. M4 saves about 8. Two more cuts
that lose nothing:

- "The swap nulls centre on zero." (6 words). The next sentence's robustness list makes the
  same point.
- "I did not correct the $p$-values for the several metrics and samples; the main results
  sit at the $10^{-4}$ floor." Keep, but it can lose "the several".

### S2. "therefore" in the white-noise paragraph

> "The band therefore carries the tail of the stochastic signal, so dividing by it
> over-corrects"

The evidence is strong but indirect. "The band therefore appears to carry…" costs two
words and removes the one place a specialist could say the note asserts a mechanism.

### S3. Say what "largest disk changes" means

> "when dropping the 100 stars with the largest disk changes ($B_1$ 1.03)"

The script ranks by the within-star sum of squared excess (leverage), which also depends
on the number of sectors. "the 100 highest-leverage stars" is exact and the same length.

---

## OPTIONAL

- **Hα subset context.** 827 of the 1,084 stars have no ESP-ELS label. That is expected:
  the low-resolution spectra only catch strong emission (median pEW here is −0.5 nm). A
  reader may take 256 of 1,084 as "only a quarter are real". If words allow: "…the 256
  stars with H$\alpha$ emission strong enough for Gaia's low-resolution spectra".
- **"always interpolated".** True for all but a handful: the bracketing-visit gap is
  156–206 d for 90% of sectors, and only 4 stars have 10 or more sectors. Leave it.
- **Discussion, first sentence.** It reads as if Labadie-Bartz et al. 2022 (432 stars) were
  a "small sample". Harmless; fix only if you are editing nearby.
- **Figure panel (b).** Still redundant with (a). Your call, as the response says. The
  recentred bins now exist, so the β-against-frequency panel could be built from
  `review_checks.csv`.

---

## Spot-check of v2 numbers

| Text | File value | Match |
|---|---|---|
| B1 0.80 [0.70, 0.91]; red 0.74; flicker 0.84 [0.73, 0.96] | `test_results.csv`: 0.803 [0.702, 0.907]; 0.740; 0.844 [0.729, 0.963] | yes |
| no cut: 2,449 stars, B1 0.59 | 2,449; 0.587 | yes |
| drop top 100: B1 1.03 | `review_checks.csv`: 1.030, 984 stars, p 10⁻⁴ both nulls | yes |
| sector fixed effects change β by < 0.03 | `review_extras.csv`: largest change +0.026 (B1) | yes |
| 256 Hα stars: B1 0.58 [0.29, 0.87], p_swap 5×10⁻⁴ | 0.575 [0.288, 0.868], 0.0005; 256 stars, 996 sectors | yes |
| < 1 c/d 0.71–0.92; > 5 c/d 0.31–0.37 | `freq_beta.csv` (unchanged from v1) | yes |
| 0.74–1.17 below 0.9 f1 | 12 bins: 0.742–1.174 | yes |
| 0.43 [0.34, 0.53] centred on f1 | 0.432 [0.342, 0.527], bin 0.891–1.122 | yes |
| 0.28–0.49 above | 12 bins: 0.282–0.487 | yes |
| A1 0.40 overall | 0.404 | yes |
| 510 stable-f1 stars: A1 0.23 [0.10, 0.36] | 0.227 [0.098, 0.363]; p_perm 0.030, p_swap 0.0035 | yes |
| 0.08 (p ≥ 0.29) where f1 > 0.8 c/d | 0.078; p 0.37 / 0.29; 445 stars | yes |
| white noise 0.33 | 0.333 | yes |
| noisiest third: flat, p ≥ 0.87 | −0.002; p 0.985 / 0.868 | yes |
| …while B1 0.63, p < 10⁻⁴ | 0.626; 10⁻⁴ both | yes |
| within-star correlation 0.47 | 0.470 (all stars; see M3) | yes |
| slope with ln B1 in fit −0.03 ± 0.09 | −0.033 [−0.123, 0.058] | yes |
| SPOC formal error 0.09, p 0.56 | 0.094; 0.561 | yes |
| ratios 0.41–0.51, p < 10⁻³; SPOC p_swap ≥ 0.13 | 0.407–0.511; SPOC 0.133–0.218 | yes |
| k = 0.21 (0.11–0.30), 125 stars | 0.213 (0.105–0.295), 125 | yes |
| SPOC corrected 0.37 / 0.88 / 0.94, p_swap ≤ 0.005 | 0.369 / 0.881 / 0.936; 0.0042 / 10⁻⁴ / 2×10⁻⁴ | yes |
| QLP 0.52 / 0.86 / 0.90; 755 and 329 stars | 0.518 / 0.861 / 0.895 | yes |
| shell stars removed: change < 0.01 | 0.400 / 0.796 / 0.840 | yes |
| 0.54 rms; ρ 0.09–0.20 | 0.544; 0.092–0.196 | yes |
| "Half the stars" below 0.1 mag | 519 of 1,084 = 47.9% | no (M1) |
| 177 stars > 0.4 mag, 71% positive B1 | `review_signs.csv`: 177; 0.712 | yes |
| 48% negative A1 slopes | 0.484 | yes |
| abstract: about 8% and about 3% per 0.1 mag | e^0.080 = 1.083; e^0.033 = 1.034 | yes |
| 3,300 candidates; 3,237; 2,541; 2,020 / 219 / 198 / 104 | Jia et al.; `sample.csv` | yes |

I also reran the Hα subset independently from the per-sector table and the Gaia cache:
256 stars, B1 0.57 ± 0.28. It agrees.

Not checked: the word counts, the two new dataset DOIs, the Huat and Semaan bibliography
entries (your response says Crossref-checked), and the UAT keyword IDs. From memory the
four IDs (142, 235, 1625, 2109) are right, but I could not verify them.

---

## Housekeeping before submission (already on your list)

- Compile in Overleaf with aastex701. It has not been compiled, and the AAS word count
  comes from the compiled file.
- Recheck the current AAS policy on AI use. The disclosure section now mentions a review by
  a separate AI agent, which is accurate.
- Verify the UAT IDs.
- Tag the repo commit that matches the submitted text.
