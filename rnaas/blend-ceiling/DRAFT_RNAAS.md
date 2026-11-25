# DRAFT: Catalogued Eclipsing Binaries Flatter Light-curve Classifiers of Long-period TESS Candidates

> **Superseded by `rnaas.tex` (2026-11-25)**, which rewords the introduction after the reference check; see
> `numbers.md`.

Jon McBee (Coracle)

> Draft status (2026-11-24): every number traces to `results/auc.txt` (amendment B, the model
> reported here), `results/amendmentA/auc.txt` (the first, pre-registered model) and
> `results/blend_fraction.txt`; the scripts are in this folder and the rules in `RULES.md`. The
> quiet training class (amendment B) was added **after** the first model failed the
> pre-registered no-transit check, and the text says so. RNAAS allows one figure or one table, so
> this draft has the figure only. Not submitted, and not checked against the RNAAS word-count tool.

## Abstract

Classifiers of TESS candidates are often tested against catalogued eclipsing binaries (EBs).
At short periods that flatters them, because many real false positives are blends. We repeat
the test at long periods (P > 10 d or single transits): 336 confirmed or known planets against
215 TFOP false positives or false alarms. A gradient-boosting classifier, never trained on TFOP
false positives, reaches AUC 0.959 (0.945–0.971) against catalogued EBs and 0.959 against
windows with no transit, but 0.649 (0.603–0.695) against TFOP false positives. Against those
whose TOI comments name a blend it scores 0.525 (0.441–0.607), consistent with chance. Comments
flag 22% of the false positives as blends, which caps any light-curve-only classifier at AUC
0.89. Within depth bins the AUC is 0.57–0.63. A first model without no-transit training windows
failed a no-transit check but gave the same picture.

## 1. Introduction

Machine-learning classifiers for TESS candidates are usually scored against catalogued EBs.
Pagés Navarrete (2026) showed that this flatters them. A short-period CNN reached AUC 0.970
against catalogued EBs but 0.500 against real TFOP false positives. In that sample 42% of the
false positives were blends, which sets a ceiling AUC_max = 1 − f/2 on any classifier that sees
only the target's light curve. That test covered P < 1.5 d.

Salinas et al. (2026) recently introduced a CNN that characterises single transits in
long-period candidates. They validated it on 120 confirmed TOIs with P > 10 d, but against no
TFOP false positives. Their network regresses depth, duration and mid-time, and gives no
detection score, so we cannot score candidates with it.

Here we run the two-negative-set test at long periods with our own light-curve-only baseline,
and measure the blend fraction among long-period TFOP false positives.

## 2. Data and method

**Sample.**
- From the ExoFOP TOI table (downloaded 2026-11-24), the TOIs with TFOPWG disposition CP/KP
  (positives) or FP/FA (negatives) and P > 10 d, or no period.
- Light curves are SPOC 2-min, else QLP, one per sector.
- Windows are 99.5 h long, binned to 30 min and centred on each predicted transit with timing
  error 3σ < 25 h, and need at least half their bins filled.
- 336 positives and 215 negatives have at least one valid window (median 8 valid windows per TOI).

**Training.**
- Positives: CP/KP windows.
- Easy negatives:
  - 410 Prša et al. (2022) EBs with P > 10 d that are not TOIs, at their primary eclipses;
  - 1146 "TCE-like" windows: the highest-SNR box dip in each light curve of those stars, at
    least 3 d from any known transit or eclipse;
  - 1182 "quiet" windows: random windows from the same stars, at least 3 d from any known
    transit or eclipse.
- TFOP FP/FA are never used in training.

**Classifier.**
- Gradient boosting on 15 features from the window flux alone: box-search SNR and depth,
  duration, noise and red-noise ratio, V/U shape, asymmetry, second-dip SNR, skewness,
  brightening and gap fraction.
- Cross-fitted over 5 folds by TIC, so no star is scored by a model that saw it.
- An object's score is its maximum over its windows.

**Blends.** An FP/FA counts as a blend under either of two rules:
- **(C)** its TOI comments name an NEB, a centroid offset or an off-target source;
- **(G)** a Gaia DR3 neighbour within 42″ could produce the TOI depth on its own (G-band flux).

We give f and the ceiling with Wilson intervals, and AUCs with 2,000-sample bootstrap 95% CIs.

**Order of work.** The rules were committed before any window was scored. The first model
(model A) had no quiet windows. It failed the pre-registered no-transit check: 50 random
out-of-transit windows from planet hosts scored a median of 0.91, against 0.99 for planets,
because every training negative was a strong dip. We then added the quiet class as a dated rule
change, with a pass criterion fixed in advance (median no-transit score below 0.5 and below the
FP/FA median), and retrained once (model B). We report B and give A's values in brackets.

## 3. Results

**Blend fraction.**
- Rule C: 47 of the 216 false positives are blends (215 have a valid window and are scored),
  f = 0.218 (0.168–0.277), so AUC_max = 0.891 (0.861–0.916). This is half the short-period
  fraction.
- Rule G alone flags 70% of the false positives but also 60% of the confirmed planets, so it
  barely tells the two apart. C or G gives f = 0.718 and AUC_max = 0.641, but we regard that as
  an overcount.

**AUCs**, CP/KP against each negative set (Figure 1; model A in brackets):
- catalogued EBs: 0.959 (0.945–0.971) [0.949];
- no-transit windows: 0.959 (0.939–0.976) [0.786];
- all TFOP FP/FA: 0.649 (0.603–0.695) [0.643];
- non-blend FP/FA (rule C): 0.684 (0.637–0.730) [0.668];
- blend FP/FA (rule C): 0.525 (0.441–0.607) [0.555].

**Checks.**
- No-transit check, passed by B: the 50 out-of-transit windows score a median of 0.46
  (10–90%: 0.21–0.74), against 0.98 for planets and 0.97 for TFOP false positives. A
  feature model without the three depth-scale features and a small 1D CNN, both trained the
  same way, fail it (medians 0.56 and 0.68); they reach 0.664 and 0.629 against all FP/FA.
- The false positives are deeper than the planets (median 4,090 vs 1,680 ppm, KS p = 2e-9)
  and slightly fainter (Tmag 10.4 vs 10.1).
- Within four fixed depth bins the AUC is 0.57–0.63. The lower limits of the 95% CIs run from
  0.49 to 0.52, so much of the 0.649 comes from depth.

## 4. Discussion

At long periods, as at short ones, catalogued EBs flatter a light-curve classifier: AUC 0.96
against them, 0.65 against real TFOP false positives. Blends are scored at chance. Teaching the
model what an empty window looks like fixed the no-transit check (0.79 to 0.96) but left the
TFOP comparison where it was. The two failures are separate: easy negatives say little about
the false positives TFOP actually rejects.

The blend ceiling from TOI comments, 0.89, sits well above what this baseline reaches. So for
it, blends are not yet the binding limit. The non-blend false positives, scored at 0.68, are.

Any long-period classifier, including a detection head added to networks like that of
Salinas et al., should be scored against TFOP FP/FA, within depth bins, and on no-transit
windows. Our sample, rules and code are public.

**Figure 1.** Left: ROC curves of CP/KP against all TFOP FP/FA, non-blend FP/FA and blend
FP/FA (comment rule), for model B. Right: AUC with 95% bootstrap CIs per negative set, for model
B (filled) and model A (open). The horizontal lines are the blend ceilings 1 − f/2 (bands:
Wilson intervals) for the comment rule and for comment or Gaia, and the dashed line is chance.

## Acknowledgments

This work used data from TESS (SPOC, QLP) via MAST, ExoFOP-TESS, and Gaia DR3 (served by the
ARI Heidelberg TAP mirror). Analysis was assisted by Claude, an AI model by Anthropic. The
author designed the study, checked the analysis and takes responsibility for the content.
(Check the current AAS policy on AI use before submission.)

## References

- Gaia Collaboration, Vallenari, A., et al. 2023, A&A, 674, A1
- Guerrero, N. M., et al. 2021, ApJS, 254, 39
- Huang, C. X., et al. 2020, RNAAS, 4, 204
- Jenkins, J. M., et al. 2016, Proc. SPIE, 9913, 99133E
- Lightkurve Collaboration 2018, ascl:1812.013
- Pagés Navarrete, (2026), arXiv:2610.42009 (to verify: author list and title)
- Pedregosa, F., et al. 2011, JMLR, 12, 2825
- Prša, A., et al. 2022, ApJS, 258, 16
- Salinas, H., et al. (2026), arXiv:2610.24882 (to verify: author list and title)
