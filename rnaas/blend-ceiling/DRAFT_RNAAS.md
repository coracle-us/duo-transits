# DRAFT: Catalogued Eclipsing Binaries Flatter Light-curve Classifiers of Long-period TESS Candidates

Jon McBee (Coracle)

> Draft status (2026-11-24): every number traces to `results/auc.txt` and
> `results/blend_fraction.txt` (scripts in this folder, rules in `RULES.md`). **The
> classifier fails the pre-registered sanity check** (section 3), and the text says so. Not
> submitted, and not checked against the RNAAS word count tool.

## Abstract

Classifiers of TESS candidates are often tested against catalogued eclipsing binaries (EBs).
At short periods that flatters them, because many real false positives are blends. We repeat
the test at long periods (P > 10 d or single transits). The sample is 336 confirmed or known
planets and 215 TFOP false positives or false alarms from the TOI table. We use a
gradient-boosting classifier trained without any TFOP false positives. It reaches AUC 0.949
(0.934–0.962) against catalogued EBs but 0.643 (0.595–0.689) against TFOP false positives,
and 0.555 (0.472–0.630), consistent with chance, against those whose TOI comments name a blend.
Comments flag 22% of the false positives as blends, which caps any light-curve-only classifier
at AUC 0.89. Within depth bins the AUC falls to 0.54–0.61. Random out-of-transit windows score
like planets, so this baseline responds to dip strength, not transit presence.

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
    least 3 d from any known transit or eclipse.
- TFOP FP/FA are never used in training.

**Classifier.**
- Gradient boosting on 15 features from the window flux alone: box-search SNR and depth,
  duration, noise and red-noise ratio, V/U shape, asymmetry, second-dip SNR, skewness,
  brightening and gap fraction.
- Cross-fitted over 5 folds by TIC, so no star is scored by a model that saw it.
- An object's score is its maximum over its windows.
- Secondary models: one without the three depth-scale features, and a small 1D CNN.

**Blends.** An FP/FA counts as a blend under either of two rules:
- **(C)** its TOI comments name an NEB, a centroid offset or an off-target source;
- **(G)** a Gaia DR3 neighbour within 42″ could produce the TOI depth on its own (G-band flux).

We give f and the ceiling with Wilson intervals, and AUCs with 2,000-sample bootstrap 95% CIs.
All rules were committed before any window was scored.

## 3. Results

**Blend fraction.**
- Rule C: 47 of the 216 false positives (one has no valid window, so 215 are scored) are blends, f = 0.218 (0.168–0.277), so
  AUC_max = 0.891 (0.861–0.916). This is half the short-period fraction.
- Rule G alone flags 70% of the false positives but also 60% of the confirmed planets, so it
  barely tells the two apart. C or G gives f = 0.718 and AUC_max = 0.641, but we regard that as
  an overcount.

**AUCs (CP/KP vs each set, Figure 1):**

| Negative set | AUC (95% CI) |
|---|---|
| Prša EBs | 0.949 (0.934–0.962) |
| All TFOP FP/FA | 0.643 (0.595–0.689) |
| Non-blend FP/FA (C) | 0.668 (0.621–0.716) |
| Blend FP/FA (C) | 0.555 (0.472–0.630) |

The shape-only model gives 0.643 and the CNN 0.605 against all FP/FA.

**Checks.**
- The false positives are deeper than the planets (median 4,090 vs 1,680 ppm, KS p = 2e-9)
  and slightly fainter (Tmag 10.4 vs 10.1).
- Within four fixed depth bins the AUC is 0.54–0.61, and every 95% CI includes 0.5. Most of the 0.643 therefore comes from depth.
- **Sanity check, failed:** 50 random out-of-transit windows from planet hosts score a median
  of 0.91, close to the planets (0.99). The easy negatives are all strong dips (median box SNR
  125 for EBs and 17 for TCE-like windows, against 9 for planets). So the model learned "a
  modest dip or none means planet", not "transit present". The rules forbade retraining after
  seeing this, and we report the model as it stands.

## 4. Discussion

At long periods, as at short ones, catalogued EBs flatter a light-curve classifier: AUC 0.95
against them, 0.64 against real TFOP false positives. Blends are scored at chance.

The blend ceiling from TOI comments, 0.89, sits well above what this baseline reaches. So for
it, blends are not yet the binding limit; its easy training negatives are. The sanity failure
makes the same point from the other side: easy negatives teach a model to reject strong
signals, not to recognise transits.

Any long-period classifier, including a detection head added to networks like that of
Salinas et al., should be scored against TFOP FP/FA, within depth bins, and on no-transit
windows. Our sample, rules and code are public.

**Figure 1.** Left: ROC curves of CP/KP against all TFOP FP/FA, non-blend FP/FA and blend
FP/FA (comment rule). Right: AUC with 95% bootstrap CIs per negative set. The horizontal lines
are the blend ceilings 1 − f/2 (bands: Wilson intervals) for the comment rule and for comment
or Gaia, and the dashed line is chance.

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
