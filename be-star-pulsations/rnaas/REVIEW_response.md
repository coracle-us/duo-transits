# Response to the Mythos review (draft v1 → v2, 2026-11-24)

**Files**
- Review: `rnaas/REVIEW_mythos.md`.
- Revised draft: `rnaas/ms.tex` (v2).
- Provenance: `rnaas/numbers.md`.
- Reruns: `src/c5_review_checks.py`.

**Outputs of the reruns** (all in `results/checks/`)
- `review_checks.csv`: every row went through `stats4.run_test`, with both nulls at 10,000 draws.
- `review_extras.csv`: bootstrap only.
- `review_signs.csv`.

## Must fix

### M1. Candidates, not Be stars

**Done.** The title, abstract, Method and caption now say "Be-star candidates". Method also states
the selection and that most stars lack spectroscopic confirmation.

**Confirmed subset.** CDS was still down on 2026-11-24, so I could not cross-match Jian et al.
2024. I used Gaia DR3 ESP-ELS instead, which classifies stars by Hα emission in their RP spectra.

| Gaia ESP-ELS label | Main-sample stars |
|---|---|
| beStar | 265 |
| HerbigStar | 0 |
| WR | 1 |

**Result on the confirmed stars** (256 stars with at least 2 usable sectors):

| Metric | β | p_swap |
|---|---|---|
| B1 | 0.58 [0.29, 0.87] | 5e-4 |
| red noise | 0.57 | 2e-4 |
| flicker | 0.59 | 6e-4 |
| A1 | 0.17 | 0.17 (n.s.) |

The abstract and Results now report this.

### M2. Lead with B1; A1 is not the pulsation

**Done.** The abstract follows the review's draft, adapted. "Power" is now "amplitude".

**A1 rerun:**

| Stars | n | A1 β |
|---|---|---|
| f1 stable within 10% | 510 | 0.23 [0.10, 0.36], p_swap 0.004 |
| stable, with f1 > 0.8 c/d | 445 | 0.08, p_swap 0.29 |
| f1 moving between sectors | 574 | 0.47 |

### M3. Dilution applies to SPOC only

**Done.** I checked the claim first:
- Huang et al. 2020, §II item 7: QLP adds the TIC-expected flux to the difference flux, and §III
  anchors each orbit's median to the TIC magnitude.
- The cached QLP `SAP_FLUX` has a median of about 1 in every sector.

**SPOC stars with k = 0.213**, all through the pipeline:

| Metric | SPOC, corrected | p_swap | QLP |
|---|---|---|---|
| A1 | 0.37 | 0.004 | 0.52 |
| B1 | 0.88 | 1e-4 | 0.86 |
| flicker | 0.94 | 1e-4 | 0.90 |

The bootstrap QLP − SPOC gap at k = 0 is +0.35 / +0.18 / +0.16 (A1 / B1 / flicker). The gap
predicted from dilution is +0.20.

### M4. White-noise control

**Done.** Terciles by each star's median white-noise level, through the pipeline:

| Tercile | white-noise β | B1 β | A1 β |
|---|---|---|---|
| noisiest | −0.00 (p_swap 0.87) | 0.63 | 0.43 |
| middle | 0.49 | — | — |
| quietest | 0.61 | — | — |

- With ln B1 in the fit, the white-noise slope on excess is −0.03 [−0.12, 0.06].
- The SPOC ratio failure (p_swap ≥ 0.13) is now stated.

### M5. "About half go the other way"

**Done.** The text now gives the fraction of stars with a positive slope, by how much their disk
excess changed:

| Excess range | A1 positive | B1 positive | flicker positive |
|---|---|---|---|
| < 0.05 mag | 47% | 50% | 50% |
| > 0.4 mag | 65% | 71% | 72% |

**Dropping the top leverage stars:**

| Stars dropped | B1 β |
|---|---|
| top 10 | 0.88 |
| top 50 | 0.94 |
| top 100 | 1.03 |

All are significant under both nulls.

## Should fix

**S1. Recentred bins.** Done. The bin centred on f1 gives 0.43 [0.34, 0.53]; above f1, the
bins give 0.28–0.49. In stable-f1 stars the centred bin is 0.30, the second-lowest of their bins.
"The same as at high frequency" is gone. The text now says "no separate response of the dominant
pulsation".

**S2. Outbursts in B1 and flicker.** Done, one sentence in the Discussion.

**S3. Sector fixed effects.** Done. Every β changes by less than 0.03 (bootstrap only). Within
star × cadence groups, through the pipeline:

| Metric | β | p_swap |
|---|---|---|
| B1 | 0.58 | 0.003 |
| flicker | 0.48 | 0.016 |
| A1 | 0.22 | n.s. |

The star × cadence rows are not in the note, to save words.

**S4. Niwano paragraph.** Done: "sample stars", and "their statistic is not my A1". No new
cross-match count is added.

**S5. Observing geometry and errors-in-variables.** Done, one clause.

**S6. Small fixes.** Done:
- the p-values;
- the injection-test sentence;
- the flux columns: M3 explains how SPOC and QLP normalise;
- the DOIs, now resolved: GSFC-ELEANOR-LITE 10.17909/j2yt-t417 and NEOWISE-R IRSA144.

Not done:
- The W1 < 8 remark is unchanged. The catalogue's 7–12.5 range was not re-checked.
- The flux column names are not written out in Method, to save words; they are in `numbers.md`.

The repo is public, and `results/full/sector_metrics.csv` is committed. The release still needs a
tag.

**S7. "Planned" control.** Done: the note says "a control".

## Optional items

| Item | Status |
|---|---|
| Title | Taken: "Be-star Candidates Are More Variable in TESS When Their WISE Disk Excess Is Larger" |
| Figure panel (b) → β vs frequency | **Not done**; left for Jon to decide |
| Phase sentence | Shortened, as suggested |
| Post hoc caveat next to the result | Done |
| Huat 2009 and Semaan 2018 citations | Added; both checked on Crossref, and the Semaan abstract supports the citation |
| Labadie-Bartz 2021, Nazé 2020, Richardson 2025 | Not cited, to keep within the word limit |
| Common name for TIC 55295028 | Not looked up (SIMBAD is down) |

## Word count

- Abstract: 146 words.
- AAS-counted text: 1,283 words, about 1,310 by the AAS count once the 2% calibration is applied.
- Limit: 1,350.

## Second review (`REVIEW_mythos_v2.md`, 2026-11-24)

Verdict: ready to submit after four wording fixes, with no analysis outstanding. All of them are
applied in `ms.tex`.

**Must fix**
- "Half the stars" is now "Nearly half (48%) of the stars".
- "−0.00" is now "0.00".
- "Elsewhere it tracks B1" is now "Across the sample it tracks B1".
- B2 and B3 are no longer defined in Method, since the note reports neither.

**Should fix**
- "therefore carries" is now "therefore appears to carry".
- The leverage drop is now described as "the 100 highest-leverage stars".
- The length margin is wider. I dropped "The swap nulls centre on zero" and shortened the phase
  sentence. Counted text is now 1,264 words, about 1,295 by the AAS count.

**Mythos's own check of the Gaia subset** (bootstrap only, not added to the note)
- All 256 stars have a Be-class probability above 0.5.
- The 227 stars with probability above 0.8 give B1 0.79 ± 0.32.

