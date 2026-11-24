# Rules for the blend-ceiling test (fixed before any data)

Written 2026-11-24, before the TOI table was downloaded, any light curve was fetched or any
window was scored. Every later change goes in the dated "Changes" section at the end, with its
reason, and the original text stays.

**Results are reported whatever they show**, including an AUC at or below chance, a ceiling
that the measured AUC already exceeds, or a sample too small to say anything.

## 1. Sample (Step 2)

- **Source.** The ExoFOP-TESS TOI table, full CSV
  (`https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv`), saved
  unedited as `data/toi_YYYY-MM-DD.csv`. The download date and the row count go in
  `data/README.md`.
- **Labels**, from the `TFOPWG Disposition` column, exact match after stripping spaces:
  - positives: `CP` or `KP`;
  - negatives: `FP` or `FA`;
  - every other value (PC, APC, blank) is out.
- **Long-period cut:** `Period (days)` > 10, or blank, or 0 (single transits and duos).
- **Light curves.** Per sector, SPOC 2-min (author SPOC, 120 s) if it exists, else QLP (any
  cadence). Found with `lightkurve.search_lightcurve` on MAST. TESS-SPOC, TGLC and 20 s
  data are not used.
- **Predicted transit times:** T = Epoch + k·P, with k ranging over every sector that has a
  light curve. For a blank or zero period, only the listed epoch counts. The timing error is
  σ_T = sqrt(σ_E² + (k·σ_P)²), from the table's error columns, with a blank error taken as 0.
  A predicted time is usable only if 3σ_T < 25 h (the epoch jitter Salinas et al. validated
  over) and it falls inside an orbit of that sector:
  - orbit start and end are from `singles/data/tess_orbit_times_latest.csv`;
  - those are UTC, converted to BTJD as JD_UTC + 69.184 s − 2457000;
  - the barycentric offset (≤ 8.3 min) is ignored at this stage.
- **A TOI stays in the sample** if it has at least one usable predicted time in a sector with
  a light curve. Every dropped TOI keeps its row in `data/sample.csv`, with a `drop_reason`.
- **Rows are TOIs.** A star with several TOIs gives several rows. The bootstrap resamples TOIs.
- **Easy negatives (optional):** long-period EBs from Prša et al. (2022, ApJS 258, 16), using
  the local copy at the repo root, `data/raw/tess_eb_prsa2022.tsv` (gitignored, so it is only in
  the main checkout). Their period is > 10 d, or blank, or 0 in
  that table, and they are not TOIs. The same light-curve and coverage rules apply, with the
  primary eclipse epoch as the epoch. These are a separate comparison only, never mixed into
  the FP/FA negatives.

## 2. Windows and score (Step 4)

- **Preprocessing.** Exactly as `computeparams/param_inference_core.py` in
  `helemysm/singletransiters` at commit `b7b9d49`:
  - PDCSAP flux (SAP if PDCSAP is missing);
  - quality == 0, NaNs removed, binned to 30 min;
  - a 99.5 h window of 200 points centred on the predicted time;
  - linear interpolation, with gaps beyond the data filled with the median;
  - divided by the window median, minus 1;
  - two input channels: scaled time and flux.

  No detrending.
- **Window validity**, checked on the binned light curve before interpolation: at least 50% of
  the 30 min bins within ±T14/2 of the predicted time hold data, and at least 50% of the 200
  window bins hold data. T14 is the table's `Duration (hours)`. Invalid windows are dropped
  and counted.
- **Per-TOI score:** the maximum of the window scores over that TOI's valid windows.
  - Secondary, reported alongside: the mean over the TOI's windows, and the main AUC with each
    TOI's first valid window only. These show whether TOIs with more windows gain from the
    maximum.
- **Window score: not yet defined.** Step 0 found the Salinas CNN gives no transit score
  (`STEP0_CNN_OUTPUTS.md`). The score is fixed by Jon's choice between the options listed there.
  It goes into this file as a dated change before any window is scored.

## 3. Blends (Step 3)

A negative (FP or FA) is tagged a **blend** if either rule holds.

**(C) Comment rule.** The TOI `Comments` field matches, case-insensitively, any of these
regular expressions:

```
\bNEBs?\b
nearby\s+(eclipsing\s+binar|EB\b|star|source|target|neighbou?r)
\bcentroid\s+(offset|shift|motion)
\boff[- ]target
\bon\s+(a\s+|the\s+)?(nearby\s+)?(neighbou?r|companion|nearby)
\bblend(ed|ing)?\b
\bcontaminat(ing|ion|ed)\s+(star|source|EB)
\bsignal\s+(is\s+)?(from|on)\s+(TIC|a\s+nearby|the\s+nearby|neighbou?r)
\bassociated\s+with\s+TIC
```

The comment is split into segments at `;` and `.` followed by a space. A match does not count
if its segment also contains any of these:
- `clear` (as in "cleared", "NEB clear");
- `\bno\b`, `\bnot\b`;
- `ruled out`, `excluded`, `unlikely`.

Every match, with its segment and whether a negation removed it, is written to
`data/blend_comment_matches.csv`. The matches are not edited by hand. A by-eye check of the
matches is reported as a count of the ones that look wrong, and it does not change the tags.

**(G) Gaia rule.** A Gaia DR3 source other than the target, within 42″ of the target, with
flux_n / (flux_t + flux_n) > depth, where:
- flux comes from G-band magnitudes, flux = 10^(−0.4 G);
- depth is the table's `Depth (ppm)` × 10⁻⁶;
- the target's Gaia DR3 source is the TIC v8 `GAIA` (DR2) ID mapped to DR3 through
  `gaiadr3.dr2_neighbourhood`, taking the match with the smallest angular distance;
- if that fails, the brightest DR3 source within 3″ of the TIC position at epoch 2016.0;
- if both fail, the TOI is flagged `gaia_unmatched` and counts as not meeting (G);
- neighbours with no G magnitude are skipped.

The Gaia rule says only that a neighbour *could* make the signal, so it overcounts.

**Diagnostic.** Rule (G) is also applied to the positives and its pass rate reported. This
gives the rate at which (G) flags real planets.

**Blend fraction and ceiling:**
- f = n_blend / n_negatives, with a 95% Wilson interval, reported for rule (C) alone and for
  (C or G);
- the ceiling is AUC_max = 1 − f/2, with its interval from the Wilson limits on f;
- this assumes blends are indistinguishable from planets (AUC 0.5 on that subset) and all other
  negatives are perfectly separable, as in Pagés Navarrete (arXiv 2610.42009).

## 4. Metrics (Step 4)

- **ROC AUC** (`sklearn.metrics.roc_auc_score`), with a 95% CI from 2,000 bootstrap resamples:
  - positives and negatives are resampled separately, with replacement;
  - percentile interval;
  - numpy seed 20261124.
- **Comparisons** (positives are always the CP/KP sample):
  - all FP/FA;
  - non-blend FP/FA;
  - blend FP/FA;
  - easy EB negatives, if built.

  The non-blend and blend sets are reported for both blend definitions.
- **Confounder checks:**
  - depth and Tmag distributions of positives and negatives: medians, quartiles and a
    two-sample KS test;
  - the main AUC (CP/KP vs all FP/FA) repeated in fixed depth bins of the table depth: < 1,000,
    1,000–3,000, 3,000–10,000 and > 10,000 ppm. A bin with fewer than 10 of either class is
    reported but not interpreted.
- **Overlap with Salinas et al.:**
  - the CNN was trained only on synthetic light curves, so it has no real training TOIs;
  - its 120 real validation TOIs (CP/KP, P > 10 d) are not listed in the paper;
  - if the list can be had, the main AUC is repeated with them removed;
  - if not, the note says the overlap could not be checked;
  - if option (b) is picked, its training set must share no TIC with the test sample, and that
    is checked and reported.
- **Sanity check:**
  - 50 out-of-transit windows from positives, drawn with seed 20261124 from the positives' light
    curves;
  - each centred at least 3 days from any predicted transit of any TOI on that star, and
    valid under the same rule;
  - their score distribution is reported next to the negatives' scores;
  - it is a check only, and it is not used to tune anything.

## 5. Changes

(none yet)
