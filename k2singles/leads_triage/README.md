# Triage of the AAS82225 revision's blind-search leads (2026-11-24)

**Yield: 0 of 47 lead stars survive as an on-target periodic candidate.** No period is claimed here.

## Inputs

- The 41 stars of `k2singles/revision2/results/new_leads.csv`: the pre-registered blind search, 56 signals (18 linked
  families, 38 BLS detections).
- The 8 post-hoc linked-period detections not in Table 1 (`results/link_calibrated.csv`): 202085278, 204776782,
  247795097, 215067200, 205966706, 247756662, 211598816 and 212805678. Two of them (211598816, 212805678) are also in
  `new_leads.csv`. For each one, the top-support period is triaged.
- In all: 47 stars and 64 signals (`leads.py`).

## What ran

- **`catalog_lookup.py` → `catalog_lookup.csv`** (row queries only, no files downloaded):
  - Gaia DR3 from the ESA archive: RUWE, non_single_star, and neighbours within 60″.
  - VizieR: VSX type and Gaia DR3 variability class, for the stars whose literature census already has those rows.
- **`triage.py` → `signals.csv`** (one row per signal) and **`triage.csv`** (one row per star). The rules are in the
  docstring and were fixed before the full run. The checks:
  - **R:** is each event confirmed in a second reduction?
  - **C:** the star's own inverted-flux control statistic against the signal's (control-level if ≥ 0.8×).
  - **X:** a light-curve centroid test from the flux-weighted centroids in the K2 llc and SPOC/TESS-SPOC files.
  - **E:** EB tests: odd/even, a dip at phase 0.5, Rp > 2 R_J, the catalogue's "secondary eclipse" comment, and
    VSX/Gaia/TESS EB classes.
  - **Flags:** Gaia binary (RUWE > 1.4 or NSS) and known variability.
- The light curves are the ones the search used. Log: `triage_run.log`.

## Counts (per star: the best verdict among its signals)

| Verdict | Stars | Pre-registered only | Post-hoc only | Both |
|---|---|---|---|---|
| on-target candidate | **0** | 0 | 0 | 0 |
| unclear | 6 | 5 | 0 | 1 |
| off-target | 1 | 1 | 0 | 0 |
| EB | 19 | 16 | 2 | 1 |
| control-level | 7 | 3 | 4 | 0 |
| systematic | 14 | 14 | 0 | 0 |
| **total** | **47** | **39** | **6** | **2** |

Per signal (64): EB 26, systematic 23, control-level 8, unclear 6, off-target 1.

- **EB (19 stars):**
  - 7 have the catalogue's own "secondary eclipse" comment.
  - 2 have a VSX eclipsing type.
  - 2 have a dip at phase 0.5 (202085278: 460 ppm, 4.7σ; 247756662: 3.7σ). Both are post-hoc detections.
  - 1 has an odd/even difference (212332380, 4.2σ).
  - 8 rest on size alone (Rp > 2 R_J from the box depth and the K2 catalogue R★). The borderline case is
    211894612 at 2.2 R_J (see below).
- **Systematic (14):** in 13 the events are not reproduced in a second K2 reduction. In 246425172 one family member
  is confirmed nowhere.
  - Most are BLS hits of 1–4 h at P ≈ 22–37 d, i.e. two epochs in one K2 campaign.
  - 9 of the 22 systematic BLS signals include the catalogued single event among their epochs. The rest of the
    "period" is noise.
- **Control-level (7):**
  - 4 post-hoc: 204776782 (support 12 vs control 11), 205966706 (6 vs 10), 215067200 (6 vs 5), 247795097 (10 vs 8).
  - 3 BLS at SDE 7.9–8.3 against own-control 6.5–7.1, near the 7.79 threshold.
- **Off-target (1):** 245929407.
  - Its K2 centroid puts the source 6.1″ ± 0.5″ from the target (one epoch). A Gaia neighbour sits 4.8″ away
    (ΔG = +2.4).
  - The TESS centroid gives 2.5″ ± 0.6″, small for 21″ pixels.
  - **Confirmed by the pixel test (below):** the K2 C12 dip is on the neighbour Gaia DR3 2433468825703524480
    (G = 14.17), not on the target.

## The six unclear stars (none is a clean candidate)

| EPIC | Signal | Why unclear |
|---|---|---|
| 201132684 | BLS 10.06 d, 500 ppm | R pass, centroid on target, but RUWE 1.53. **Not new:** a known K2 planet host ("apparent multi-planet system"; Mayo+2018 table 3; Exominer++ 9.965 d). |
| 211498244 | family of 3 events, 5.8 %, 7.2 h | R pass, K2 and TESS centroids on target, Rp ≈ 0.74 R_J on an M dwarf (G 52-4). But **no single period fits all three events** (n_allowed = 0), and Gaia **RUWE = 10.0**. The two TESS events are 1043.30 d = 2 × 521.65 d apart, which is the Exominer++/SPOC 521.64 d; the K2 event is off that ephemeris. Not new as a period. The star has two more SNR > 50 events outside the family (K2 C18: 1.9 %, 14.7 h; TESS S71: 6.8 %, 3.9 h). Durations run 3.9–14.7 h with no common period. Worth a look, but not a candidate here. |
| 211598816 | post-hoc 2.972 d, 159 ppm, 8.1 h | A chemically peculiar rotating variable (VSX/Gaia ACV/roAp), RUWE 3.72; centroid not informative at this depth. |
| 212549089 | BLS 27.4 d, 1 h | One epoch carries the BLS peak; second reduction untestable; centroid excursion (47″); Gaia non_single_star = 2. |
| 212685467 | BLS 27.3 d, 1 h | One epoch carries the BLS peak; untestable. |
| 248912804 | BLS 27.2 d, 1 h | A 1-h box at 30-min cadence cannot be re-measured; its epochs include the catalogued 9.8-h event. |

**Closest to a candidate: 211894612 (EB by the size rule).**
- Three events: K2 C5, K2 C18 and TESS S44.
- 6 aliases from 55.7 to 111.4 d.
- K2 and TESS centroids on target, confirmed in pixels in K2 C5 and C18 (below). RUWE 1.006, and no Gaia
  neighbour within 60″.
- Rp = 2.2–2.8 R_J (family box depth to catalogued depth, R★ = 1.31 R☉). That points to a low-mass stellar
  companion of an F star (an EBLM-type system), not a planet.

## Pixel difference images (2026-11-24, with Jon's OK)

Two stars, the only ones where a difference image could change the outcome.

- **Data** (`pixels_download.py`, about 80 MB, in `data/pixels/` → `duo-transits-cache/leads_triage/`, outside git):
  - K2 long-cadence TPFs: 245929407 C12; 211894612 C5 and C18.
  - 15×15 TESScut cutouts: S69 (245929407) and S44 (211894612). Both TESS events are in FFI-only sectors.
- **Method** (`pixels_diff.py`, rules in its docstring):
  - **K2:** a per-pixel deficit image, as in `independent_211924561/s09_k2_pixels.py`. It gets a free Gaussian fit,
    then a fit with PSF templates at the Gaia positions of the target and its neighbours.
  - **TESS:** the `scripts/diff_image.py` difference image and Gaussian fit.
- **Output:** `pixels.csv`, `figs/pix_*.png`.

| Star | Data | Fitted source vs target | Template fit | Verdict |
|---|---|---|---|---|
| 245929407 | K2 C12 | 6.6″ ± 0.6″ (11.5σ) from the target's Gaia position; 5.8″ from the out-of-event image centroid | target 0.4σ (depth 0.6 %); neighbour at 4.8″, ΔG +2.44: **9.0σ**, implied depth ≈ 10 % | **off-target** |
| 245929407 | TESS S69 | 4.1″ ± 0.6″; the nearest Gaia star to the fit is the 4.8″ neighbour (1.7″ away) | – | unresolved (< 0.5 px), consistent with K2 |
| 211894612 | K2 C5 | 0.5″ ± 0.2″ (2.7σ) | target 5.9σ, implied depth 2.3 % (no neighbour in the stamp) | **on-target** |
| 211894612 | K2 C18 | 2.4″ ± 0.35″ from the WCS position, but 0.1″ from the out-of-event image centroid; the C18 WCS is off by 2.5″ | target 8.8σ, implied depth 3.1 % | **on-target** |
| 211894612 | TESS S44 | 2.6″ ± 0.4″; the nearest Gaia star to the fit is the target; target pixel diff/err = 72 | – | unresolved (< 0.5 px), on target in practice: no Gaia star within 60″ |

**What this changes:** no verdict.
- 245929407 stays off-target, now with the source named. Its K2 C12 event is an eclipse of a G = 14.2 star 4.8″
  away (ΔG ≈ 2.4, needing roughly 10 % or more depth there). That is most likely an eclipsing binary.
- 211894612 stays on target and EB by size. The yield remains 0 of 47.

**Two changes made while running, before anything was committed:**
- The TESS fit now starts at the peak of diff/err. `diff_image.py` starts at the largest absolute difference. In
  S44 that was a saturated G = 9.8 star on the cutout edge (diff/err ≈ 3), and the free fit converged there
  (151″ off). The target pixel has diff/err ≈ 72.
- A TESS offset that is significant but under half a pixel is labelled "unresolved", not "on-target". The
  shared `scripts/diff_image.py` is unchanged.
- TESS offsets of 2.6–4.1″ at 6–7σ also appear for the isolated 211894612. TESS's formal errors are below its
  systematic floor at this scale.

## What could not run, and limits

- **Pixel difference images: only for 245929407 and 211894612** (section above). For the other 45 stars, X is a
  light-curve (moment-centroid) test. It is informative for 34 of the 64 signals, and on target for 33 of those.
- **SPOC DV reports:** cached only for two Table 1 stars. None was fetched for the leads.
- **Odd/even and phase 0.5:** run only for signals with a single period: the BLS and post-hoc signals, and
  families with one allowed period. They can't run on two-event alias lists.
- **1-h BLS signals** (24, nearly all in K2 long cadence) cover about 2 cadences per epoch. 16 fail R and 8 are
  untestable.
- **The Rp > 2 R_J rule** uses box depths and catalogue radii. Eight EB verdicts rest on it alone.
- **The Gaia ESA archive is down from 2026-11-25 for DR4** (released 2026-12-02). Re-running `catalog_lookup.py`
  needs the archive up.

## Rule changes (all logged in commits)

1. **Fixes after the three-star test, before the full run** (commit `1c24521`):
   - The per-epoch σ is now β·σ_white/√n_in (depth/SNR was infinite at zero depth).
   - The odd/even and phase-0.5 σ include half the spread between products. K2 detrendings differ in depth by up to
     40%: 212805678 has 31 ppt in PDC against 44 ppt in K2SFF.
2. **After the first full run** (commit `b824f2e`): an implied centroid offset ≥ 5 px is an "excursion", not
   off-target. A source that far off cannot shift a moment centroid in proportion to the depth.
   - This moved 212549089 (47″, a one-epoch 1-h BLS) from off-target to unclear.
   - No other star's verdict changed.

## Reproduce

```
python k2singles/leads_triage/catalog_lookup.py   # network: ESA Gaia TAP, VizieR CfA mirror
python k2singles/leads_triage/triage.py           # offline; needs the symlinks below
python k2singles/leads_triage/pixels_download.py  # network: MAST (K2 TPFs, TESScut)
python k2singles/leads_triage/pixels_diff.py      # network: ESA Gaia TAP for the field stars
```

The search's data cache must be linked in (none of it is edited):
- `k2singles/revision2/data` → `duo-transits-cache/revision2`
- `data/raw/lightcurves` and `data/raw/tesscut` → `duo-transits-cache/raw/…`
- `k2singles/cache` → `duo-transits-cache/k2cache`
- `k2singles/leads_triage/data` → `duo-transits-cache/leads_triage` (the downloaded pixels)
