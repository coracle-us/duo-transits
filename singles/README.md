# Single-transit community candidates (pilot)

This is the follow-on project to the duotransit work; see `docs/problem_set_survey.md` for why.
The aim is to revisit stale, never-promoted single-transit CTOIs with all TESS data to date.
Each candidate gets one of three outcomes:
- a period or alias set, from repeat transits;
- a source location, from difference imaging at the event;
- "still single", with the period range the data excludes.

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 1 | `scripts/build_ctoi_singles.py` | Filters the ExoFOP CTOI table to stale singles, then drops confirmed hosts and stars that have a TOI. Output: `data/ctoi_singles.csv`. |
| 2 | `scripts/select_pilot.py` | Brightest targets with TESS sectors after the CTOI's last update, and no community period on ExoFOP. Output: `data/pilot_targets.csv`. |
| 3 | `scripts/pilot.py` | Per target: confirm the event, search all data for repeats, infer candidate periods, map the excluded periods, and make difference images at the event and repeats. |

Step 3 reuses the validated tools in `../scripts/`: `alias_pipeline.test_window`,
`diff_image.analyse` and `fetch_data._fetch_one`.

**Target list (step 1, 2026-10-11):**
- 450 open single/duo CTOIs whose owners are not actively working them.
- 392 CTOIs on 378 stars remain after the filters. Dropped: 46 with missing parameters,
  8 whose host star already has a confirmed planet, and 4 stars that already have a TOI.
- Stale entries are real: CTOI 1042868.01 is HD 88986 b, confirmed at 146.05 d by
  Heidari et al. 2024 (2024A&A...681A..55H), but the CTOI was never updated.

## Validation of single-event mode (`results/validation/`)

The test set is the 26 Hawthorn duotransit candidates with known periods. The pipeline is
given **only the first transit** (Tc1, depth, duration) and must find the rest.

| Check | Result |
|---|---|
| Catalogued event recovered | 26/26 |
| Second transit (Tc2) found as a repeat | 23/26 |
| Repeats lying on the true ephemeris | 69/71; the other 2 are TOI-6695, which has known TTVs |
| True period in the candidate-period list | **26/26** |
| Uniquely recovered | **19/26** |
| Wrong unique answers | **0** |

**Change made during validation (disclosed).** The period-consistency step first projected
the period from the nearest repeat, which left 4 targets with no candidates. Repeat times are
only known to the grid step (±D/8), and projecting to distant repeats multiplied that error.
It now fits a linear ephemeris to all events for each trial cycle assignment, with residuals
required below D/3. Detection thresholds were not changed: they are the alias pipeline's
frozen values.

**Caveats:**
- The validation stars are duotransit hosts, which are bright and have deep, already-confirmed
  signals. Stale CTOIs may be shallower or noisier.
- Difference imaging was validated earlier on TOI-5542 b and TOI-2589 b
  (`results/validation_notes.md`).

## Pilot results (20 brightest eligible targets, 2026-10-11)

Run: `pilot.py` then `pilot.py --redo-diffs` (after the difference-image reference-point fix),
then `summarize_pilot.py`, giving `results/pilot/classified.csv`. Classification rules were
fixed in code before the outcomes were seen.

| Outcome | Count | Targets |
|---|---|---|
| Still single, no repeat in new sectors | 12 | 28–84% of log-period space 5–1500 d now excluded |
| **New duotransit** | 1 | **TIC 147035057** (see below) |
| Period narrowed | 1 | TIC 51079186: 21.17 / 10.59 / 5.3 d; weak dips (event SNR 5.9), unvetted |
| Short-period periodic signal, not a long-period planet | 4 | 53847727, 89620053, 140298355, 235060821 |
| Catalogued event not seen in our data | 1 | TIC 152283270 (10 ppt CTOI) |

Notes on the short-period signals:
- **53847727** is the known EB P = 3.0113 d (Prša+2022). We rediscovered it independently:
  3.011 d, secondary eclipse at phase 0.499.
- **89620053:** the event is 51″ off-target at 113σ, on a G = 11.6 neighbour. Likely a nearby
  eclipsing binary.
- **140298355** and **235060821:** their 6–7 dips also align at about 1 d, but that alignment is
  weak because few dips plus short trial periods line up by chance.

Since the pilot, the target builder also drops known TESS EBs (Prša+2022): 24 CTOIs on 19 stars.

**Go/no-go gate: FAILED as specified.** Only 7/20 targets (35%) have a conclusive difference
image, against the ≥60% bar. Picking the brightest stars first selected mostly shallow events
(0.3–0.7 ppt), whose difference images in 30-min/10-min FFIs are too weak. The full list has a
median depth of 3 ppt, so selecting by expected difference-image SNR (depth × brightness) rather
than by brightness alone should do much better. That is a hypothesis to test, not a result.

### TIC 147035057: new duotransit candidate

- **Star:** TYC 5963-01913-1, Tmag 8.49, Teff 6480 K, R* ≈ 2.0 R☉.
- **CTOI history:** 147035057.01, posted 2021-02 by the Single Transit Finder group as a single
  transit in S33. ExoFOP has no follow-up of any kind.
- **Second transit found:** S7 (BTJD 1515.065), in addition to the catalogued S33 event (BTJD 2224.051).
  - Depths agree: χ² 0.08/1. The two depths are 3.19 and 3.28 ppt; the model depth is 3.6 ppt.
  - QLP, TESS-SPOC and SPOC 2-min depths agree to within 4%.
  - No centroid shift; nothing near momentum dumps or segment edges.
- **Period:** ΔT = 708.976 d, and 15 aliases survive all TESS data to date (37.3–709 d). This is
  a warm-Jupiter-size candidate (CTOI radius 11.8 R⊕).
- **Difference images:**
  - S33 is on target (2.9″ ± 1.0, SNR 15).
  - S7 is 6.3″ ± 1.1 (0.3 px). The pre-registered rule therefore classifies the source as
    "mixed". In both images the missing light is a compact blob on the target pixel, and the
    G = 13.5 neighbours at 17″ and 24″ (each needing a ~25% eclipse) are not favoured. The
    Cycle-1 control TOI-2589 showed 3″ offsets, so a 6″ systematic in S7 is plausible but not proven.
- **Not in Hawthorn+24**, although it matches their S1–13 + S27–39 two-transit pattern.
- **Files:** `results/vetting/147035057_n1.*`, `results/pilot/147035057*`.

#### Update, same day: period solved; host is RV-variable

- **S108 eclipse (BTJD 4313.689) at the predicted time for P = 37.315 d.** S108 was released weeks ago and has no SPOC/QLP light curve, so we used our own TESScut photometry (`scripts/tesscut_photometry.py`).
  - The same code recovers the known S33 eclipse at 3.21 ± 0.07 ppt.
  - The S108 depth is constant across 3/9/21-pixel apertures (2.77/2.78/2.69 ppt) and the background is flat.
  - The difference image puts the source 2.2″ ± 0.8″ from the target.
  - Of the four aliases that predicted an S108 eclipse, only n = 19 shows one; n = 14 falls on flat data.
- **Ephemeris:** P = 37.31497 ± 0.00003 d, T0 = BTJD 1515.0661 ± 0.0011. The three events have O−C ≤ 1.3 min (χ² 0.97, 1 dof). Files: `results/vetting/147035057_ephemeris.json` and `147035057_next_eclipses.csv`; next eclipse well placed from Paranal is 2027-01-20 03:11 UTC ± 3.4 min.
- **No secondary eclipse at phase 0.5:** −0.05 ± 0.07 ppt.
- **Gaia DR3: strongly RV-variable** (rv_chisq_pvalue 0, rv_renormalised_gof 99.8, rv_amplitude_robust 49.8 km/s, 27 epochs; RUWE 1.12; no NSS orbit). This points to a massive, probably stellar companion. An eccentric grazing EB is favoured; a planet or BD in a binary is not excluded. RV monitoring is needed.
- **Write-up and draft email for observers:** `writeups/tic147035057/`.

## Batch 2 (20 targets ordered by expected SNR, depth 0.02–5%; 2026-10-11)

`select_pilot.py --by-snr --exclude data/pilot_targets.csv --out batch2_targets`, then
`pilot.py data/batch2_targets.csv batch2` and `summarize_pilot.py batch2`. The selection excluded 15 CTOIs with
implausible listed depths (e.g. 95–99%, a units or entry error) or depths below 200 ppm.

- **Gate:** difference images conclusive for **65%** of targets, which passes the 60% bar (the pilot reached 35%).
- **Period verdicts:**

  | Verdict | Count |
  |---|---|
  | Unique period | 3 |
  | Narrowed | 2 |
  | Many aliases (new or extended duotransits) | 3 |
  | Repeats inconsistent | 1 |
  | Catalogued event not seen | 4 |
  | Still single | 7 |

- **Source verdicts:** 8 on-target, 3 off-target, 2 mixed, 7 inconclusive.

**Unique periods, all eclipsing binaries.** Gaia DR3 SB1 orbits match the TESS periods independently:

| TIC | TESS period | Gaia SB1 period | e | K1 | Notes |
|---|---|---|---|---|---|
| 231969683 | 21.978 d (5 eclipses) | 21.969 ± 0.005 d | 0.36 | 32.6 km/s | CTOI by steuer, 2021 |
| 230130556 | 22.259 d (33 eclipses since S15) | 22.245 ± 0.024 d | 0.32 | 12.3 km/s | companion ≈ 0.2 M☉; CTOI by vogl, 2022 |
| 6610653 | 16.134 d | — | — | — | off-target in all 4 difference images (13–18″, 13–32σ); nearby EB |

**Lesson.** Gaia binarity flags (NSS orbit, RV variability, RUWE) are cheap and catch these before
any light-curve work, so `scripts/gaia_flags.py` now computes them for every target
(`data/gaia_binarity_flags.csv`). TIC 147035057 shows the same pattern (RV-variable F star, P of
tens of days) without a published Gaia orbit, which strengthens the eclipsing-binary
interpretation there.

**Still to vet from batch 2:**
- TIC 270450687 (4% deep, ΔT 680 d)
- TIC 206666063 (narrowed to 231.7/115.9/77.2/57.9/29.0 d; S2 difference image off 6.4″)
- TIC 30836645 (ΔT 2525 d)
- TIC 458419988 (ΔT 734.6 d)
- TIC 100776118 (also Hawthorn; Gaia SB1)

**Outreach log.** 2026-10-11: email with the TIC 147035057 write-up sent by Jon to S. Gill and D. Bayliss (Warwick), asking them to forward it to the CTOI submitter. Draft: `writeups/tic147035057/email_draft.md`. The earlier Hawthorn-catalogue email (`../results/outreach/`) has not been sent as of this date.

2026-10-11 (later): Jon sent the TIC 159614322 write-up to the Planet Hunters TESS team at Oxford: S. Aigrain, cc O. Barragán. The addresses come from the Oxford Physics people pages. The email gave P = 89.9651 or 179.9302 d and the next transit, 2026-12-28 06:42 UTC ± 45 min. Draft: `writeups/tic159614322/email_draft.md`.

The TIC 147035057 write-up was revised the same day (the eclipse is grazing, so the companion is likely stellar). The version Warwick received said "≈1.1 R_Jup, larger if grazing". The email itself already called the companion "probably a low-mass star". A short follow-up with the revised PDF is optional.

## Batch 3 (20 Gaia-clean targets, depth 0.2–1.2%, by expected SNR; 2026-10-11)

New since batch 2: pilot.py measures unprocessed recent sectors (≥ 100) with our own TESScut
photometry, and selection keeps only Gaia-clean stars (`gaia_binarity_flags.csv`: no NSS
solution, not RV-variable, RUWE ≤ 1.4).

- **Gate:** difference images conclusive for 50% of targets.
- **Unique periods (4), all off-target:**

  | TIC | Period | Source evidence | Depths (ppt) |
  |---|---|---|---|
  | 197736350 | 29.55 d | 77–104″ offset, 7–23σ | 2.2–5.3 (inconsistent) |
  | 290686569 | 29.16 d | 19–21″ offset, 6–20σ | |
  | 356868485 | 27.91 d (23 repeats) | 19–20″ offset, 23–42σ | alternate 1.4× and 0.6× (EB primary/secondary pattern) |
  | 427344083 | 82.73 d | 13–14″ offset, 6–9σ | 7.3–12 (inconsistent) |

- **New duotransits:**
  - 58946313 (ΔT 719.8 d) is on target (0.6″, 3.1″) but the two depths differ at 8σ in the same SPOC 2-min product (11.1 vs 7.5 ppt). EB-like.
  - 814205011 (ΔT 1499 d) is off-target (36–38″, 40–70σ).
- **Other:** 12 still single; 1 event not recovered; 1 with inconsistent repeats.

### Running tally after 60 targets

Every repeating signal found so far is a binary star or a blend:
- **Gaia-confirmed or Gaia-flagged EBs on the target:** 53847727, 231969683, 230130556, 458419988,
  270450687, 30836645, 206666063, and probably 147035057.
- **Off-target (nearby-star EB):** 89620053, 6610653, 431945694, 77437461, 197736350, 290686569,
  356868485, 427344083, 814205011, 190989264, 173639849.

No planet-like periodic signal has survived. The useful product is a classification of these
stale "single-transit planet candidates" (EB, blend, or still single with excluded period ranges).

## Recheck list: tentative detections in unprocessed sectors (own photometry)

| TIC | What | Why tentative | Recheck when |
|---|---|---|---|
| 268534931 (also Hawthorn+24, open, N_alias 31) | Possible third transit, S107, BTJD 4268.870, depth ~6–7 ppt, from our own TESScut photometry. If real, the aliases narrow from 31 to 4: 92.60 / 46.30 / 30.87 / 23.15 d. Host: Gaia-clean K dwarf (Teff 5179 K, R* 0.78 R☉); candidate ~6 R⊕ | See below | SPOC/TESS-SPOC/QLP products for S107 appear on MAST; rerun `pilot.py` (it prefers pipeline products over OWN) |
| 139726266 | Second transit, S107 (own photometry), BTJD 4289.85, depth scale 1.07 ± 0.07 vs S32. K dwarf (R* 0.75), so Rp ≈ 2.8 R⊕; 27 aliases 2109.0/n | Repeat seen only in our photometry; diff images SNR < 5; a G = 16.1 star at 6.7″ would need a 23% eclipse | Same: pipeline products for S107 |
| 168702709 | Second transit, S107 (own photometry), BTJD 4268.9 | Partial: a gap through mid-transit; 2-min centroid shift 3.9σ in S31 | Same |

Why 268534931 is tentative:
- Only an ingress plus the bottom are covered; the egress falls in a 5-h gap at a momentum dump (+7.2 h).
- The flux level jumps by ~5 ppt across the gap.
- S107 has strong scattered-light trends (±7 ppt).
- The difference image is 9.5″ off at 4.5σ.

Plot: `results/vetting/268534931_s107.png`.


## Full run (all remaining targets; 2026-10-11)

`run_remaining.py`, 3 workers:
- 278 targets: 246 analysed, 32 skipped because ExoFOP already carries a community period (`results/full/skipped.csv`);
- `classify_final.py` over all runs gives **306 stars**.

**Classification:**

| Class | Stars |
|---|---|
| Single | 193 (median 52% of log P excluded) |
| Event not recovered | 34 |
| Repeats inconsistent | 23 |
| EB on target (Gaia) | 21 |
| Off-target EB (firm) | 18 |
| Duotransit / periodic, hand-vetted | 17 |

The hand vetting is in `results/vetting/leads_full_run.md`.

**Two leads survive:**
- **TIC 159614322:** S14 + S41; ≈ 0.96 R_J; P = 719.657/n, n ∈ {1, 2, 4, 5, 8}.
- **TIC 229055790:** a PHT candidate with CHIRON spectra already on ExoFOP; S1 + S68; ≈ 3.8 R⊕; 38 aliases ≥ 26.7 d.

Both are:
- on target;
- consistent in depth across epochs and products;
- re-detected with independent code.

Two more (139726266, 168702709) join the recheck list.

## Partial-transit sweep (all 272 stars with a recovered event; 2026-10-11)

**Why.** TIC 159614322's third transit was a partial one: the ingress in S21, with the egress lost
in a downlink gap. The repeat search never counts such events, because it needs baseline on both
sides. `scripts/partial_sweep.py` therefore tests every data-segment edge in every product.

**Rules.** They were fixed in advance and are listed in the script docstring:
- a one-sided step model at the right depth;
- red-noise SNR ≥ 5 and a flat transit floor;
- detection in ≥ 2 products;
- a brightening control (steps up) run with identical rules.

Three implementation details were changed after the positive control (159614322 S21) failed and
before any other star was examined. Details are in the docstring.

**Result: no new convincing partial transits.**

| Count | Two-product dips | Two-product brightenings (control) |
|---|---|---|
| All stars | 48 | 37 |
| Without TIC 236770234 (a variable star, 21 dips) | 27 | 24 |

As a population, gap-edge dips do not exceed the false-alarm rate. One reason is that SPOC 2-min
and TESS-SPOC share pixels, so the 2-product test is weaker than intended. Individual hits fall
into three groups:
- **Shared timing.** Many hits coincide with hits on other stars at the same orbit edge
  (S7, S45, S52 and S57 starts and ends), so they are systematics.
- **Already-classified binaries.** Many are on stars already classed as EBs or off-target blends
  (53847727, 258799245, 384882500, 82550215, 350274116 and 427344083). These are more eclipses of
  known binaries.
- **Noisy stars.** Some are on stars that also show brightenings at similar SNR.

The two left over were inspected by hand. Both are start-of-sector ramps, and QLP trims exactly
those cadences:
- **126278294:** S44 start (2514.6); if real, it would leave P = 55.5 or 27.8 d.
- **253553257:** S13 start.

The positive control is the only clean, high-SNR, on-target, alias-selecting partial transit:
**TIC 159614322 S21 (SNR 24)**.

Two other stars show dips at the same S21 orbit end, but on a different detector (camera 4 /
CCD 1, while 159614322 is on camera 3 / CCD 3):
- 356868485 is a real eclipse of a known off-target EB;
- 236770234 is the variable star.

So there is no sign of a field-wide systematic at that edge.

**Outputs:** `results/partial_sweep/` (one JSON per star, `summary.csv`) and
`results/partial_sweep_stdout.txt`.

## New-sector re-runs (monthly, from 2026-10-11)

`scripts/new_sector_rerun.py`:
- finds stars observed in a released TESS sector (ended more than 7 d ago) that have no
  pipeline light curve in our cache;
- re-fetches from MAST and re-runs `pilot.run_target`;
- writes `results/rerun_<date>/` (`report.md`, `changes.csv`, JSONs);
- flags NEW/LOST repeats, period-list changes and event-verdict changes.

A scheduled task (`monthly-new-sector-rerun`) runs it on the 5th of each month at 10:00. The
task commits the results but does not vet them. `classify_final.py` includes `rerun_*`
directories, and the newest result wins. Hand-vetting outcomes go into `REJECT_RERUN` in
`classify_final.py`; unvetted re-run changes are provisional. The RNAAS note uses the frozen
`results/final_classification_2026-10-11_rnaas.csv`.

**Catch-up run 2026-10-11:** 45 stars had pending sectors (S100–108). Only 6 had new files on
MAST so far, all of them own TESScut photometry; 39 had no new data yet. TIC 152359413 also
comes from the 2-star dry run. All flags were vetted by hand, and **none is a real new event**:

| TIC | Flag | Verdict |
|---|---|---|
| 152359413 | 3 S106 "repeats" | Rejected: 2–11 d apart, no single P ≥ 5 d, scattered-light systematics |
| 228463329 | 3 S107 "repeats" | Rejected: 0.8–6.7 d apart, systematics |
| 98794593 | 1 S108 "repeat" (Tmag 7.7) | Rejected: post-gap ramp; S108 own photometry has ±3 ppt trends; diff image 19″ off |
| 100776118, 206666063 | more eclipses | Known Gaia binaries; the class is unchanged |

The three rejected stars are listed in `REJECT_RERUN`. The classification is unchanged from the
frozen copy.

**Lesson:** own TESScut photometry in scattered-light-heavy sectors (S106–108) produces
repeat-like dips at SNR 8–12. Treat OWN-only repeats as tentative until pipeline products confirm
them, as for the S107 recheck list.
