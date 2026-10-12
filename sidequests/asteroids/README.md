# Side quest B: asteroid rotation periods from TESS full-frame images (scoping + pilot)

This is separate from the duotransit work. It writes only inside `sidequests/asteroids/`. Downloads
are in `cache/` (gitignored: 3.4 GB of TESScut cutouts, the LCDB zip and SBDB/Horizons tables);
all of it can be re-fetched by the scripts.

## Bottom line

- **Verdict: dead end for Jupiter Trojans and Hildas.** The pipeline works, but only where the
  periods are already known. No new periods came out of the pilot.
- **The method works at V ≲ 16.5 and fails at V ≳ 17.5.** Across 22 LCDB U=3 validation asteroids,
  12 passed the automatic acceptance rule. **All 12 recovered the published period**, to a median
  0.06% and at most 0.33%, over 5–164 h. None was wrong. Completeness falls from 7 of 9 at V < 16
  to 2 of 7 at V 17–18, and 0 of 1 above 18.
- **Pilot (30 Jupiter Trojans crossing Sector 72, selected before any photometry):**
  - 3 passed and were accepted by eye: (588) Achilles 7.307 h, (4063) Euforbo 8.849 h and
    (659) Nestor 15.98 h. All three match their LCDB U=3 periods to 0.04% or better.
  - A fourth, (5123) Cynus at 9.898 h, is clearly real by eye and matches LCDB to 0.15%, but the
    contamination cut rejected it.
  - **0 of the 10 pilot Trojans that lack an LCDB period got a period.** All 10 objects at
    V 18–19, which is where the unmeasured Trojans are, gave noise. Their residual rms was about
    0.5 mag per 30-min bin, and 60–90% of frames were star-contaminated.
- **The gap does not overlap the reachable range.** Take every numbered Trojan and Hilda on
  silicon for at least 10 d in any sector from S56 to S107 with predicted V < 16.5. That is 47
  Trojans and 30 Hildas, and **every one already has an LCDB U ≥ 2 period.** Extending to
  V < 17.5 adds 11 Trojans and 5 Hildas without LCDB periods (`data/yield_estimate.csv`). The method
  succeeds on fewer than half of objects at that brightness, and several of those 16 are likely
  already in ZTF (Li et al. 2025) or Ďurech & Hanuš 2026.
- **The only open niche is main-belt slow rotators at V ≲ 16.5 in sectors after S48.** This pilot
  did not test it, and it is a crowded race against the Konkoly TSSYS EDR2 team and the NASA Ames
  `tess-asteroids` team. Not recommended for a hobby project. See the recommendation below.

## 1. Literature and the gap (`notes/asteroid_literature.md`)

- TESS Year 1 has been mined three times: Pál et al. 2020 (TSSYS-DR1, 9,912 light curves),
  McNeill et al. 2023 (3,492 reliable periods) and Vavilov & Carry 2025 (4,521).
- The Konkoly group holds a processed data set through S48 (EDR2) and is publishing Trojan and
  Hilda results from it.
- Ďurech & Hanuš 2026 (arXiv:2610.06082) give 247 Trojan and 157 Hilda spin states from nine surveys
  including TESS. ZTF adds 216 Trojan periods (Li et al. 2025). Gaia DR4 arrives on 2 Dec 2026.
- **LCDB is stale.** The public file is still the 2023-Oct-01 release.
- Our cross-match of numbered asteroids (`scripts/lcdb_coverage.py` → `data/lcdb_coverage.csv`)
  gives these fractions with an LCDB U ≥ 2 period:

  | population | numbered | with LCDB U ≥ 2 |
  |---|---|---|
  | Jupiter Trojans | 10,633 | 376 (3.5%) |
  | Hildas | 3,821 | 184 (4.8%) |
  | outer MB | 292,903 | 2.6% |

  By H, Trojans run from 99% (H < 10) to 71% (10–11), 31% (11–12), **7% (12–13)** and 1% (13–14).
  The gap is at H ≳ 11.5, which means **V ≈ 18–20** at TESS geometry.
- Chosen population: Jupiter Trojans, because they have the lowest coverage and a high
  slow-rotator fraction. The pilot sector is S72, which has the most Trojans on silicon of
  S56–S107 (1,132 with H ≤ 13.5 for ≥ 10 d; `data/sector_counts_trojan_H13.5.csv`). S72 is also
  after the end of EDR2.

## 2. Data and tools

- **Cutouts:** `astroquery.mast.Tesscut.download_cutouts(objectname="<number>;", moving_target=True,
  sector=S, size=15)`.
  - The object ID must use Horizons syntax: `"624;"` works, while `"624"`, `"(624)"` and
    `"624 Hektor"` fail.
  - Moving-target cutouts work up to at least S100. One 27-d, 200-s cutout of 15×15 px is about
    30–70 MB and takes 1–3 min.
  - **Limits we hit:**
    - Parallel requests trigger HTTP 429/502. `fetch_cutout` now retries with backoff, and we
      download serially.
    - The 55-d Sector 98 gave 502 on every attempt.
    - TESScut listed (1600) in S76 but returned an empty zip.
    - The service raised errors for several numbered objects (e.g. 5182, 3068, 503) in
      `get_sectors`.
    - An object that crosses CCDs comes back as separate files. We keep the largest.
    - The cutout keeps following the target into the collateral columns (see §3).
  - HDU 2 of the cutout carries the full-CCD SIP WCS. From it, the CCD position of every stamp is
    recovered: `corner = round(WCS(TGT_RA, TGT_DEC) − TGT_X/Y)`, integer to ≤ 0.5 px.
- **Ephemerides:** JPL Horizons from TESS (`location="500@-95"`) gives r, Δ, α and V at 3-h steps.
- **Target selection:** `scripts/select_targets.py`.
  - Two-body propagation of SBDB elements on a 0.5-d grid, plus tess-point for the given sector.
    Sector dates come from `cache/TESS_orbit_times.csv` (tess.mit.edu).
  - Checked against the Hektor cutout, the RA/Dec error is about 0.13° after 4 years of
    propagation. That is fine as a pre-filter.

## 3. Method (`scripts/mtphot.py`, `scripts/period.py`)

1. **Background:** per-frame background is the median of the faintest 40% of stamp pixels. This
   is the rule of `singles/scripts/tesscut_photometry.photometry`. That function writes into
   `data/raw/`, so its rule is reused but the function is not called.
2. **Static scene:** each stamp pixel is mapped to its CCD pixel. The static sky is the per-CCD-pixel
   median over all good frames in which the asteroid was more than 3.5 px away. It is then
   subtracted.
3. **Aperture photometry:** pixel centres within 1.6 px of the Horizons position, minus the median
   residual sky in a 3.5–7 px annulus.
4. **Frames kept:**
   - QUALITY = 0.
   - Stamp entirely on science pixels (cols 44–2091, rows 0–2047).
   - At least 8 static samples per aperture pixel.
   - Annulus rms below 2.5× the sector median (scattered light).
   - Not star-contaminated: static-scene flux in the aperture below 3× the asteroid's median flux.
5. **Light curve:** 30-min bins (median, at least 3 frames), light-time-corrected times, and
   reduced magnitude H(α) = m − 5 log10(rΔ). A linear phase-angle term is fitted together with the
   Fourier series.
6. **Period search:**
   - Lomb–Scargle over 0.05–24 c/d.
   - 4th-order Fourier fits with the period free, at P_LS and 2·P_LS.
   - **Choosing between P and 2P:** an F-test of the odd harmonics at 2·P_LS, scaled for red noise
     (time-averaging β). Significant odd harmonics mean the two halves differ, so the rotation is
     2·P_LS. Otherwise 2·P_LS is adopted by convention and flagged.
   - **False-alarm level:** 200 shuffles of 1-d blocks of the light curve. This is conservative for
     P < 1 d, because a short-period signal survives inside each block.
7. **Acceptance rule** (`scripts/classify.py`, set on validation and frozen for the pilot):
   - boot_fap ≤ 0.005, i.e. the observed peak beats all 200 shuffles.
   - Contaminated fraction below 0.4.
   - At least 2 cycles of the adopted period.
   - Then a by-eye decision (`data/pilot_manual.csv`).
8. **Reused by import:** `scripts/diff_image.py` `gaussian` (2-D PRF model), used for the
   stacked-difference-image centroid check. Not reused: `diff_image.gaia_field`. Its cone is one
   15-px cutout wide, while a track is about 500–1000 px long, so contamination is measured from
   the static scene itself instead. `alias_pipeline` and `independent_bls` are transit-specific;
   only their quality == 0 convention is kept.

**Changes made during validation (disclosed):**
- (a) Red-noise β was redefined, because the first version confused underestimated errors with red
  noise (Hektor).
- (b) An off-chip flag was added after (12198) ran into the collateral columns.
- (c) The block-bootstrap FAP and the acceptance rule were added after seeing the validation
  results.
- (d) Tried and **not adopted:** a soft-edged aperture re-centred on the stacked source. It lowered
  the rms for 5 of 7 objects but broke (698)'s period (`scripts/compare_aperture.py` →
  `data/aperture_comparison.csv`).
- (e) LCDB parse bug fixed: 601 rows with a `*` flag had been dropped, including (1143) Odysseus.

## 4. Validation (`data/validation_results_classified.csv`, plots in `plots/`)

The selection rule is in `scripts/pick_validation.py`: LCDB U=3, non-binary, amplitude ≥ 0.15 mag,
two random objects per period bin, plus a Trojan and a Hilda, in a 200-s sector (S56+), V ≤ 17.
Then a fainter V 16.4–18.2 set was declared before running.
- (1600) was dropped: TESScut returned an empty cutout.
- (1345) was dropped: S98 failed with 502.
- (16421) moved to its next sector after an empty S62 cutout.
- (624) Hektor was the development object, so it is not independent.

| object | S | V | published P (h) | recovered P (h) | Δ | amp (mag) | auto |
|---|---|---|---|---|---|---|---|
| (698) Ernestina | 67 | 15.6 | 5.0363 | 5.0366 | 0.01% | 0.58 | pass |
| (12198) 1980 PJ1 | 84 | 15.9 | 5.494 | 5.5011 | 0.13% | 0.82 | pass |
| (733) Mocia | 66 | 14.5 | 11.374 | 11.3756 | 0.01% | 0.32 | pass |
| (5) Astraea | 71 | 11.0 | 16.806 | 16.7952 | −0.06% | 0.29 | pass |
| (314) Rosalia | 92 | 15.5 | 20.465 | 20.4775 | 0.06% | 0.22 | pass |
| (17681) Tweedledum | 69 | 16.8 | 75.2 | 75.24 | 0.05% | 0.82 | pass |
| (460) Scania | 92 | 15.7 | 164.1 | 163.96 | −0.08% | 0.38 | pass |
| (624) Hektor (dev) | 59 | 14.5 | 6.924 | 6.9204 | −0.05% | 0.44 | pass |
| (4348) Poulydamas (Trojan) | 83 | 16.2 | 9.908 | 9.9151 | 0.07% | 0.30 | pass |
| (3670) Northcott | 72 | 17.8 | 6.992 | 6.9700 | −0.31% | 0.63 | pass |
| (13018) Geoffjames | 70 | 17.6 | 4.334 | 4.3330 | −0.02% | 0.59 | pass (P/2P by convention) |
| (5254) Ulysses | 71 | 16.4 | 28.72 | 28.815 | 0.33% | 0.33 | pass |
| (1218) Aster | 72 | 17.4 | 3.1581 | 3.1614 | 0.10% | 0.38 | fail (boot FAP 0.010), correct |
| (12) Victoria | 72 | 12.8 | 8.660 | (8.646 at P_LS) | — | 0.04 | fail: only 5 d on one CCD, low amplitude |
| (181) Eucharis | 80 | 13.5 | 52.23 | — | — | 0.06 | fail: 12 d, amplitude 0.06 |
| (16421) Roadrunner | 84 | 17.3 | 174 | — | — | — | fail: 11.5 d, 53% contaminated |
| (3395) Jitka | 71 | 16.7 | 18.29 | — | — | — | fail (43% contaminated) |
| (18582) 1997 XK9 | 71 | 17.4 | 114 | — | — | — | fail |
| (1513) Matra | 71 | 17.0 | 117.6 | (122.3) | 4% | — | fail (59% contaminated) |
| (2235) Vittore | 81 | 16.5 | 32.1 | — | — | — | fail (75% contaminated) |
| (4791) Iphidamas | 83 | 17.1 | 9.696 | — | — | — | fail |
| (8550) Hesiodos (Hilda) | 92 | 18.2 | 6.719 | — | — | 0.09 published | fail |

- 12 of 22 pass, and all 12 are correct. No false positive passed the rule.
- Formal period errors are 5–15× too small, judged against the LCDB values (which have their own
  errors). Expect ~0.3% real accuracy from one sector.
- The 30-min-bin residual rms is about 0.04 mag at V < 16, 0.07 at V 16–17 and about 0.18 at
  V 17–18 (`data/success_vs_V.csv`, `plots/success_vs_V.png`). Star residuals dominate faint
  objects: the median contaminated fraction rises from 0.10 to 0.36 to 0.65.

## 5. Pilot: 30 Jupiter Trojans in Sector 72 (`data/pilot_vs_lcdb.csv`)

Selection rule (`scripts/pick_pilot.py`, fixed before photometry): ≥ 15 d on silicon in S72; the 20
brightest by predicted V (15.2–17.2), plus 10 spread evenly over V 18–19. LCDB status was not used.
The LCDB periods of the sample were printed at selection time, so the by-eye vetting was not blind.

| | N | auto-pass | accepted by eye | matches LCDB | new |
|---|---|---|---|---|---|
| V < 16 | 3 | 2 | 2 (Achilles, Euforbo) | 2 | 0 |
| V 16–17 | 11 | 1 | 1 (Nestor), plus Cynus (contamination-rejected, correct) and Meges (marginal, correct) | 1 (+2) | 0 |
| V 17–18 | 6 | 0 | 0 | — | 0 |
| V 18–19 | 10 | 0 | 0 | — | 0 |
| **total** | **30** | **3** | **3 (+1 known)** | **3** | **0** |

- The 10 pilot objects with no LCDB entry are (4835) Asaeus, (22014) and 8 objects at V 18–19.
  **None yielded a credible period.** Asaeus shows slow variation, but it tracks orbit-start jumps
  and contaminated stretches, so it is not credible.
- Every plot was examined; verdicts and notes are in `data/pilot_manual.csv`.
- S72 is a crowded field: 40–90% of frames are star-contaminated for V > 16.5.
- Typical failure modes are day-scale systematics after gaps and at sector starts (scattered light),
  and residuals of stars the asteroid passes. These look like slow rotators, which is exactly the
  regime one would want TESS for.

## 6. Caveats

- Validation is small (22 objects) and the acceptance rule was tuned on it. The pilot is a single
  sector in a crowded field. A sparser sector would do somewhat better at V 16–17.
- The photometry is basic: hard aperture, median static scene, no PSF fitting and no per-frame
  re-centring. A better pipeline (e.g. the Ames `tess-asteroids` or Konkoly's) would push the limit
  fainter by perhaps 0.5–1 mag. That would not reach V 18–19 for 0.1–0.3 mag amplitudes in 30-min
  bins.
- The yield estimate uses Keplerian positions (about 0.1° error) and H–G magnitudes with G = 0.15
  (about ±0.3 mag), and only checks LCDB. Adding ZTF, Ďurech & Hanuš 2026 and Konkoly would shrink
  the remaining set further.

## 7. Recommendation

**No-go for Trojans and Hildas.** TESS FFIs reach V ≈ 16.5 reliably, and everything there already
has a period. The unmeasured objects at V 18–20 are beyond reach, and Rubin and Gaia DR4 will
take them.

The only defensible TESS asteroid project would be a **post-S48 main-belt slow-rotator catalogue**
(P ≈ 24–300 h, V ≲ 16.5). LCDB covers only 25–70% of outer/middle-belt objects at H 12–14 (`data/lcdb_coverage.csv`), and the
validation shows we recover 75–164 h rotators cleanly. It would need:
- tumbling and binary checks;
- a many-sector run (TESScut throughput is about 1 object per 2 min, serial);
- an ALCDEF deposit.

It runs straight into Konkoly's promised EDR2 catalogue and the Ames pipeline. If attempted, scope
it to slow rotators the 2023 LCDB lacks, cross-check against Ďurech & Hanuš 2026 and Gaia DR4
(December 2026), and publish in the **Minor Planet Bulletin**: free, ADS-indexed and the standard
venue, with light curves to **ALCDEF**. An RNAAS note fits only a single clean population-level
result, which this pilot does not have.

## Files

| path | what |
|---|---|
| `notes/asteroid_literature.md` | literature check, references |
| `scripts/lcdb_coverage.py` | LCDB × SBDB coverage by population and H → `data/lcdb_coverage.csv`, `data/lcdb_slow_rotators.csv` |
| `scripts/select_targets.py` | Kepler + tess-point crossing tables → `data/sector_counts_*.csv` |
| `scripts/pick_validation.py`, `pick_pilot.py` | pre-registered samples → `data/validation_candidates.csv`, `data/validation_list.csv`, `data/pilot_sample.csv`, `data/pilot_list.csv` |
| `scripts/mtphot.py` | TESScut moving-target photometry → `data/lightcurves/<num>_s<SSS>.csv.gz` |
| `scripts/period.py` | periodogram, Fourier fits, P/2P test, bootstrap FAP → `plots/<num>_s<SSS>.png` |
| `scripts/run_batch.py` | runs both on a list → `data/validation_results.csv`, `data/pilot_results.csv` |
| `scripts/classify.py` | acceptance rule → `*_classified.csv` |
| `scripts/compare_aperture.py` | aperture variants → `data/aperture_comparison.csv` |
| `scripts/summarize.py` | `data/success_vs_V.csv`, `data/pilot_vs_lcdb.csv`, `plots/success_vs_V.png` |
| `scripts/yield_estimate.py` | Trojans/Hildas V < 17.5 in S56–107 vs LCDB → `data/yield_estimate.csv` |
| `data/pilot_manual.csv` | by-eye verdicts for all 30 pilot objects |
| `data/validation_tesscut_failures.csv` | validation objects lost to TESScut |
