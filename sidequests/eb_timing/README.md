# Side quest A: eclipsing binaries with hidden third bodies (eclipse timing variations)

This is a scoping study, a 20-system southern pilot, a calibrated timing-error model, and a blind test against the
Borkovits+2026 northern CVZ results. It is separate from the duotransit work and writes only inside
`sidequests/eb_timing/`. Light curves are cached in `cache/lc/` (gitignored, 1.4 GB for 141 stars).

## Bottom line

- **The gap is real but narrow, and it is in the south.** The northern CVZ has been done twice by the
  Borkovits/Mitnyan group, for >3,500 EBs: S14–60 in Mitnyan+2024, and S14–86 in Borkovits+2026, which
  found 168 triple candidates. I found **no systematic ETV study of southern CVZ EBs** using the
  7.4-yr baseline TESS now has (S1–98, up to ~43 sectors). The only partial coverage is
  Marcadon & Prša 2024 (CVZ EBs from both poles, S1–60, ~100 timed, 5 solved, all northern) and
  Esmer & Daylan 2025 (152 detached EBs all-sky, few eclipses each). Details are in `notes/eb_literature.md`.
  The obvious risk is that the Borkovits group does the southern CVZ next.
- **The pipeline is validated.** Against Borkovits+2026, which uses the same baseline as ours, we recover P3 to
  ≤0.5% and the LTTE amplitude to ≤4% for all three LTTE systems. P3 is within 0.4% of the
  published value for all 7 validation triples, including two dynamically dominated southern CVZ triples
  (P_out = 77 d and 331 d). Three method changes were made during validation and are listed below.
- **Error bars are now honest.** A per-system O−C noise model replaces the formal errors. It has a per-sector common jitter,
  an antisymmetric per-sector "spot" term, and a wild cluster bootstrap over sectors. Against 21 B26
  group-L1 systems (42 P3 and A pulls), the pull spread drops from **RMS 7.2 / robust sd 2.8** (formal) to
  **RMS 1.3 / robust sd 1.0**. The calibration factor needed is only **1.02**. Held-out checks agree: 15 B26 D1
  periods give RMS 0.98, and the 5 validation systems (8 pulls) give RMS 1.04.
- **Northern CVZ blind test (S14–86, B26's own baseline).** Of the 35 B26 "certain" (L1) triples, 13 cannot be
  tested from MAST light curves (≤2 sectors; B26 made their own FFI photometry). Of the **22 testable, P3 is
  recovered for 21 (95%) and 16 (73%) are flagged automatically**. The 5 unflagged ones fail our
  primary/secondary-consistency rule. Among **49 random NCVZ EBs that B26 did not list**, the pre-registered
  thresholds flag 6. **4 are pipeline false alarms (8%)**, and 2 look like real long-period signals missing from
  B26. A third, TIC 233193944, is real but unflagged. Three extra rules written after vetting those false alarms
  (P3 ≥ 5 P1, A ≥ 3 × sector jitter, season jackknife) remove all 4. That is in-sample, so the 0% is optimistic. The rules
  keep 16/22 L1, 11/20 D1, and all 4 southern candidates.
- **Pilot (top 20 by eclipse count, pre-registered):** 1 convincing new LTTE orbit, 3 large
  long-period ETVs with only about one cycle covered, 6 dominated by spots or anticorrelated signals, 7 flat, and
  3 too noisy. All 4 third-body signals were re-derived by an independent timing method
  (correlation 0.90–0.998), and none appears in the 27 triple/ETV papers checked or in its SIMBAD
  reference list.
  - **TIC 32000625** (TYC 9155-22-1, P = 0.7018 d): **P3 = 749 ± 5 d, A_LTTE = 371 ± 14 s, e3 = 0.43 ± 0.06**
    (calibrated noise-model errors),
    with 3.6 cycles covered and primary and secondary moving together. f(m3) = 0.13 M☉, so **m3,min ≈ 0.9–1.1 M☉**
    for m12 = 1.4–2.0 M☉. Gaia RUWE 11.5 supports an astrometric companion.
  - TIC 279323859 (P3 ≈ 2,700 d ≈ span, A = 789 ± 5 s), TIC 278861826 (P3 = 2,100 ± 370 d, A = 620 ± 140 s),
    and TIC 294327312 = HD 56910 (Am-type, Gaia SB2; P3 = 2,240 ± 17 d, A = 692 ± 31 s; already flagged
    as an astrometric/proper-motion binary in Hipparcos-based lists, but no orbit has been found). These are real ETVs with only one
    cycle seen, and **tess-point (0.9.5) shows no southern-CVZ sectors in S99–134**, so TESS cannot pin these periods before ~2028.
- **Verdict: still a modest go, not a big one, and the NCVZ test strengthens it.** The pipeline now has measured
  recall (73% flagged and 95% period recovery on testable B26 L1 triples), a measured false-alarm rate (8% with
  the pre-registered thresholds), and calibrated errors. The yield is good (4/20 with third-body ETVs, 1/20 with
  a well-constrained orbit), but the result is a candidate list. There is one new practical limit: **faint stars
  (T ≳ 13.5) often have no MAST FFI light curves after S26**. A full southern run therefore needs either our own TESScut
  photometry (slow: ~1 min and 10–40 MB per star-sector) or a magnitude cut. Turning
  it into a paper needs dynamical/LTTE modelling of the best systems and ideally RVs. Recommended outlet:
  **an RNAAS note** listing southern CVZ triple candidates from a full run on the ~114 selected
  southern CVZ EBs (or all 1,191 that pass the cuts). Contact the Borkovits group **before** scaling up.
  AAVSO VSX is not the right venue: it catalogues variability, not ETV orbits.

## Data

- **Catalogue:** Prša et al. 2022 TESS EB catalogue (VizieR J/ApJS/258/16; 4,584 rows, 4,578 TICs),
  from `scripts/get_catalog.py` → `data/prsa2022_tess_ebs.csv`. It has **no** eclipse times or ETV flags.
- **Coverage:** tess-point on-silicon test for S1–107 (`scripts/coverage.py`). Good-data days
  per sector are measured from SPOC 2-min files already in `data/raw` (`scripts/sector_table.py`).
  S97/S98 are ~42 good days each.
- **Sample** (`scripts/make_sample.py` → `data/sample_all.csv`, `data/sample_top200.csv`). The cuts were fixed
  before any timing: 0.4 ≤ P ≤ 30 d, morphology 0–0.6 (detached/semi-detached), primary
  depth ≥ 1%, span ≥ 2,000 d. Systems are ranked by expected primary-eclipse count, N_ecl = Σ good_days / P.
  1,191 EBs pass. **The top 200 are all southern** (81 in the SCVZ proper). The northern CVZ falls out
  because its S14–86 span is ≤1,979 d (and it is already covered by Borkovits+2026). Top-200 range: N_ecl
  135–1,971, span ≥ 2,176 d, median 29 sectors, P = 0.41–7.6 d, median Tmag 11.1.
  **13 of the top 200** are in `data/known_multiples.csv`. That list holds 481 TICs harvested from the LaTeX of
  27 triple/quadruple/ETV papers by `scripts/known_triples.py`.
- **Light curves** (`scripts/fetch_lc.py`): one product per sector, SPOC 120 s > TESS-SPOC > QLP,
  in the same format as `scripts/fetch_data.py`, which is reused for column conversion. Files go to `cache/lc/`.
  `fetch_data._fetch_one` itself was not used because it downloads every product including 20-s data
  into `data/raw`. 27 stars and 1,072 sector files were fetched, with 0 failures. Loading uses
  `scripts/alias_pipeline.load_lightcurves` (imported, pointed at the cache): quality == 0, sector-normalised,
  QLP segment starts trimmed.

## Method (`scripts/etvlib.py`, `scripts/run_etv.py`)

1. **Seed ephemeris:** catalogue P and T0. A PDM-refined period is also tried (`etvlib.refine_period`).
   The catalogue seed is kept unless the PDM seed times >10% more primaries or has less than half the O−C scatter.
2. **Template:** a full-phase template (median in ~90-s phase bins, light 3-bin smoothing) built from the
   2-min data. It includes ellipsoidal and reflection shape. For FFI points the template is **averaged over the exposure**
   (sub-sampled box), so 200/600/1800-s data are fitted with a matching template.
3. **Per eclipse:** the window is ±1.5 eclipse half-widths (capped at 0.3 P). The model is a + b(t−T) + c·(template(t−dt) − 1),
   where c absorbs dilution differences between products and sectors. For each trial dt the linear parameters are solved,
   giving a profile χ²(dt): a coarse grid, then a fine grid, then a parabolic vertex. σ comes from the curvature with χ² scaled to reduced χ² = 1.
   Each side of the minimum must have ≥60% of the expected cadences.
4. **Iteration:** the template is rebuilt from points shifted by the running-median O−C of *their own* eclipse
   type, and the eclipses are re-timed with that running median as a prior (up to 3 passes). This keeps hour-scale ETVs and
   apsidal motion from smearing the template.
5. **Cleaning and errors:** 5σ outliers from a 15-eclipse running median are flagged. Errors are scaled by the robust
   scatter of consecutive-eclipse differences within a sector. This measures white noise only. Red noise from
   spots remains, so formal P3/A errors are optimistic (see Caveats).
6. **O−C models** are fitted jointly to primaries and secondaries, with separate zero points: linear; quadratic; "aps"
   (secondary slope differs, i.e. apsidal motion); and LTTE (Irwin 1952, common to both eclipse types) on top
   of the linear *or* the quadratic ephemeris. The LTTE search uses a circular-orbit periodogram over P3 = 20 d to 1.5 × span,
   and the best 4 peaks are refined with eccentric fits from 36 starts each. Models are compared with a BIC in which errors are
   inflated so that the best model has reduced χ² = 1 ("ΔBICc"). Diagnostics: the correlation ρ between the primary and secondary
   sinusoids at P3 (+1 means LTTE or dynamical; −1 means spots or apsidal motion), and the per-sector O'Connell asymmetry as a spot proxy.
7. **Per-product offsets** (`scripts/summarize.py`): the mean LTTE residual per product (SPOC/TESS-SPOC/QLP).
   All are consistent with zero, e.g. TIC 279323859 QLP +0.2 ± 2.2 s and TIC 32000625 QLP −5 ± 7 s. No
   cadence-dependent bias is seen.
8. **Gaia/TIC** (`scripts/gaia_check.py`, reusing `singles/scripts/gaia_flags.py`): RUWE, non_single_star,
   NSS orbit, RV variability, TIC mass → `data/gaia_tic_flags.csv`. **Literature** (`scripts/lit_check.py`):
   SIMBAD id/type/bibcodes and VSX → `data/lit_check.csv`.
9. **Vetting by eye:** every O−C figure (`figs/<TIC>_oc.png`) was inspected. The verdicts are in `data/manual_classification.csv`.
10. **Independent re-derivation** (`scripts/independent_check.py`): this shares only the loader. Each sector is folded on
    the linear ephemeris, and the vertex of a parabola through the binned fold minimum gives that sector's O−C. The result is compared
    with the pipeline's sector means → `results/independent_check.csv`, `figs/<TIC>_indep.png`.

M_min uses m12 = 2 × TIC mass (twin assumption; the TIC mass is a single-star estimate from the blended Teff),
or 2.0 M☉ when the TIC has no mass. It is meaningful only for LTTE-dominated ETVs.

## Validation (`scripts/validate.py` → `data/validation_table.csv`)

Published references: M24 = Marcadon & Prša 2024 (2-min data, S14–60); B26 = Borkovits+2026 (FFI, S14–86,
the **same baseline as ours** for the NCVZ); B20 = Borkovits+2020 (photodynamical). The A values for B26 are
computed from their a·sin i, e and ω.

| TIC | ref | P3 pub (d) | P3 ours (d) | A pub (s) | A ours (s) | e pub | e ours | note |
|---|---|---|---|---|---|---|---|---|
| 219900027 | B26 | 1631 | 1634.7 | 134 | 139 | 0.49 | 0.57 | B26: "very uncertain". M24 gives 1159 d / 113 s on a shorter baseline |
| 229771234 | B26 / M24 | 359.1 / 359.8 | 359.3 | 181 / 180 | 176 | 0.50 / 0.34 | 0.50 | |
| 424461577 | B26 | 1013.7 | 1008.3 | 462 | 466 | 0.48 | 0.47 | M24 (S14–60): 996 d / 317 s / 0.17 |
| 356896561 | B26 / M24 | 475.3 / 477.5 | 476.5 | – / 51 | 51 | 0.52 / 0.005 | 0.33 | GZ Dra; LTTE + dynamical |
| 259006185 | B26 / M24 | 200 / 199.6 | 200.2 | – / 39 | 33 | 0.12 / 0.08 | 0.26 | dynamical + apsidal motion (ρ = 0.16) |
| 167692429 | B20 | 331.5 | 329.9 | – | – | 0.56 | 0.87 | dynamical ETV, ±3 h, eclipse depth varies; σ = 170 s |
| 220397947 | B20 | 77.08 | 76.8 | – | – | 0.23 | 0.50 | dynamical (2+1+1) |

P3 is recovered for 7/7 systems. A agrees with B26 to 1–4% for all three pure-LTTE systems. For e and for the
dynamically perturbed systems only P3 is meaningful: our model is pure LTTE, and in those systems the e and A values it
returns are not physical. For TIC 220397947 it gives f(m3) = 8.5 M☉, an obvious flag. **Rule:** if P3/P1 is small or
the implied f(m3) is absurd, the ETV is dynamical and the LTTE mass is not reported.

**Changes made during validation and the pilot.** All runs reported here use the final code.
1. The P3 upper bound was cut from 3× to 1.5× the data span, and when the BIC-preferred variant sits at that bound, the
   other LTTE variant (linear vs quadratic ephemeris) is used instead. For TIC 219900027 the quadratic+LTTE fit ran off to
   P3 = 5,937 d (a degeneracy). The linear+LTTE fit gives 1,635 d, which agrees with B26.
2. Templates are aligned per eclipse type with an iterated prior. Without this, TIC 167692429 (apsidal motion, ±3 h
   ETVs) smeared the template and gave σ ≈ 3 h.
3. The PDM period seed was added in the pilot. The catalogue P of TIC 140661916 is off by 4×10⁻⁵ d, which drifts by 4 h over 7 yr.
   Using PDM unconditionally broke TIC 167692429, hence the "keep the catalogue seed unless PDM is clearly better" rule.

## Error model (`scripts/ocfit.py`, `scripts/noise_fit.py`, `scripts/calibrate.py`)

The timing step is unchanged. The O−C fit is redone on the saved eclipse times (`results/oc*/<TIC>_oc.csv`).

- **Model:** O−C = ephemeris (linear or quadratic, separate p/s zero points) [+ LTTE] + **k·δ_s** + **γ_s** + ε.
  - δ_s is a free per-sector *antisymmetric* offset (k = +1 primary, −1 secondary). It is used only for sectors with ≥3 of each
    eclipse type. Starspots shift primary and secondary minima in opposite directions, while a light-time effect shifts
    both equally. So δ_s absorbs spot noise without absorbing an orbit.
  - γ_s is a per-sector *common* random offset, N(0, τ²), implemented as ridge rows in a GLS fit.
  - ε is white noise: the formal σ times a scale factor s_w ≥ 1.
  - τ and s_w are estimated by the method of moments from the residuals of the best model, then held fixed for every model.
- **Fit:** the linear parameters are profiled with a precomputed QR projector. P3 comes from a projected sinusoid
  periodogram (20 d to 1.5 × span), followed by eccentric refinement of the best 4 peaks from 8 starts each.
- **Uncertainties:** a wild cluster bootstrap. Each sector's residual block gets a random sign (B = 150), and the model is
  refitted from the best solution. σ is half the 16–84% interval. This keeps any within-sector correlation and the sector
  jitter. It does **not** include the choice between a linear and a quadratic ephemeris, which matters when P3 ≳ span.
- **Calibration:** pulls are (ours − B26)/√(σ² + σ_B26²) for P3 and A of the 21 B26 L1 systems whose P3 we recover.

| errors | n | RMS | robust sd | frac. within 1σ | frac. within 2σ |
|---|---|---|---|---|---|
| formal (run_etv, Jacobian × √χ²_r) | 42 | 7.24 | 2.84 | 31% | 50% |
| noise model + bootstrap, raw | 42 | 1.32 | 1.02 | 62% | 88% |
| calibrated (× 1.02) | 42 | 1.30 | 1.00 | 64% | 88% |
| held out: B26 D1 P3 (formal → calibrated) | 15 | 1.27 → 0.98 | 1.58 → 1.05 | 47% → 73% | 87% → 93% |
| held out: 5 validation systems vs B26, full baseline | 8 | 1.04 | 0.75 | 88% | 88% |

  The tails (RMS 1.3 against a robust sd of 1.0) come from systems where B26 add terms we do not fit: a quadratic, a cubic, or a
  fourth body. Figure: `figs/pulls.png`. CAL = 1.02 is applied to all southern numbers (`results/noise_fit.csv`).
- **Classification, pre-registered ("v1"):** a system is flagged if ΔBIC(LTTE vs the best of linear/quadratic/apsidal) > 30,
  A/σ_A > 4, P3 is not at the bound, and, when both eclipse types are present, ρ(pri, sec) ≥ 0.8 and A_pri/A_sec is within 1.5×.
  The class is "LTTE" if P3/span ≤ 0.7, otherwise "LTTE-long". The SNR, ρ and amplitude-ratio thresholds were tuned on
  the southern pilot (my by-eye labels) **before** the NCVZ test. "v2" adds three rules written after vetting the NCVZ
  false alarms (below). On the southern pilot, v1 and v2 both flag exactly the 4 by-eye candidates and nothing else.

## Precision ranking (`scripts/rank_precision.py` → `data/sample_ranked_precision.csv`)

The 1,191 selected EBs are re-ranked by the expected LTTE amplitude error, σ_A ≈ σ_T √(2/N_eff).
- The per-eclipse timing error is σ_T = K σ_1h(Tmag) √(W · 1 h)/(2δ) for a V-shaped eclipse, with depth δ and
  duration W from the catalogue.
- K = 1.09 is calibrated on our 27 timed systems. Measured and predicted σ_T scatter by 0.35 dex.
- The pre-registered `sample_top200.csv` and `pilot_list.csv` are untouched.
- Only 73 of the old top 200 stay in the new top 200. The new top 200 are brighter (median T = 9.5) and have longer
  periods (median P = 1.9 d). The old pilot's 20 systems rank between 8 and 801.
- Spots are not modelled, because the catalogue has no Teff.

## Northern CVZ test (`scripts/b26_tables.py`, `scripts/ncvz_sets.py`, `scripts/score_ncvz.py`)

- **Reference:** `data/b26_triples.csv` holds the B26 LaTeX tables: L1 35, L2 28, L3 64, D1 31, D2 10, for 168 TICs. Each row
  has P2, a·sin i, e, ω and τ with errors, plus any second (fourth-body) solution.
- **Sets fixed before timing:**
  - Recall set (commit 689e0e6): all 35 L1. The 31 D1 were added in commit 5e0420a, after the first download showed that 13 L1 stars have
    ≤2 sectors of MAST light curves, but still before any timing.
  - False-alarm set: 50 of the 66 Prša EBs in the NCVZ (β > 78°) that pass the southern cuts, have T ≤ 15, and are not in
    `known_multiples.csv`. They were drawn with `default_rng(20261012)`.
  - Seed ephemerides use B26 P1, because 30 of the 35 L1 are not in Prša.
  - Timing was restricted to S14–86 (`run_etv.py --max-sector=86`). Outputs: `results/oc_ncvz/`, `figs/ncvz/`,
    `results/noise_ncvz.csv`.
- **Recall** (`results/calibration.csv`, `results/ncvz_summary.json`):

| B26 group | n | testable (≥10 sectors) | P3 within 10% | flagged v1 | flagged v2 |
|---|---|---|---|---|---|
| L1 (certain pure LTTE) | 35 | 22 | **21 (95%)** | **16 (73%)** | 16 |
| D1 (certain LTTE + dynamical) | 31 | 20 | 15 (75%) | 12 | 11 |

  - The 5 L1 systems with the right P3 but no flag fail the primary/secondary consistency rule: three have ρ = 0.14, −0.24
    and 0.76, and two have amplitude ratios of 0.50 and 2.27. B26 note that for most systems they fitted *averaged* p/s ETVs,
    because the secondaries are noisier or spotted.
  - The one testable L1 miss, TIC 232632890, has τ ≈ 30 min of sector jitter.
  - The D1 misses are mostly compact systems where our pure-LTTE model locks onto a harmonic.
  - Thresholds were **not** retuned on the recall set.
- **False alarms** (`results/ncvz_fa_results.csv`; by-eye verdicts in `data/ncvz_fa_vetting.csv`):
  - TIC 377307592 turned out to be a B26 L3 system, because `known_triples.py` misses 15 B26 TICs written with LaTeX
    superscripts. It is excluded, which leaves 49.

| | flagged | pipeline false alarms | real signals not in B26 |
|---|---|---|---|
| v1 (pre-registered) | 6 / 49 | **4 (8%)**: 362104293 (P3 ≈ P1 sampling alias in a 29-d eccentric EB), 377192658 (one bad season), 237116051 and 199716496 (spots) | 2: 353854078, 441804888 |
| v2 (rules below) | 1 / 49 | 0 (in-sample) | 1: 353854078 |

  - One more real long-period signal, TIC 233193944 (a common ±5 min O−C), is not flagged because its best LTTE sits at the
    P3 bound.
  - So about 3/49 (6%) of the NCVZ EBs that B26 did not list show convincing long-period common ETVs. B26's catalogue is not
    complete for P3 ≳ span, and these three are worth mentioning to them.
- **Changes made during this test (disclosed):**
  1. v2 rule: P3 ≥ 5 P1 (hierarchical stability). This removes 362104293.
  2. v2 rule: A ≥ 3τ (the orbit must exceed the sector-to-sector jitter). This removes 237116051 (A/τ = 1.3) and 199716496 (2.3).
     Every recovered L1 system has A/τ ≥ 6.8.
  3. v2 rule: a season jackknife. A common non-linear O−C (quadratic or LTTE vs linear) must keep ΔBIC > 30 with each observing
     season left out in turn. This removes 377192658, but also 441804888, which looks real, so v2 is deliberately conservative.
     A first version required the LTTE to beat a parabola, which killed every LTTE-long system; it was replaced before the
     final run.
  4. D1 was added to the recall set before timing (see above). The FA list was not changed, and the generator still reproduces it.
  5. Because v2 was derived from these 49 systems, its false-alarm rate (0/49) is in-sample. Use the v1 rate of 8% as
     the honest estimate until a fresh FA sample is run.
- **Bolt-on for the B26 minima** (`scripts/compare_minima.py`):
  - Input: their times of minima as a CSV with columns TIC, T (BJD or BTJD), type and err.
  - It matches each minimum to our eclipse by type and epoch, and reports the offset, RMS difference and pull RMS.
  - It also fits their times with the same noise model, for a like-for-like P3/A comparison.
  - Self-test on GZ Dra: our own times plus a 30 s offset and 20 s noise, treated as "external", give back the 30 s offset and
    20 s scatter (`results/minima_compare_selftest.csv`).
  - Usage: `python scripts/compare_minima.py --ext=b26_minima.csv --oc=results/oc_ncvz`.

## Pilot results (`results/pilot_summary.csv`, `results/noise_fit.csv`, figures `figs/<TIC>_oc.png`)

Selection was fixed before timing (`data/pilot_list.csv`): the first 20 rows of the top-200 after removing
validation systems and known multiples. None of ranks 1–20 needed removing. In total there are 33,608 eclipse times. The
median primary timing error is 15 s, and the range is 5–270 s.

| class | n | TICs |
|---|---|---|
| **LTTE candidate** (orbit well covered) | 1 | 32000625 |
| **LTTE, long period** (≈1 cycle; P3 lower limit) | 3 | 279323859, 278861826, 294327312 |
| spots (O−C tracks the O'Connell asymmetry) | 4 | 140661916, 33864821, 25078924, 220474209 |
| anticorrelated pri/sec (spots or apsidal motion; M dwarfs) | 2 | 38586438, 308451526 |
| flat (no ETV; best-fit A ≤ 4 s for 6 of 7) | 7 | 393491149, 349832824, 350622185, 350443417, 271596290, 374859582, 364302118 |
| too noisy (σ_pri > 100 s) | 3 | 349911034, 167123437, 349642881 |

Third-body signals:

| TIC | P (d) | P3 (d) | P3/span | A_LTTE (s) | e3 | f(m3) (M☉) | m3,min (M☉) | Gaia RUWE | independent check (corr pri / sec) | known? |
|---|---|---|---|---|---|---|---|---|---|---|
| 32000625 | 0.7018 | 749 ± 5 | 0.28 | 371 ± 14 | 0.43 ± 0.06 | 0.13 | 0.89–1.08 (m12 1.4–2.0) | 11.5 | 0.985 / 0.948 | no ETV/triple paper; SIMBAD refs: RAVE, catalogues, Prša+22 |
| 279323859 | 1.0800 | 2,728 ± 12 (lin variant 2,692) | 1.0 | 789 ± 5 | 0.63 ± 0.01 | 0.15 | ≳1.4 (m12 3.0) | 1.07 | 0.998 / 0.856 | no |
| 278861826 | 1.0676 | 2,110 ± 370 (earlier fits 2,580 to >3,900) | 0.8 | 620 ± 140 | 0.47 ± 0.14 | 0.07 | ≳1.0 (m12 2.5) | 2.25 | 0.998 / 0.997 | no |
| 294327312 (HD 56910) | 0.9185 (true 1.837, Gaia SB2) | 2,240 ± 17 (quad variant 2,170) | 0.82 | 692 ± 31 | 0.63 ± 0.07 | 0.13 | ≳1.5 (m12 3.6) | 5.8 | 0.90 / – | flagged as an astrometric/proper-motion binary (SIMBAD refs: Makarov & Kaplan 2005, Frankowski+2007, Kervella+2019, Brandt 2021); no ETV or orbit paper found |

The errors come from the calibrated noise model (`results/noise_fit.csv`, CAL = 1.02). The automatic classification, in
both its v1 and v2 forms, flags exactly these four: 32000625 as "LTTE" and the other three as "LTTE-long". For the three
long-period systems, the bootstrap error leaves out two things: the choice between a linear and a quadratic ephemeris,
and the fact that only about one cycle is covered. Read their P3 as "about the baseline or longer", not as a measurement. TIC 32000625's m3,min ≈ 1 M☉ for a K-type
(Teff ≈ 4,960 K) system implies the tertiary contributes substantial third light. That fits its shallow 8% eclipses
for a semi-detached geometry, but this was not checked further.

## Caveats

- **Errors:** the formal errors (`results/oc/*_fit.json`, `pilot_summary.csv`) are 3–7× too small. Use the noise-model
  errors in `results/noise_fit.csv`, which are calibrated to B26 (robust pull sd 1.0, RMS 1.3). They still leave out the
  choice of ephemeris model when P3 ≳ span.
- **Data availability:** 13 of the 35 B26 L1 stars and 11 of the 31 D1 stars have ≤2 sectors of MAST light curves. B26's own FFI
  photometry reaches fainter stars than the SPOC, TESS-SPOC and QLP products do.
- **Model scope:** pure LTTE plus a quadratic term. There is no dynamical (three-body) term and no apsidal-motion model. Compact
  triples with P3 ≲ 100× P1 come out with correct periods but wrong e and A.
- **Catalogue periods** can be half the true period (TIC 294327312, 271596290 and 25078924 have VSX/Gaia periods twice the
  catalogue value). The O−C is still valid, because both eclipse types shift together, but labels like "primary" are mixed.
- **TIC mass:** m12 is a crude twin assumption, and 5 of the 20 pilot stars have no TIC mass.
- **Literature check** used SIMBAD bibliographies, VSX, and TIC IDs harvested from 27 papers. It did not include the TESS
  Ten Thousand catalogue's free-text "pETV" comments or a full ADS full-text search. ADS was not available.
- The top-200 ranking favours short-period near-contact systems (median P ≈ 1.1 d), which are the most spot-prone.
  The precision ranking (`data/sample_ranked_precision.csv`) fixes this for any future run.
- The v2 false-alarm rate is in-sample. An unbiased number needs a fresh random sample, northern or southern.

## Recommendation

- **Go, at modest scale, with a phone call first.** The southern CVZ (S1–98, 7.4 yr) is the one place where
  a systematic TESS ETV search has not been published. The pipeline reproduces the reference (B26) solutions,
  and a blind 20-system pilot gave one clean new LTTE triple and three large long-period ETVs. Scaled to the
  1,191 EBs that pass the cuts (~17 GB of light curves, ~1–2 CPU-days at ~2 min per system), the NCVZ-like
  rate (~4%) suggests **~40–60 candidates, ~10–20 with well-covered orbits**.
- **Before scaling:**
  1. Email the Borkovits/Mitnyan group to ask whether a southern CVZ paper is in preparation, and offer the candidate list or a
     collaboration instead of racing them. Also mention the 3 NCVZ long-period signals that are not in B26 (353854078,
     441804888, 233193944).
  2. Done: calibrated noise model, precision ranking, measured recall and false-alarm rate.
  3. Decide what to do about faint stars: a magnitude cut (T ≲ 13.5), or our own TESScut photometry.
- **Expected performance, measured in the NCVZ:** about 73% recall for B26-quality LTTE triples, and an 8% false-alarm rate before
  vetting by eye (v1). Most of those false alarms are removed by the v2 rules plus vetting by eye.
- **Outlet:** RNAAS ("Third-body candidates among TESS southern CVZ eclipsing binaries from 7 years of
  eclipse timing", one table). A full A&A/MNRAS catalogue paper needs dynamical modelling and is the Borkovits group's
  territory. AAVSO VSX is not suitable. For TIC 32000625 alone, a short note or a message to the group is enough.
- **Not worth it** if the Borkovits group confirms a southern CVZ study is in progress. Then hand over the 4 candidates and stop.

## Files

- `notes/eb_literature.md`: literature check. `notes/arxiv_search_raw.txt`: raw arXiv API listing.
- `data/`: catalogue, coverage, sector table, samples, known multiples, pilot/validation lists,
  Gaia/TIC flags, literature check, manual classification, validation table.
- `results/oc/<TIC>_oc.csv`: every eclipse time (BTJD), errors, flags and residuals. `results/oc/<TIC>_fit.json`: fits.
- `results/pilot_summary.csv`, `results/validation_summary.csv`, `results/independent_check.csv`.
- `figs/<TIC>_oc.png` (vetting figures) and `figs/<TIC>_indep.png` (independent check).
- Error model and NCVZ test: `results/noise_fit.csv` (south + validation, calibrated), `results/noise_ncvz.csv`
  (NCVZ, raw), `results/calibration*.csv/json`, `results/ncvz_summary.json`, `results/ncvz_fa_results.csv`,
  `figs/pulls.png`, `figs/ncvz/<TIC>_oc.png`, `data/b26_triples.csv`, `data/ncvz_*.csv`,
  `data/sample_ranked_precision.csv`.
- Reproduce: `get_catalog → sector_table → coverage → known_triples → make_sample → fetch_lc TICs →
  run_etv TICs → gaia_check / lit_check TICs → validate → summarize data/pilot_list.csv results/pilot_summary.csv →
  independent_check TICs`. South noise model: `noise_fit --cal=1.0206 TICs`; then `rank_precision`.
  NCVZ: `b26_tables → ncvz_sets → fetch_lc --max-sector=86 TICs → run_etv --ephem=sidequests/eb_timing/data/ncvz_ephem.csv
  --max-sector=86 --tag=ncvz TICs → noise_fit --oc=results/oc_ncvz --out=results/noise_ncvz.csv TICs → calibrate →
  score_ncvz → plot_pulls`.
