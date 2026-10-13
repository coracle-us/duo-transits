# Pre-registration: southern eclipse-timing run, v1

Written 2026-10-13. It must be committed **before** any southern target outside the 20-system pilot is timed.
The run itself waits for the Borkovits/Mitnyan reply (see `outreach/`). If they are already doing the southern CVZ, this
plan is not executed.

## 1. Frozen code

- **Freeze commit: `f3af3ac`.** All detection code is as in that commit. The detection and noise-model code is unchanged
  since `be2601b`. Between the two commits only `independent_check.py` changed: it gained `--oc/--out/--max-sector` options,
  with unchanged defaults. Git blob hashes at the freeze:

| file | blob |
|---|---|
| `sidequests/eb_timing/scripts/etvlib.py` | e8c71207094c0cce6e4c345c10112b129b4fd913 |
| `sidequests/eb_timing/scripts/run_etv.py` | e5be1f379fdcdb79f6adfe0bdf82f5013ea43432 |
| `sidequests/eb_timing/scripts/ocfit.py` (noise model + classification, "v2") | 0bd5cf34619972537d0d9885ffab1f0858c00e9b |
| `sidequests/eb_timing/scripts/noise_fit.py` | 7c9a0642547a0ce3639714be2ab117238616dd20 |
| `sidequests/eb_timing/scripts/fetch_lc.py` | 73e36048da90e85a9360d214b8693b29ecdb05dd |
| `sidequests/eb_timing/scripts/independent_check.py` | 1eef39b61edbc9d867aa7bd512066e82d922cb3c |
| `sidequests/eb_timing/scripts/candidate_lit.py` | a8318f89b819b421be62adea38ac02cc73e3acee |
| `sidequests/eb_timing/scripts/ads_check.py` | d346cd3fbe154b90260b803d41f1e10d73930e79 |
| `scripts/alias_pipeline.py` (shared loader, imported) | 557722d64a958ea64f2e2b919897143364c4f513 |
| `scripts/fetch_data.py` (shared, imported) | e6d915aeddf5652a79608a0ba2e6e65574324928 |

- The two shared modules live outside `sidequests/eb_timing/` and can be changed by other work. Before the run, check
  `git ls-tree HEAD` against the table above. If any blob differs, run from a `git worktree` checked out at `f3af3ac`.
- Environment: Python 3.12 venv, numpy 2.5.4, scipy 1.18.1, pandas 2.3.3, lightkurve 2.6.0, astropy 8.0.1.
- Calibration factor **CAL = 1.0206**: `noise_fit.py --cal=1.0206`, derived in `results/calibration_summary.json`
  from the B26 L1 pulls.

## 2. Targets (`data/southern_run_targets.csv`, made by `scripts/make_southern_targets.py`)

- **Source:** `data/sample_ranked_precision.csv`. These are the 1,191 Prša+2022 EBs with 0.4 ≤ P ≤ 30 d, morphology 0–0.6,
  primary depth ≥ 1%, and TESS span ≥ 2,000 d, ranked by the expected light-travel-time amplitude error.
- **Magnitude cut: Tmag ≤ 13.5.** In `data/lc_availability_bins.csv`, all 112 Prša stars fetched so far with Tmag ≤ 13.5
  have MAST light curves for at least 96% of their on-silicon sectors (median 100%). The 4 stars at Tmag 13.5–14 have
  a median of 77% (only 1 of 4 reaches 80%), and the 5 stars at 14–16 have 64–79%. Without our own FFI photometry,
  13.5 is the faintest limit the measurements support.
- **≥ 3 observing seasons** (on-silicon sector groups separated by > 100 d), the minimum for the season jackknife.
- **Result: 1,088 targets.** 109 are in the southern CVZ proper, 18 are pilot systems (the other 2 pilot stars are fainter than 13.5),
  and 37 are in `known_multiples.csv`. All are processed. Known multiples are labelled, not removed.
- **Cost:** 11,602 sector files, about **3.5 GB** at 0.31 MB/file, and about **40 CPU-hours** (8.3 s per file for timing,
  45 s per star for the noise model; both measured). With 8 parallel processes that is about 5–6 h wall time, plus about 3 h
  of MAST downloads with 4 workers.
- **Order:** by `rank_precision`, so a partial run is still a well-defined, precision-limited sample.

## 3. Procedure (per target, no manual tuning)

1. `fetch_lc.py` (≤ S107; SPOC 120 s > TESS-SPOC > QLP, one product per sector).
2. `run_etv.py` with the catalogue ephemeris. Outputs go to `results/oc_south/`, `figs/south/`.
3. `noise_fit.py --oc=results/oc_south --cal=1.0206 --out=results/noise_south.csv`.
4. No parameter, threshold or code change during the run. If a crash forces a code fix, it is committed separately,
   explained in the README, and **all** targets are re-run with the fixed code.

## 4. Detection rules (frozen "v2", `ocfit.py`)

- **Noise model:** O−C = ephemeris (linear or quadratic, separate primary/secondary zero points) [+ LTTE]
  + per-sector antisymmetric spot offset (sectors with ≥3 of each eclipse type) + per-sector common random offset N(0, τ²)
  + white noise with scale s_w ≥ 1. τ and s_w come from the best model and are fixed for all models.
- **Uncertainties:** wild cluster bootstrap over sectors (B = 150, seed 0), 16–84% half-width, × CAL.
- **Flag** (all must hold):
  - ΔBIC = BIC(best of linear/quadratic/apsidal) − BIC(LTTE) **> 30**
  - A / σ_A(raw bootstrap) **> 4**
  - P3 not at the search bound (P3 ≤ 0.98 × 1.5 × span; search 20 d to 1.5 × span)
  - with ≥10 primaries and ≥10 secondaries: ρ(pri, sec) **≥ 0.8** and 1/1.5 ≤ A_pri/A_sec **≤ 1.5**
  - **P3 ≥ 5 P1**
  - **A ≥ 3 τ**
  - **season jackknife:** a common non-linear O−C (quadratic or LTTE vs linear) keeps ΔBIC > 30 with each observing season
    (gap > 100 d) left out in turn.
- **Categories:**
  - **"LTTE"** (flagged, P3/span ≤ 0.7)
  - **"LTTE-long"** (flagged, P3/span > 0.7: period ≈ baseline or longer, not a measurement)
  - **"no-flag"**
  - The v1 class (without the last three rules) is also recorded, for comparison with the NCVZ tests.

## 5. Vetting protocol (every flagged system; also every "no-flag" system with ΔBIC > 300, as a check)

1. **By eye:** the O−C figure (`figs/south/<TIC>_oc.png`): fold and template, periodogram, O−C with sector means, residuals,
   primary−secondary, O'Connell. Verdicts go in `data/southern_vetting.csv`: real / spots / systematics / apsidal / ambiguous.
2. **Independent re-timing:** `independent_check.py` (shares only the loader; fold-vertex timing per sector).
   Required: correlation ≥ 0.9 with the pipeline sector means for the primaries, and peak-to-peak agreement within 20%.
3. **Season jackknife:** already part of v2. The per-season fits are recorded for the candidate table.
4. **Product/cadence check:** the mean LTTE residual per product (SPOC/TESS-SPOC/QLP) must be consistent with zero
   (`summarize.py`).
5. **Gaia:** `gaia_check.py` (RUWE, NSS, RV variability), reported, not used to reject.

## 6. Per-candidate literature check (every system that passes §5)

`candidate_lit.py` must be run, and its output read, for each candidate:
- SIMBAD (cross-identified object): identifiers and the full reference list. Any reference title about timing, minima,
  period change, LTTE or a third body is read.
- Tokovinin MSC (J/ApJS/235/6, 20″), WDS (B/wds, 10″), VSX.
- `known_multiples.csv` (27 papers) and the B26/Mitnyan tables.
- **ADS full text** (`full:` over the TIC and every SIMBAD identifier). **"Already timed" check (CATS-style):** any hit
  with O−C, minima, timing, period analysis or study, or ephemeris in title or highlight is read. Example: TIC 32000625's two 2023
  CATS conference abstracts, which plot an O−C but report no cyclic term.
- The **Kreiner O−C atlas** and BAV/IBVS/OEJV minima lists are checked for named or bright stars (V* designations, HD/HIP).
- **Reporting:** prior O−C or timing work is always cited. A candidate is "new" only if no prior paper reports a cyclic O−C,
  a third body, or a tertiary in MSC/WDS consistent with the ETV. Otherwise it is "known: confirmed/updated".

## 7. What gets reported

- **Candidate table columns:** TIC, other IDs, Tmag, P1, N_sectors, N_eclipses, span, class, P3 ± σ, A_LTTE ± σ, e3 ± σ,
  ω3, T3, f(m3), m3,min (with the stated m12 assumption), ΔBIC, A/σ_A, ρ, τ, jackknife min ΔBIC, independent-check
  correlation, Gaia RUWE, literature status (new / known, with citation), vetting note.
- **Counts:** targets processed, timed, flagged v1, flagged v2, passed vetting, new, known.
- **False-alarm rate,** quoted from the out-of-sample NCVZ test (FA set 2, `results/ncvz_fa2_summary.json`):
  the v2 pipeline false-alarm rate k/N with its Clopper–Pearson 95% interval. It is applied as "expected number of false
  flags ≈ rate × N_targets before vetting by eye". The in-sample FA set 1 rate (0/49) is **not** quoted as a rate.
- **Recall** is quoted from the NCVZ test: 16/22 testable B26 L1 flagged; 21/22 periods recovered.

## 8. Stopping rules and marginal cases

- **Success:** at least 5 vetted "LTTE" candidates (not only "LTTE-long") with literature status "new", or at least
  10 vetted flags of either class. Then the RNAAS note is written, and the list goes to the Borkovits group first.
- **Null result:** fewer than 5 vetted flags means a short null-result note (rates per 1,000 EBs, with the recall
  from the NCVZ test), or no publication, as Jon decides.
- **Marginal cases:** flagged but failing one §5 check, ΔBIC 30–100, or "ambiguous" by eye. These are listed in a separate
  "marginal" table with the failing check named. They are not counted as candidates and not promoted after the fact.
- **Abort:** if the dry run (§9) or the first 50 targets show a pipeline failure rate (crash / no O−C) above 20%, stop and
  re-plan instead of patching mid-run.
- No second pass with changed thresholds will be presented as the pre-registered result. Any v3 is a new, separately
  labelled analysis.

## 9. Dry run (done before committing this file)

The frozen code was re-run on 3 committed pilot systems: 32000625 (LTTE candidate), 350622185 (flat, deep eclipses)
and 393491149 (flat, near-contact). It was compared with the committed `results/oc/` and `results/noise_fit.csv` by
`scripts/dryrun_compare.py`, with output in `results/dryrun_compare.csv`. Result: **identical.** For all three systems, every eclipse time and its error are reproduced exactly (2,121, 1,476 and
2,561 eclipses; max |ΔT| = 0, max |Δσ| = 0). The noise-model P3, A, e, their σ, ΔBIC and class are identical, because the
bootstrap is seeded. The dry run used the code at `be2601b`, whose detection blobs are identical to `f3af3ac`.

## 10. False-alarm rate measured for this plan (NCVZ FA set 2, out of sample)

- **Sample:** `data/ncvz_fa2_list.csv`, fixed before timing (commit 685ef7a). It holds all 38 remaining Prša EBs in the
  Mitnyan/B26 footprint (≥ 8 Year-2 sectors) with the same cuts as FA set 1. Excluded: FA set 1, the recall sets, the known
  multiples, and every TIC in the B26/Mitnyan LaTeX (superscripts handled). Only 38 were eligible, not 50.
- **Results** (`results/ncvz_fa2_summary.json`, verdicts in `data/ncvz_fa2_vetting.csv`):
  - **v2 (the frozen rules): 1/38 flagged; pipeline false alarms 0/38 = 0% (Clopper–Pearson 95% CI 0–9.3%).** The one flag,
    TIC 255921197, is a clean, independently re-timed LTTE orbit that is not in the literature (see README).
  - v1: 2/38 flagged; pipeline false alarms 1/38 = 2.6% (95% CI 0.1–13.8%).
  - Combined with FA set 1 for v1 (both sets out-of-sample for the v1 rules): 5/87 = 5.7% (95% CI 1.9–12.9%).
- **Quoted rate for the southern run:** the v2 pipeline false-alarm rate **0/38 (95% upper limit 9.3%)**. With 1,088
  targets that allows **up to ~100 false flags before vetting by eye**, which is why §5 vetting is mandatory. The v1 rate
  (5/87) is a cross-check.
- Recall cost seen in FA set 2: one clearly real signal (TIC 237290473) is not flagged because of the amplitude-ratio rule
  (spotted secondaries). This matches the 5/22 B26 L1 systems lost to the same rule. It is accepted: the rule stays frozen.
