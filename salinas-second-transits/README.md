# Second transits for the Salinas et al. candidates

Goal: find a second transit, in TESS Sector 61 or later, for the new candidates in
Salinas et al., "Single-Transit Parameter Estimation with Calibrated Uncertainties"
(arXiv 2610.24882). Each hit becomes a new duotransit whose period aliases can be
resolved the usual way, and a set of them could make a Research Note of the AAS.

Round 1 ran on Windows (commit 73ff358); round 2 (FOLLOWUP_BRIEF_2026-10-27) ran on
macOS and reproduced round 1 before the round 2 changes. `tests/test_synthetic.py`
checks every step end to end on fake light curves.

## What's here

| File | What it is |
|---|---|
| `candidates.csv` | Tables 4 to 6 of the paper, parsed: 173 rows, 169 stars (108 single, 57 multi-transit with a period, 4 multi-planet systems with 2 rows each). `t0` is BTJD (the paper's column says BJD but the values are BJD minus 2457000). |
| `salinas_2610.24882_fulltext.txt` | Paper text the tables were parsed from. |
| `parse_salinas_tables.py` | Rebuilds `candidates.csv` from that text. |
| `st2_pipeline.py` | The pipeline (steps below). |
| `st2core.py` | Box search, alias and period-limit maths, no network. |
| `tests/test_synthetic.py` | Injects transits into fake sectors, runs every step with MAST mocked, checks it finds them. |
| `run_all.sh` | Makes a venv, installs, runs the test, then the full pipeline. |
| `period_hypotheses.csv` | Periods to tabulate epoch by epoch in the alias step. |
| `cnn_extra.csv` | Extra transits for the CNN step. |
| `tesscut_phot.py` | TESScut cutouts: aperture photometry, difference images, centroids. |
| `round2/` | One-off round 2 analyses (sections 4, 5, 6, 8 of the brief); outputs in `results/round2/`. |

## Run it

```bash
cd salinas-second-transits
./run_all.sh                 # most of the time goes to MAST downloads
open results/shortlist.md
```

Or step by step, which is rerunnable (finished candidates are skipped unless `--force`):

```bash
python st2_pipeline.py coverage --candidates candidates.csv
python st2_pipeline.py search   --candidates candidates.csv
python st2_pipeline.py toi      --candidates candidates.csv
python st2_pipeline.py aliases  --candidates candidates.csv
python st2_pipeline.py report   --candidates candidates.csv
```

By default only `single` and `multiplanet` rows are searched; the multi-transit ones
already have a BLS period. Use `--types all` to include them. If MAST throttles you,
use `--workers 1`.

On macOS make the venv with Python 3.12 (`uv venv --python 3.12 .venv`, then
`uv pip install -r requirements.txt`); `run_all.sh` uses an existing `.venv`. numba
only has its non-thread-safe workqueue layer there, so wotan calls are serialised
behind a lock; downloads still run in parallel.

## What each step does

1. coverage: TIC lookup (coordinates, R*, M*), every TESS light curve MAST has
   (SPOC, then TESS-SPOC, then QLP per sector) and tess-point pointings. Writes
   `work/coverage.csv`, including `pointed_no_lc_ge61`: sectors 61+ where the star
   fell on silicon but no light curve exists yet (would need TESSCut or TGLC).
   "Sector 1751" is the 7.2 d 3I/ATLAS special pointing (2026 Jan 15-22): real data,
   listed under `special_sectors`, used by the alias step but not searched as new.
2. search: remeasures the paper's event with the same detrending, then slides a box
   of that duration across each new sector (biweight detrend, window 4x duration).
   Noise is the scatter of duration-long bins, so red noise counts. An event passes
   if SNR >= 7, depth and duration are within 0.5 to 2x the reference, and SNR beats
   the strongest *brightening* in that sector by 1.5x (a cheap systematics check).
   Writes `work/events_all.csv`, every top event whether it passes or not.
3. toi: downloads the ExoFOP TOI and CTOI lists and flags any candidate that has
   since become a TOI or CTOI, and whether that TOI's ephemeris passes through the
   paper's event.
4. aliases: T2 is the highest-SNR *clean* passing event (at least 0.5 d from a gap
   and 80% covered; top SNR only if none is clean). Every P = (T2 - T1)/n down to
   1 d is checked against all TESS sectors of that star. Only clean predicted epochs
   count: an alias is `excluded` if one has a box depth below
   max(depth0 - 3 sigma, 0.5 depth0), `supported` if one matches depth0 (a third
   transit), `allowed` otherwise, and `allowed_below_pmin` if shorter than the
   shortest circular orbit that can give that duration for the TIC star. Partly
   covered or near-edge epochs are `inconclusive`. Coverage counts each point's own
   cadence, so 30-min early sectors are used. Also: a common-period test for
   candidates with 2+ distinct events (`work/common_period.csv`, flags `periodic`
   and `chance`), per-epoch tables `work/aliases/<cand>_epochs.csv` for the period
   hypotheses and surviving aliases, `work/aliases/hypotheses_summary.csv`, and plots
   of each excluding epoch in `results/epoch_plots/`.
   Writes `work/aliases/<cand>.csv` and `work/aliases_summary.csv`.
5. report: `results/shortlist.csv`, `results/shortlist.md` and a plot per hit
   (paper event, new event, full new sector).
6. cnn (optional): refits each new event with the Salinas CNN for depth, duration
   and mid-time with calibrated errors, so both events are measured the same way as
   the paper. Needs their repo and its own environment (Python 3.10, numpy<2, torch):

   ```bash
   git clone https://github.com/helemysm/singletransiters
   uv venv --python 3.10 .venv-cnn && uv pip install --python .venv-cnn/bin/python -r singletransiters/requirements.txt
   .venv-cnn/bin/python st2_pipeline.py cnn --candidates candidates.csv --cnn-repo singletransiters
   ```

   It fits every passing event, the paper event of each such candidate, and anything
   in `cnn_extra.csv`. QLP files from S94 on carry a read-only column that broke
   their `preprocess_lc` (`LightCurve.bin()`: "output array is read-only", 16 failures
   in round 1); the wrapper now rebuilds the light curve from copied arrays first.

## Reading the results

- Box depths and durations are screening numbers. They run a little shallow and
  short of the truth, but the paper event and the new event are measured the same
  way, so their ratio is what matters. Use the CNN step or your own fit for numbers
  that go in a note.
- The paper already mentions later data for three systems, and the report flags them:
  TIC 302377820 (a transit in S68 matching candidate b), TIC 404736082 (more
  transits in S70) and TIC 318706605 (extra dips in S42).
- `edge_dist_d` under 0.5 d means the dip sits near a gap or sector edge. Check those
  against momentum dumps and scattered light before trusting them.
- Before calling anything a duotransit, vet it as you normally would: centroids or a
  difference image (TESSCut or TGLC), nearby stars in Gaia, odd/even and secondary
  eclipse checks across the two events, and the paper's own notes on blends.

## Known limits

- The alias exclusion rule treats an epoch more than 3 sigma shallower than the
  reference as "no transit". For high-SNR signals measured with different pipelines
  (TESS-SPOC reference, QLP epochs) that excludes true periods whose epochs hold
  transits at 60-80% of the reference depth. See the round 2 summary before
  trusting an `excluded` status for a strong signal.
- TESScut times are barycentre-corrected for another direction than the target (up
  to ~25 s off); `round2/s6_qlp_times.py` shows how to correct them with the JPL
  Horizons TESS ephemeris.
