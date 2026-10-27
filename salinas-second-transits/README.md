# Second transits for the Salinas et al. candidates

Goal: find a second transit, in TESS Sector 61 or later, for the new candidates in
Salinas et al., "Single-Transit Parameter Estimation with Calibrated Uncertainties"
(arXiv 2610.24882). Each hit becomes a new duotransit whose period aliases can be
resolved the usual way, and a set of them could make a Research Note of the AAS.

Nothing here has touched real TESS data yet. The cloud session that built it could
not reach MAST or ExoFOP, so the pipeline was only checked end to end on fake light
curves (`tests/test_synthetic.py`). Run it on your Mac.

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

## What each step does

1. coverage: TIC lookup (coordinates, R*, M*), every TESS light curve MAST has
   (SPOC, then TESS-SPOC, then QLP per sector) and tess-point pointings. Writes
   `work/coverage.csv`, including `pointed_no_lc_ge61`: sectors 61+ where the star
   fell on silicon but no light curve exists yet (would need TESSCut or TGLC).
2. search: remeasures the paper's event with the same detrending, then slides a box
   of that duration across each new sector (biweight detrend, window 4x duration).
   Noise is the scatter of duration-long bins, so red noise counts. An event passes
   if SNR >= 7, depth and duration are within 0.5 to 2x the reference, and SNR beats
   the strongest *brightening* in that sector by 1.5x (a cheap systematics check).
   Writes `work/events_all.csv`, every top event whether it passes or not.
3. toi: downloads the ExoFOP TOI and CTOI lists and flags any candidate that has
   since become a TOI or CTOI, and whether that TOI's ephemeris passes through the
   paper's event.
4. aliases: for each hit, every P = (T2 - T1)/n down to 1 d, checked against all
   TESS sectors of that star. An alias is `excluded` if a predicted transit falls in
   data and isn't there, `supported` if one is there (a third transit, so worth a
   look first), `allowed` otherwise, and `allowed_below_pmin` if shorter than the
   shortest circular orbit that can give that duration for the TIC star.
   Writes `work/aliases/<cand>.csv` and `work/aliases_summary.csv`.
5. report: `results/shortlist.csv`, `results/shortlist.md` and a plot per hit
   (paper event, new event, full new sector).
6. cnn (optional): refits each new event with the Salinas CNN for depth, duration
   and mid-time with calibrated errors, so both events are measured the same way as
   the paper. Needs their repo and its own environment (Python 3.10, numpy<2, torch):

   ```bash
   git clone https://github.com/helemysm/singletransiters
   python st2_pipeline.py cnn --candidates candidates.csv --cnn-repo singletransiters
   ```

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

## Untested here

MAST, ExoFOP and the CNN step could not be run in the cloud session. The light curve
download, TIC query and ExoFOP column handling follow the libraries' documented
interfaces but haven't seen real responses. If a step errors, it's most likely there.
