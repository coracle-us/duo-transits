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
