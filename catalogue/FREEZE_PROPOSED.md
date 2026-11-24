# The freeze (step 2 of 2): what is proposed

**Status: PROPOSED, not approved. No candidate's light curve has been read.** On approval the files
below are hashed into `FREEZE.json`, the approval goes in `APPROVALS.md`, and the run may start.

## Files that would be frozen

`prereg_catalogue.yaml`, `defs.csv`, `defs_not_scored.csv`, `tic_info.csv`, `check_periods.csv`,
`validation/v2_truth.csv`, `validation/v3_truth.csv`, `validation/arxiv_2611.03019/published_periods.csv`,
`validation/arxiv_2611.03019/source/anc/table2_full.csv`, `build_defs.py`, `build_check_periods.py`,
`validation/build_truth.py`, `catrun.py`, `run_catalogue.py`, `sheets.py`, `make_sheets.py`,
`VETTING.md`, `merge_verdicts.py`, `report.py`, and every file of `calibration/mono/round3/FREEZE3.json`.

## Checks done before the freeze

`driver_check/DRIVER_CHECK.md`: the driver reproduces round 3's rows on 116 of 116 stars, and its
duotransit table agrees with a blind re-implementation on 24,966 of 24,966 alias statuses.

## Choices made while building, where the pre-registration left room

Each is fixed in the code above. None was made after seeing a candidate's result.

1. **"multi".** More than one scored row on the star, or more than one TOI, or more than one CTOI
   (ExoFOP tables of 2026-11-22). A star with one TOI and one CTOI is not multi on that ground: they
   are usually the same signal. 123 of 1,183 rows are multi. (Counting TOIs and CTOIs together would
   make 187.)
2. **Excluded fraction of ln P.** Two numbers per row: over 1 to 1000 d (the quantity the calibration
   measured), and over 1 d to the longest period the data can test (the pre-registration's wording).
3. **Stated periods.** A period stated by a catalogue, a paper or this project is tested exactly as the
   calibration tested a true period (any epoch that excludes), after the run, never as an input.
4. **V2, "not an alias".** A solved period that is not a whole fraction of the gap between the two
   catalogued times is counted separately, not as an exclusion by the table.
5. **Depth error of a Hawthorn second transit.** That of the pipeline's event within half a duration of
   Tc2, passing or not; else a refit at Tc2; else the reference's.
6. **A third event on a Hawthorn row.** Two things are kept: for each alias of the two-transit table,
   whether it puts an epoch on the third event within that epoch's T1 window; and a second table
   through Tc1 and the third event.
7. **Events at another catalogued transit on the same star.** Marked on the sheet and in the table.
   Not a reason to reject: acceptance does not claim the same body.
8. **Sheets are anonymous.** File names are hashes; no star name, period or alias table is shown.
9. **Neighbour panel.** Up to five of the 10 nearest stars within 0.5 degree and 1 mag in Tmag that have
   a light curve of the same product type in the sector, on the same camera and CCD.
10. **Difference images.** Only from cut-outs already on disk; nothing is downloaded for them. Most
    sheets will say "no difference image".
11. **Vetter check (V1).** False events are taken from arm A1 and arm B, true events from arm A2 and
    arm B. Arm-B sheets show the injected planet in every product of the target.
12. **Catalogued depths out of range.** 14 CTOI rows have a catalogued depth above 25% or below
    50 ppm. They are run as they are and flagged.
13. **Network failures.** A star whose light curves fail to load is retried and, if it still fails,
    reported as an error, not as "no light curve".

## Known before the run

- V3 has 3 published periods on our stars (TOI-2433.01, TOI-6981.01, and TOI-6669.01 as Hawthorn H72),
  so its criterion is "none excluded". 138 of our rows are in that paper's table.
- K2 rows (164) and this project's own periods (21 rows) are not scored.
