# AAVSO and ASAS-SN ground check, 2026-11-01

Command: `python ground_check.py --aavso aavso_kic8462852.csv --asassn` (log: `ground_check_aavso_run.log`).
Jon downloaded the AAVSO file while signed in: 22,578 observations, JD 2459400 to 2461343 (2026 Oct 29).
The raw file stays local and isn't committed. Nothing was posted, submitted or signed in to.

Data: American Association of Variable Star Observers (AAVSO). (2026). AAVSO International Database.

**Result: no trial tc in any window reaches the 5 sigma line for "excluded" or "consistent". AAVSO does
not rule out the 1.1% transit anywhere in window 3, and it shows no dip there either.**

22,234 AAVSO points and 815 ASAS-SN points passed the cuts. ASAS-SN's curve still ends at JD 2460836
(2025 June). IRSA ZTF reset the connection, so ZTF isn't included; the brief lists it as optional.
TESS: TESSCut lists no sector after S82, so none covers window 3.

## Per-window summary (as printed)

    W3 n=+2 (Oct 2026)       excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    W3 n=+1 (Mar 2023)       excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    A776 n=+1 (Oct 2021)     excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    A776 n=+2 (Dec 2023)     excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    A790 n=+1 (Nov 2021)     excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    A790 n=+2 (Jan 2024)     excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    A792 n=+1 (Nov 2021)     excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded  30% consistent   1%
    A792 n=+2 (Jan 2024)     excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    W4 n=+1 (May 2023)       excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    W5 n=+1 (Jun 2023)       excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%

"Untested" covers both no data and data too weak to reach the line. Positive dchi2 means the data disfavour
a transit; negative means they favour one. The line is +/-25.

| Window | tcs with a core point | Highest dchi2 | Lowest dchi2 |
|---|---|---|---|
| W3 n=+2 (Oct 2026) | 225 of 521 (43%) | 13.3 | -3.4 |
| W3 n=+1 (Mar 2023) | 40 of 261 (15%) | 19.4 | 19.4 |
| A776 n=+1 (Oct 2021) | 0 of 47 | | |
| A776 n=+2 (Dec 2023) | 35 of 93 (38%) | -1.1 | -1.3 |
| A790 n=+1 (Nov 2021) | 0 of 21 | | |
| A790 n=+2 (Jan 2024) | 0 of 41 | | |
| A792 n=+1 (Nov 2021) | 119 of 137 (87%) | 8.2 | -1.2 |
| A792 n=+2 (Jan 2024) | 0 of 273 | | |
| W4 n=+1 (May 2023) | 1 of 2 | 4.4 | 4.4 |
| W5 n=+1 (Jun 2023) | 80 of 101 (79%) | 4.2 | -2.7 |

The closest to a verdict is W3 n=+1 at dchi2 19.4, or about 4.4 sigma, on 15% of its tcs.

The A792 n=+1 per-point figures (30% excluded, 1% consistent) are not a result. That statistic treats
points from one night as independent, which the main statistic exists to correct.

The groups sum as if independent, but one observer's bands on one night share weather and comparison
stars. Using Johnson V alone, W3 n=+2 tops out at dchi2 10.7 and still has no verdict.

## Who covers window 3 (2026 Oct 18 to 23 UTC)

19 single visits on 6 UTC evenings or early mornings, from 3 observers. There are no time series, so
the hours of continuous coverage round to zero. The residual is against the group's median over +/-15 d;
positive is fainter. The transit would be about +12 mmag.

| UTC | Observer | B (mmag) | V (mmag) | R (mmag) |
|---|---|---|---|---|
| Oct 19 02:18 | SJMJ (DSLR) | +17 | -2 | -12 |
| Oct 19 19:50 | DFS | | +17 | |
| Oct 21 02:10 | SJMJ (DSLR) | -35 | +1 | +36 |
| Oct 21 19:48 | DFS | +14 | -4 | +6 |
| Oct 21 19:57 | DUBF | +4.5 | -12 | -12 |
| Oct 22 19:48 | DFS | -4 | -16 | -10 |
| Oct 23 19:44 | DFS | +12 | -18 | +27 |

Night-to-night scatter in the baselines is 7 to 24 mmag depending on the group, comparable to or larger
than the transit depth, so no single visit decides anything. The one visit worth noting is DFS V on
Oct 19 at +17 mmag against 7.4 mmag scatter. It's a single point, with no other observer that night.

## Changes

- `load_aavso` reads AAVSO's newer export (`jd`, `mag`, `uncertainty`, `band`, `observer`, `fainterthan`)
  as well as the old column names. `test_ground_check.py` still passes, but it only exercises the old names.
- `main` forces the time and magnitude columns to floats. A manual `--asassn --tess --ztf` run on 2026-11-01
  (not the scheduled task) crashed in A776 n=+2 because ZTF's empty result made the error column objects
  (`ground_check_2026-11-01.log`). That run changed no results. `main` now also drops empty sources (c6f9d01).
- `ground_photometry_merged.csv`, `epoch_coverage.csv` and the `win_*.png` plots now include AAVSO.
- The scheduled recheck needed no change: it finds `aavso_kic8462852.csv` in this folder and passes `--aavso`.
