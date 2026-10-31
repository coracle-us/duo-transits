# ASAS-SN ground check (2026-10-31, scheduled recheck)

Command: `.venv/bin/python ground_check.py --asassn --tess` (log: `ground_check_asassn_run.log`).
No AAVSO file was found, so this is ASAS-SN only. Nothing was posted, submitted or signed in to.

**Result: every trial tc in all 10 windows is "untested". This run excludes nothing and finds no dip.
Window 3 n=+2 (Oct 2026) has no data at all.**

## Service state
`asassn_status.py` exits 0 and lists 6 block servers, and the download worked. The server banner
still shows the 13 Oct 2026 maintenance message. The light curve stops at 2025-06-10. I can't tell
whether that is a side effect of the maintenance or the current end of the data for this star.

## Per-window summary (as printed)

    W3 n=+2 (Oct 2026)       excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    W3 n=+1 (Mar 2023)       excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    A776 n=+1 (Oct 2021)     excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    A776 n=+2 (Dec 2023)     excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    A790 n=+1 (Nov 2021)     excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    A790 n=+2 (Jan 2024)     excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    A792 n=+1 (Nov 2021)     excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    A792 n=+2 (Jan 2024)     excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    W4 n=+1 (May 2023)       excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%
    W5 n=+1 (Jun 2023)       excluded    0%  consistent    0%  untested  100%   | per-point statistic: excluded   0% consistent   0%

## ASAS-SN points
ASAS-SN id 60130747136: 939 points loaded, 815 kept (quality G, mag_err < 0.05). The 124 rejected
points are all quality B or have large errors, and none falls in a window's core range.

| Group | Kept points | Nights | First | Last |
|---|---|---|---|---|
| ASASSN:V:bd | 211 | 210 | 2015-02-24 | 2018-11-13 |
| ASASSN:g:bt | 283 | 277 | 2018-04-15 | 2025-06-10 |
| ASASSN:g:bD | 321 | 313 | 2019-03-28 | 2024-10-22 |

A separate one-off query of the unfiltered curve confirmed that the June 2025 end comes from the
server, not our filter: 0 raw points after JD 2460837.

## W3 n=+2 (Oct 2026), JD 2461326-2461338
**No group has points: 0 nights, 0 hours.** The last point is about 16 months before the window.
TESSCut lists sectors 14, 15, 41, 54, 55, 75, 81 and 82, so TESS doesn't cover it either. This
window remains untested.

## Why the other windows are untested
Trial tcs total 1497, and only 35 of them were actually compared with data. All 35 are in A776 n=+2.

| Window | Trial tcs | With a core point | Outcome |
|---|---|---|---|
| W3 n=+1 | 261 | 12 | one g:bD point; skipped, only 5 baseline points (needs 10) |
| A776 n=+2 | 93 | 35 | one g:bt point, used; dchi2 between -1.27 and 0 |
| A790 n=+1 | 21 | 21 | one g:bD point; skipped, only 7 baseline points (needs 10) |
| A776 n=+1, A790 n=+2, A792 n=+1, A792 n=+2, W4, W5 | 601 | 0 | no point within 0.20 d of any tc |
| W3 n=+2 | 521 | 0 | no data (see above) |

- **Sampling is too sparse.** There is about one kept point every 4.6 days, and the transit core
  is only 0.4 d wide.
- **More coverage of this kind wouldn't settle it.** With one night in the core, a single group
  needs a nightly scatter below 0.0024 mag to reach the |dchi2| > 25 threshold for a 0.012 mag
  transit. The whole-curve scatter (MAD) of the g groups is 0.008 and 0.014 mag, which gives
  dchi2 of about 2 at best. As configured, sparse ASAS-SN points like these can't exclude a tc.

## "Consistent" tcs
None. The most negative dchi2 anywhere is -1.27, far from the -25 threshold.

The A776 n=+2 point (JD 2460282.611, g:bt) plots about 0.03 mag faint in `win_A776_n2.png`. It is
not a detection. The g:bt points on JD 2460279.6, 2460280.6 and 2460283.6 are similarly faint, so
the offset lasts about 4 days, and other g:bt nights nearby (for example JD 2460271.65) show
similar scatter. It comes from one group only.

## Plots
`win_W3_n2.png`, `win_W4_n1.png` and `win_W5_n1.png` have no points, because nothing falls within
3 d of those windows. The other seven have 1 to 5 points each.

## What was and wasn't tested
- **Tested:** the ASAS-SN download and id check end to end (the 1" cone returns exactly one id,
  60130747136), the TESSCut sector list, and 35 trial tcs in A776 n=+2, with no useful constraint.
- **Not tested:** window 3 n=+2 and, in practice, every other window.
- **AAVSO:** still the only source named so far that could test window 3. No CSV was available
  for this run.

## Scheduled task
`boyajian-asassn-recheck` could **not** be disabled from this run: the update needs a person's
approval. It stays enabled until someone turns it off, and it will rerun this check on its next
firing.
