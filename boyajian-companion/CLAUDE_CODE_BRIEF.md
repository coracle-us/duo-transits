# Brief for Claude Code on Jon's machine: Boyajian's star ground photometry check

Goal: find out whether public ground photometry shows (or rules out) the companion's transit in
window 3 (2026 Oct 18 to 24 UTC) and at the past epochs of the periods Kepler still allows. The
cloud session can't reach AAVSO, ASAS-SN, IRSA or MAST's API, so this runs locally.

Everything here is read-only. Don't post, submit or email anything. Don't create accounts or sign in
anywhere on Jon's behalf.

## Setup
```
cd duo-transits            # or wherever Jon keeps this; make a folder boyajian-companion/
python -m pip install -r requirements.txt
```
Copy in `ground_check.py`, `kepler_alias_scan.py`, `requirements.txt` from this folder.

## Steps
1. Kepler check (sanity, should match the cloud run):
   `python kepler_alias_scan.py 770 845`
   Expected surviving ranges: 776.17-776.60, 790.50-790.66, 792.08-792.92, 793.25-793.40 d.
2. ASAS-SN and TESS coverage (no login needed):
   `python ground_check.py --asassn --tess`
   Note which TESS sectors TESSCut lists. If any sector after 82 covers 2026 Oct 18 to 24,
   say so: that would be the best possible test of window 3 (pull its FFI cutout and look).
3. AAVSO (Jon does this part himself): sign in to AAVSO, open Data Download, star "KIC 8462852",
   JD 2459400 to today, all bands, CSV. Save as `aavso_kic8462852.csv` here. Then:
   `python ground_check.py --aavso aavso_kic8462852.csv --asassn --ztf`
   (ZTF is optional: the star is near ZTF saturation, so it may add little.)
4. Report back, pasting into the astrophysics-research channel:
   - the printed per-window summary (excluded / consistent / untested fractions),
   - which groups (observers, bands) cover window 3, nights and hours,
   - any tc flagged "consistent" (a dip). If one shows up, attach its `win_*.png` and check it isn't a
     single observer's bad night or a known long-term dimming episode.
   - the TESS sector list.

## What the numbers mean
- `excluded`: data inside the transit core sit at the normal level, so a 1.1% transit at that
  mid-time is ruled out at >5 sigma.
- `consistent`: data prefer a 1.1% dip there. Not a detection by itself (dust dips happen), but
  it's the thing to look at.
- `untested`: no usable points in the core.
Window 3 is 5.1 days wide and the transit lasts ~21 h, so single-site evening coverage leaves gaps.
A clean result is "X% of the window 3 period range is excluded", plus whatever remains.
