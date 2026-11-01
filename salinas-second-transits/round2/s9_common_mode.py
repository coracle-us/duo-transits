#!/usr/bin/env python3
"""Sector-wide dips: find them, measure them, and re-search the sectors they affect.

Top events cluster in time across unrelated stars in a few sectors. The trigger
was a QLP S94 dip at BTJD 3874.8 that passed the search in six stars.

  python round2/s9_common_mode.py scan      # run on light curves saved BEFORE the mask
  python round2/s9_common_mode.py tesscut   # S94 in our own TESScut photometry
  python round2/s9_common_mode.py research  # re-search the masked sectors only, then rebuild events_all.csv
                                            # with: python st2_pipeline.py search --candidates candidates.csv

The re-search keeps every other sector's events as they were, so products MAST has
added since can't change the comparison.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import st2_pipeline as pl  # noqa: E402
import st2core as core  # noqa: E402

WORK, OUT = HERE / "work", HERE / "results" / "round2"
BIN_D = 0.5
FLAG_MIN_STARS, FLAG_MIN_FRAC = 8, 0.7


def sector_files():
    ev = pd.read_csv(WORK / "events_all.csv").dropna(subset=["author"])
    author = ev.drop_duplicates(["tic", "sector"]).set_index(["tic", "sector"]).author.to_dict()
    by = {}
    for p in sorted((WORK / "lc").glob("*_s*.npz")):
        m = re.fullmatch(r"(\d+)_s(\d+)\.npz", p.name)
        if m and (int(m[1]), int(m[2])) in author:
            tic, s = int(m[1]), int(m[2])
            by.setdefault(s, []).append((tic, author[(tic, s)], p))
    return by


def scan():
    """Bins where at least half the FFI stars in a sector dip by more than 3 sigma."""
    rows = []
    for s, stars in sorted(sector_files().items()):
        stars = [x for x in stars if x[1] != "SPOC"]
        if len(stars) < 4 or not pl.is_new_sector(s):
            continue
        lcs = [(d["t"], d["f"]) for d in (np.load(x[2]) for x in stars)]
        lo, hi = min(t.min() for t, _ in lcs), max(t.max() for t, _ in lcs)
        for b in np.arange(np.floor(lo), hi, BIN_D):
            sig, cad = [], []
            for t, f in lcs:
                w = (t >= b) & (t < b + BIN_D)
                cad.append(int(w.sum()))
                if w.sum() >= 8:
                    sig.append((1 - np.median(f[w])) / (core.robust_std(f) / np.sqrt(w.sum())))
            if len(sig) < 4 or np.mean(np.array(sig) > 3) < 0.5:
                continue
            frac = float(np.mean(np.array(sig) > 3))
            rows.append(dict(sector=s, n_ffi_stars=len(stars), bin_start_btjd=b, n_stars_with_data=len(sig),
                             frac_dip_gt_3sigma=frac, median_sigma=float(np.median(sig)),
                             median_cadences=int(np.median(cad)),
                             flagged=len(sig) >= FLAG_MIN_STARS and frac >= FLAG_MIN_FRAC))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "common_mode_scan.csv", index=False)
    print(df.round(2).to_string(index=False))

    # every star in the flagged windows, SPOC included
    rows = []
    files = sector_files()
    for s, windows in core.COMMON_MODE_WINDOWS.items():
        for lo, hi in windows:
            for tic, author, p in files[s]:
                d = np.load(p)
                t, f = d["t"], d["f"]
                w = (t >= lo) & (t <= hi)
                o = ((t > lo - 3) & (t < lo)) | ((t > hi) & (t < hi + 3))
                rows.append(dict(sector=s, window_start=lo, window_end=hi, tic=tic, author=author, n_window=int(w.sum()),
                                 dip_ppm=(1 - np.median(f[w])) * 1e6 if w.sum() >= 8 else np.nan,
                                 err_ppm=core.robust_std(f[o]) / np.sqrt(max(w.sum(), 1)) * 1e6))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "common_mode_windows_by_star.csv", index=False)
    print(df.round(0).to_string(index=False))


def tesscut():
    """S94 in our own aperture photometry from the same FFIs, per star and over time."""
    import tesscut_phot as tp
    cov = pd.read_csv(WORK / "coverage.csv").set_index("tic")
    rows, prof = [], {}
    bins = np.arange(3872.0, 3878.0, 0.25)
    for tic, author, _ in sector_files()[94]:
        ra, dec = cov.ra[tic], cov.dec[tic]
        t, f, _, _ = tp.sector_lightcurve(tic, ra, dec, 94, WORK / "tesscut")
        t, f = np.asarray(t), np.asarray(f)
        o = ((t > 3871) & (t < 3873.5)) | ((t > 3876.6) & (t < 3879))
        r = f / np.polyval(np.polyfit(t[o], f[o], 1), t)
        name = next((WORK / "tesscut").glob(f"**/tess-s0094-*_{ra:.6f}_{dec:.6f}_15x15_astrocut.fits")).name
        cam, ccd = re.match(r"tess-s0094-(\d)-(\d)", name).groups()
        w = (t > 3874.5) & (t < 3875.1)
        rows.append(dict(tic=tic, lc_author=author, camera=int(cam), ccd=int(ccd),
                         tesscut_dip_ppm=(1 - np.median(r[w])) * 1e6,
                         scatter_ppm=core.robust_std(r[o]) * 1e6))
        prof[tic] = [np.median(r[(t >= b) & (t < b + 0.25)]) - 1 for b in bins]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s94_tesscut_by_star.csv", index=False)
    print(df.round(0).to_string(index=False))
    med = pd.DataFrame(prof, index=bins).median(axis=1) * 1e6
    med.rename("median_residual_ppm").rename_axis("bin_start_btjd").to_csv(OUT / "s94_tesscut_median_profile.csv")
    print(med.round(0).to_string())


def research():
    """Replace the events of the masked sectors; every other sector's events stay as they are."""
    (WORK / "lc").mkdir(exist_ok=True)
    done = 0
    for p in sorted((WORK / "per_cand").glob("*.json")):
        r = json.loads(p.read_text())
        secs = {e["sector"] for e in r["events"]} & set(core.COMMON_MODE_WINDOWS)
        if not secs:
            continue
        sr, prods = pl.lc_products(r["tic"])
        best = pl.best_per_sector(prods)
        keep = [e for e in r["events"] if e["sector"] not in secs]
        for s in sorted(secs):
            old = [e for e in r["events"] if e["sector"] == s]
            if best[s]["author"] == "SPOC":
                keep += old
                continue
            new = pl.search_sector(r["tic"], s, best[s], sr, WORK / "lc_cache", WORK,
                                   r["dur_ref_d"], r["depth_ref"])
            keep += new
            print(f"{r['cand_id']} S{s} {best[s]['author']}: "
                  f"{[round(e['tc'], 2) for e in old if 'tc' in e]} -> {[round(e['tc'], 2) for e in new if 'tc' in e]}")
            done += 1
        r["events"] = sorted(keep, key=lambda e: e["sector"])
        p.write_text(json.dumps(r, default=float, indent=1))
    print(f"re-searched {done} candidate-sectors")


if __name__ == "__main__":
    {"scan": scan, "tesscut": tesscut, "research": research}[sys.argv[1]]()
