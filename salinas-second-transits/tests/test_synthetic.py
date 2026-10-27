#!/usr/bin/env python3
"""End-to-end check on fake light curves, with MAST mocked out.

Run from the package folder:  python tests/test_synthetic.py
"""
import shutil
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import st2_pipeline as pl  # noqa: E402

rng = np.random.default_rng(42)
CAD = 10 / 1440
P_TRUE, T1 = 157.9, 2100.0
T2 = T1 + 7 * P_TRUE  # 3205.3, in S70
DEPTH, DUR = 4000e-6, 6 / 24

SECTORS = {  # tic -> {sector: (tstart, tend)}
    1: {30: (2090, 2116), 45: (2530, 2556), 50: (2568, 2594), 70: (3195, 3221)},
    2: {33: (2190, 2216), 64: (3000, 3026)},
    3: {35: (2250, 2276)},
}
TRANSITS = {1: [T1 + k * P_TRUE for k in range(8)], 2: [2200.0], 3: [2260.0]}


def fake_lc(tic, sector):
    a, b = SECTORS[tic][sector]
    t = np.arange(a, b, CAD)
    t = t[~((t > (a + b) / 2 - 0.6) & (t < (a + b) / 2 + 0.6))]  # orbit gap
    f = 1 + 2000e-6 * np.sin(2 * np.pi * t / 3.1) + rng.normal(0, 1000e-6, t.size)
    for tc in TRANSITS[tic]:
        f[np.abs(t - tc) < DUR / 2] -= DEPTH
    return t, f


def main():
    work = Path(tempfile.mkdtemp())
    out = work / "results"
    pd.DataFrame([
        dict(tic=1, type="single", t0=T1, duration_h=6.0, depth_ppm=4000, depth_err=300),
        dict(tic=2, type="single", t0=2200.0, duration_h=6.0, depth_ppm=4000, depth_err=300),
        dict(tic=3, type="single", t0=2260.0, duration_h=6.0, depth_ppm=4000, depth_err=300),
    ]).to_csv(work / "c.csv", index=False)

    pl.tic_info = lambda tic: dict(ra=100.0, dec=-60.0, tmag=11, teff=5800, rstar=1.0, mstar=1.0, contratio=0)
    pl.lc_products = lambda tic: (None, [dict(sector=s, author="TESS-SPOC", exptime=200.0, idx=(tic, s),
                                              tmin=a, tmax=b) for s, (a, b) in SECTORS[tic].items()])
    pl.load_sector_lc = lambda sr, row, cache: fake_lc(*row["idx"])
    # local TOI/CTOI lists so no network is needed
    pd.DataFrame({"TIC ID": [2], "TOI": [9999.01], "TFOPWG Disposition": ["PC"],
                  "Period (days)": [50.0], "Epoch (BJD)": [2459200.0]}).to_csv(work / "exofop_toi.csv", index=False)
    pd.DataFrame({"TIC ID": [5], "CTOI": [5.01], "Period (days)": [3.0]}).to_csv(work / "exofop_ctoi.csv", index=False)

    args = SimpleNamespace(workers=1, snr_min=7.0, depth_lo=0.5, depth_hi=2.0, dur_lo=0.5, dur_hi=2.0,
                           anti_ratio=1.5, p_floor=1.0, force=False, out=str(out))
    cands = pl.load_candidates(work / "c.csv")
    pl.step_coverage(cands, work, 1)
    ev = pl.step_search(cands, work, args)
    pl.step_toi(cands, work, args)
    al = pl.step_aliases(cands, work, args)
    pl.step_report(cands, work, args)

    hits = ev[ev.passes == True]
    print(hits[["tic", "sector", "t2", "depth_ppm", "dur_h", "snr", "anti_snr"]])
    assert set(hits.tic) == {1}, "only TIC 1 should have a second transit"
    h = hits.sort_values("snr").iloc[-1]
    assert abs(h.t2 - T2) < 0.05, h.t2
    assert 0.7 < h.depth_ratio < 1.4, h.depth_ratio
    at = pd.read_csv(next((work / "aliases").glob("1_*.csv")))
    true = at.iloc[(at.period_d - P_TRUE).abs().argmin()]
    print(at.status.value_counts().to_string())
    print("true alias:", true.to_dict())
    assert true.status == "supported", "S50 holds a third transit at the true period"
    assert (at.status == "excluded").sum() > 0, "S45/S50 should exclude some aliases"
    toi = pd.read_csv(work / "toi_match.csv", dtype=str).fillna("")
    assert toi.set_index("tic").loc["2", "toi"] != ""
    assert (out / "shortlist.md").exists() and any((out / "plots").glob("*.png"))
    print("\nALL CHECKS PASSED; outputs in", work)
    return work


if __name__ == "__main__":
    w = main()
    if "--keep" not in sys.argv:
        shutil.rmtree(w)
