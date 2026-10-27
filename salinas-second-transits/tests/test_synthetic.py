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


def check_alias_event_choice():
    """A clean lower-SNR event beats a high-SNR edge dip; with no clean event, top SNR wins."""
    ev = pd.DataFrame([
        dict(cand_id="a", tic=1, sector=95, t2=3899.75, snr=12.0, edge_dist_d=0.15, frac=0.54),
        dict(cand_id="a", tic=1, sector=68, t2=3171.99, snr=8.5, edge_dist_d=3.0, frac=1.0),
        dict(cand_id="b", tic=2, sector=70, t2=3200.0, snr=9.0, edge_dist_d=0.2, frac=0.9),
        dict(cand_id="b", tic=2, sector=71, t2=3230.0, snr=7.5, edge_dist_d=2.0, frac=0.6),
    ])
    h = pl.alias_events(ev).set_index("cand_id")
    assert h.loc["a", "sector"] == 68 and "not clean" in h.loc["a", "alias_event_reason"]
    assert h.loc["b", "sector"] == 70 and h.loc["b", "alias_event_reason"].startswith("no clean event")


def check_edge_epochs_inconclusive():
    """A flat epoch 0.2 d from a gap is inconclusive; the same epoch mid-data excludes."""
    t = np.concatenate([np.arange(0, 10, CAD), np.arange(12, 22, CAD)])
    f = 1 + rng.normal(0, 300e-6, t.size)
    v = pl.core.epoch_verdicts([9.8, 5.0, 11.0], t, f, DUR, DEPTH, 300e-6)
    assert list(v["verdict"]) == ["inconclusive", "excludes", "no data"], v["verdict"]
    assert v["verdict_any_coverage"][0] == "excludes"


def check_common_period():
    """Four later dips on a 23.75 d ephemeris through T1 give that period; a stray dip
    lets a short chance period win overall, but not among P >= 5 d."""
    times = [2396.13 + 23.7535 * k + d for k, d in ((30, 0.02), (31, -0.03), (33, 0.01), (40, 0.0))]
    cp = pl.core.common_period(2396.13, times, 4 / 24)
    assert cp["k"] == 5 and cp["N"] == 5 and abs(cp["P"] - 23.7535) < 0.01, cp
    lp = pl.core.common_period(2396.13, times + [3000.4], 4 / 24, p_floor=5.0)
    assert lp["k"] == 5 and lp["N"] == 6 and abs(lp["P"] - 23.7535) < 0.01, lp


def main():
    check_alias_event_choice()
    check_common_period()
    check_edge_epochs_inconclusive()
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
    pd.DataFrame({"tic": [1], "period_d": [157.9]}).to_csv(work / "hyp.csv", index=False)

    args = SimpleNamespace(workers=1, snr_min=7.0, depth_lo=0.5, depth_hi=2.0, dur_lo=0.5, dur_hi=2.0,
                           anti_ratio=1.5, p_floor=1.0, force=False, out=str(out),
                           hypotheses=str(work / "hyp.csv"))
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
    at = pd.read_csv(next(p for p in (work / "aliases").glob("1_*.csv") if not p.name.endswith("_epochs.csv")))
    true = at.iloc[(at.period_d - P_TRUE).abs().argmin()]
    print(at.status.value_counts().to_string())
    print("true alias:", true.to_dict())
    assert true.status == "supported", "S50 holds a third transit at the true period"
    assert (at.status == "excluded").sum() > 0, "S45/S50 should exclude some aliases"
    hy = pd.read_csv(work / "aliases" / "hypotheses_summary.csv").iloc[0]
    print("hypothesis:", hy.to_dict())
    assert abs(hy.period_used_d - P_TRUE) < 0.01 and hy.status_now == "supported"
    ep = pd.read_csv(next((work / "aliases").glob("1_*_epochs.csv")))
    assert {"paper event", "alias event (T2)", "supports"} <= set(ep.verdict)
    toi = pd.read_csv(work / "toi_match.csv", dtype=str).fillna("")
    assert toi.set_index("tic").loc["2", "toi"] != ""
    assert (out / "shortlist.md").exists() and any((out / "plots").glob("*.png"))
    print("\nALL CHECKS PASSED; outputs in", work)
    return work


if __name__ == "__main__":
    w = main()
    if "--keep" not in sys.argv:
        shutil.rmtree(w)
