"""
Quiet no-transit training windows (RULES.md amendment B1): up to 2 random out-of-transit windows per CP/KP and EB
star, seed 20261125, from the cached binned light curves that make_windows.py used.

Input:  data/features.csv (for the star list and their light curves), data/sample.csv, data/eb_sample.csv,
        data/toi_2026-11-24.csv, the binned-light-curve cache
Output: data/features_quiet.csv, data/windows_flux_quiet.npz
Usage:  .venv/bin/python -I make_quiet.py
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir; our modules live here

import numpy as np
import pandas as pd

import lcwin
from make_windows import D, HERE, exclusion_times, load

SEED_B = 20261125


def main():
    F = pd.read_csv(os.path.join(D, "features.csv"), usecols=["obj_id", "tic", "kind", "sector", "author"])
    s = pd.read_csv(os.path.join(D, "sample.csv"))
    s = s[s.drop_reason.isna()]
    eb = pd.read_csv(os.path.join(D, "eb_sample.csv"))
    eb = eb[eb.drop_reason.isna()]
    prsa = pd.read_csv(os.path.join(HERE, "..", "..", "sidequests", "eb_timing", "data", "prsa2022_tess_ebs.csv"))
    eb = eb.assign(phis=eb.tic.map(prsa.drop_duplicates("TIC").set_index("TIC")["Phis-pf"]))
    toi = pd.read_csv(os.path.join(D, "toi_2026-11-24.csv"), usecols=["TIC ID", "Epoch (BJD)", "Period (days)"])
    toi = toi.rename(columns={"TIC ID": "tic", "Epoch (BJD)": "epoch", "Period (days)": "period"})

    # light curves "used above" for each training star: those behind its pos / eb windows (valid or not)
    w = pd.read_csv(os.path.join(D, "windows.csv"))
    w = w[w.toi.isin(s[s.disposition.isin(["CP", "KP"])].toi)]
    ew = pd.read_csv(os.path.join(D, "eb_windows.csv"))
    ew = ew[ew.tic.isin(eb.tic)]
    lcs = pd.concat([w[["tic", "sector", "author"]], ew[["tic", "sector", "author"]]]).drop_duplicates()
    # T14 per star: the longest TOI duration among its CP/KP TOIs, else the EB duration
    dur = s[s.disposition.isin(["CP", "KP"])].groupby("tic").duration_h.max()
    dur = pd.concat([dur, eb.set_index("tic").duration_h[~eb.tic.isin(dur.index).values]])
    stars = sorted(set(F[F.kind.isin(["pos", "eb"])].tic))
    print(f"training stars: {len(stars)}")

    rng = np.random.default_rng(SEED_B)
    rows, fluxes, short = [], [], 0
    for tic in stars:
        g = lcs[lcs.tic == tic]
        data = {}
        for r in g.itertuples():
            t, f = load(r.tic, r.sector, r.author)
            if t is not None and len(t) >= lcwin.L:
                data[(int(r.sector), r.author)] = (t, f)
        keys = sorted(data)
        kept, draws = [], 0
        while len(kept) < 2 and draws < 200 and keys:
            draws += 1
            if kept and draws <= 100 + kept[0][3] and len(keys) > 1:
                pool = [k for k in keys if k[0] != kept[0][0][0]]  # different sector first (B1)
            else:
                pool = keys
            key = pool[rng.integers(len(pool))]
            t, f = data[key]
            T = rng.uniform(t.min() + 2, t.max() - 2)
            ex = exclusion_times(tic, toi, eb, t.min(), t.max())
            if len(ex) and np.min(np.abs(ex - T)) < 3:
                continue
            if kept and key == kept[0][0] and abs(T - kept[0][1]) * 24 < lcwin.WIN_H:
                continue
            if not lcwin.valid(t, T, dur.get(tic, np.nan)):
                continue
            x, filled = lcwin.window(t, f, T)
            ft = None if x is None else lcwin.features(x, filled)
            if ft is None:
                continue
            kept.append((key, T, ft, draws))
            rows.append(dict(obj_id=f"QUIET{tic}_{len(kept)}", tic=tic, kind="quiet", label=0, disposition="quiet",
                             sector=key[0], author=key[1], t_c=T, **ft))
            fluxes.append(np.where(filled, 0.0, x).astype(np.float32))
        short += len(kept) < 2
    Q = pd.DataFrame(rows)
    Q.to_csv(os.path.join(D, "features_quiet.csv"), index=False)
    np.savez_compressed(os.path.join(D, "windows_flux_quiet.npz"), x=np.stack(fluxes), t_rel_h=lcwin.T_REL_H)
    print(f"quiet windows: {len(Q)} on {Q.tic.nunique()} stars; stars with fewer than 2: {short}")
    print(f"box SNR of quiet windows: median {Q.snr.median():.1f} (10-90% {Q.snr.quantile(.1):.1f}-"
          f"{Q.snr.quantile(.9):.1f})")


if __name__ == "__main__":
    main()
