"""
Period constraint from K2 alone for the K2-observed single TOIs with no K2 event (TOI-5523.01,
TOI-5696.03): the fraction of log P (pilot.P_GRID) excluded because a predicted transit
T = epoch - k P would fall in a flat K2 window (k2lib v3 window test, K2PDC; repeated with K2SFF
as a second reduction), and the shortest period NOT excluded.
TOI-5560.01 is included as a control (its K2 transits must leave P = 140.494 and 70.247 d allowed).
Output: toi_singles_k2/derived/k2_period_limits.csv
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "toi_singles_k2"
sys.path.insert(0, str(OUT / "scripts"))
import s03_k2_search as s3  # noqa: E402

pilot = s3.k2lib.pilot


def status(g, epoch, D, Ps):
    Tg = g.T_pred.to_numpy()
    flat = (g.verdict == "excluded").to_numpy()
    out = []
    for P in Ps:
        ks = np.arange(np.ceil((Tg.min() - epoch) / P), np.floor((Tg.max() - epoch) / P) + 1)
        ex = False
        for k in ks[ks != 0]:
            i = np.searchsorted(Tg, epoch + k * P)
            for ii in (i - 1, i):
                if 0 <= ii < len(Tg) and abs(Tg[ii] - (epoch + k * P)) < D / 8 and flat[ii]:
                    ex = True
        out.append(ex)
    return np.array(out)


def main():
    ov = pd.read_csv(OUT / "k2_overlap.csv")
    rows = []
    for r in ov.itertuples():
        D, delta = r.dur_h / 24, r.depth_ppm / 1e6
        prods = {"K2PDC": s3._segments(s3.k2lib.k2_lightcurves(int(r.EPIC))), "K2SFF": s3.load_hlsp(int(r.EPIC), "K2SFF")}
        for name, lc in prods.items():
            g = s3.k2lib.grid(lc, D, delta, r.epoch_btjd)
            st = status(g, r.epoch_btjd, D, pilot.P_GRID)
            row = dict(TOI=r.TOI, EPIC=r.EPIC, product=name, P_grid=f"{pilot.P_GRID.min():.1f}-{pilot.P_GRID.max():.0f}",
                       frac_logP_excluded=round(st.mean(), 3),
                       shortest_allowed_P=round(float(pilot.P_GRID[~st].min()), 2) if (~st).any() else np.nan,
                       frac_excl_P_lt_30=round(st[pilot.P_GRID < 30].mean(), 3),
                       frac_excl_P_30_100=round(st[(pilot.P_GRID >= 30) & (pilot.P_GRID < 100)].mean(), 3),
                       frac_excl_P_gt_100=round(st[pilot.P_GRID >= 100].mean(), 3))
            if r.TOI == 5560.01:
                ctl = status(g, r.epoch_btjd, D, np.array([140.494, 70.247]))
                row["control_140_70_excluded"] = " ".join(map(str, ctl))
            rows.append(row)
            print(row, flush=True)
    pd.DataFrame(rows).to_csv(OUT / "derived" / "k2_period_limits.csv", index=False)


if __name__ == "__main__":
    main()
