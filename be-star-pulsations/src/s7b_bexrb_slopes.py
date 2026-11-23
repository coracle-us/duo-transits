"""Stage 7b: per-star slopes for the Be X-ray binaries found by s7_bexrb_xmatch.py.

For each matched star: within-star OLS slope of ln(metric) on disk_excess over its
wise_usable sectors (no crowding cut, so all four are shown; the 'main' column says which
are in the main sample). Context: the percentile of each slope among all main-sample stars
with >= 2 usable sectors, and the fraction of those stars with a negative A1 slope.
Then the stage-4 A1 test for the main sample with the matched stars removed.
Outputs: results/<run>/bexrb_slopes.csv, printed summary.
"""
import numpy as np, pandas as pd
from common import PILOT, SEED
import s4_test as T, stats4

METRICS = ["A1", "B1", "red_noise", "flicker_rms"]


def per_star_slopes(sm, keep, metric):
    d = sm[keep & sm.wise_usable.astype(bool)].copy()
    d["y"] = np.log(d[metric])
    d = d[np.isfinite(d.y) & np.isfinite(d.disk_excess)]
    out = {}
    for tic, g in d.groupby("tic"):
        if len(g) < 2:
            continue
        x = g.disk_excess.values - g.disk_excess.mean()
        sxx = np.sum(x * x)
        if sxx <= 0:
            continue
        out[tic] = dict(slope=float(np.sum(x * (g.y.values - g.y.mean())) / sxx),
                        n_sec=len(g), x_range=float(np.ptp(g.disk_excess.values)))
    return pd.DataFrame(out).T


def main():
    sm, sw, disks, dos = T.load()
    vs = T.variants(sm)
    xm = pd.read_csv(PILOT / "bexrb_xmatch.csv")
    hits = xm[xm.in_catalogue.fillna(False).astype(bool)]
    tics = hits.tic.astype(int).tolist()
    main_mask = vs["main (crowd-cut)"]
    all_mask = np.ones(len(sm), bool)
    rows = []
    for m in METRICS:
        ref = per_star_slopes(sm, main_mask, m)
        allst = per_star_slopes(sm, all_mask, m)
        frac_neg = float(np.mean(ref.slope < 0))
        print(f"{m}: main-sample stars with >=2 usable sectors {len(ref)}, "
              f"negative slope {frac_neg:.2f}, median slope {ref.slope.median():.2f}")
        for _, h in hits.iterrows():
            t = int(h.tic)
            if t not in allst.index:
                continue
            s = allst.loc[t]
            rows.append(dict(name=h["name"], tic=t, metric=m, slope=s.slope,
                             n_sec=int(s.n_sec), x_range=s.x_range,
                             in_main=bool(main_mask[sm.tic.values == t].any()),
                             pct_in_main=float(np.mean(ref.slope < s.slope) * 100),
                             main_frac_negative=frac_neg))
    out = pd.DataFrame(rows)
    out.to_csv(PILOT / "bexrb_slopes.csv", index=False)
    print(out.to_string(index=False))
    # stage-4 A1 test without the matched stars
    sid, t, x = sm.tic.values, sm.mid_mjd.values, sm.disk_excess.values
    lab, wu = sm.phase_label.values.astype(object), sm.wise_usable.values.astype(bool)
    sel = main_mask & ~np.isin(sid, tics)
    y = np.where(sel, np.log(sm.A1.values), np.nan)
    r = stats4.run_test(sid, x, y, lab, wu, t, disks, dos, T.evaluate, n_null=10000,
                        n_boot=2000, seed=SEED)
    print("A1, main sample without matched BeXRBs:",
          {k: r.get(k) for k in ("n_stars", "beta", "beta_lo", "beta_hi", "p_perm", "p_swap")})


if __name__ == "__main__":
    main()
