"""Stage 4b: leave-one-star-out jackknife of beta (leverage check) for the pilot."""
import numpy as np, pandas as pd
import s4_test as T, stats4
from common import PILOT

sm, sw, disks, dos = T.load()
rows = []
for vn, sel in list(T.variants(sm).items())[:2]:
    for m in ["A1", "B1", "B2", "B3", "red_noise", "flicker_rms"]:
        y = np.where(sel & sm.wise_usable.values.astype(bool), np.log(sm[m].values), np.nan)
        ids, gid = stats4.compact(sm.tic.values, np.isfinite(y))
        x, yy, tic = sm.disk_excess.values[ids], y[ids], sm.tic.values[ids]
        lc = stats4.code(sm.phase_label.values[ids])
        full, _, _, _ = stats4._stats(gid, x, yy, lc)
        jk = {}
        for t in np.unique(tic):
            k = tic != t
            _, g2 = np.unique(gid[k], return_inverse=True)
            jk[t] = stats4._stats(g2, x[k], yy[k], lc[k])[0]
        s = pd.Series(jk)
        # median-of-star slopes (each star's own LS slope), a leverage-robust alternative
        sl = []
        for t in np.unique(tic):
            k = tic == t
            if np.ptp(x[k]) > 0.02:
                sl.append(np.polyfit(x[k], yy[k], 1)[0])
        rows.append(dict(variant=vn, metric=m, beta=full, jk_min=s.min(), jk_min_drop=s.idxmin(),
                         jk_max=s.max(), jk_max_drop=s.idxmax(), n_star_slopes=len(sl),
                         median_star_slope=np.median(sl), frac_star_slopes_pos=np.mean(np.array(sl) > 0)))
out = pd.DataFrame(rows)
out.to_csv(PILOT / "jackknife.csv", index=False)
print(out.round(3).to_string())
