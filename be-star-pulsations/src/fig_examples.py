"""figures/pilot_examples.png: 6 stars, WISE W1 + spline with TESS sectors, and the sectors'
amplitude spectra stacked and coloured by disk_excess. Picks the 6 main-sample stars with
the widest disk_excess range among their usable sectors."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from scipy.interpolate import splev
from common import PILOT, FIGURES, WISE_DIR, SPEC_DIR, figpath
from wise import DiskCurve
from metrics import LOG_GRID


def main():
    sm = pd.read_csv(PILOT / "sector_metrics.csv")
    u = sm[sm.usable.astype(bool) & sm.crowd_ok.fillna(False).astype(bool)]
    g = u.groupby("tic").disk_excess
    rngx = (g.max() - g.min())[g.count() >= 2].sort_values(ascending=False)
    tics = rngx.index[:6]
    norm = colors.Normalize(0, max(0.05, u[u.tic.isin(tics)].disk_excess.max()))
    cmap = cm.viridis
    fig, axs = plt.subplots(6, 2, figsize=(12, 17), gridspec_kw=dict(width_ratios=[1.3, 1]))
    for (axw, axs_), tic in zip(axs, tics):
        v = pd.read_csv(WISE_DIR / "visits" / f"{tic}.csv")
        d = DiskCurve(v.mjd.values, v.w1mpro_ep.values)
        tt = np.linspace(d.t0, d.t1, 1500)
        axw.errorbar(v.mjd, v.w1mpro_ep, v.w1sigmpro_ep, fmt="o", ms=3, c="0.4")
        axw.plot(tt, splev(tt, d.tck), "k-", lw=1)
        ss = sm[sm.tic == tic].sort_values("mid_mjd")
        off = 0
        for _, s in ss.iterrows():
            ok = bool(s.usable) and bool(s.crowd_ok) if pd.notna(s.crowd_ok) else False
            c = cmap(norm(s.disk_excess)) if ok else "0.8"
            axw.axvline(s.mid_mjd, color=c, lw=2, alpha=0.9)
            if ok:
                a = np.load(SPEC_DIR / f"{tic}_{s.sector}.npy")
                axs_.plot(LOG_GRID, a + off, color=c, lw=1)
                axs_.text(24, off, f" S{s.sector} {s.phase_label[:4]}", fontsize=6, va="bottom")
                off += np.nanmax(a) * 0.8
        axw.invert_yaxis()
        axw.set_ylabel("W1 (mag)")
        axw.set_title(f"TIC {tic}: disk_excess range {rngx[tic]:.2f} mag "
                      f"({ss['product'].iloc[0]})", fontsize=9)
        axs_.set_xscale("log")
        axs_.set_ylabel("amplitude (ppt), offset")
    axs[-1, 0].set_xlabel("MJD (vertical lines: TESS sectors, grey = not usable)")
    axs[-1, 1].set_xlabel("frequency (c/d)")
    sm_ = cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm_, ax=axs[:, 1], label="disk_excess (mag)", shrink=0.4)
    fig.savefig(figpath("examples"), dpi=100, bbox_inches="tight")


if __name__ == "__main__":
    main()
