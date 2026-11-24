"""Check 4: RNAAS figure 1 (single column, greyscale-safe).

Panels: within-star demeaned ln(B1) and ln(flicker_rms) against demeaned disk_excess, main
crowd-cut WISE-usable sample (stars with >= 2 sectors, as stage 4). Faint markers: individual
sectors by disk phase; black squares: means in 10 quantile bins of x with standard errors;
line and band: pooled beta and its 68% bootstrap interval from results/full/test_results.csv
(main, k = 0), with p_perm and p_swap. Outputs figures/rnaas_fig1.png and .pdf.
"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import PILOT, FIGURES
import stats4, c1_phase_terms as C1

PANELS = [("B1", r"$\ln B_1$"), ("flicker_rms", r"$\ln\,$flicker rms")]
PHASES = [("building", ["building"], "^", "0.35"),
          ("plateau or quiescent", ["plateau", "quiescent"], "o", "0.65"),
          ("dissipating", ["dissipating"], "v", "0.15")]


def ptxt(p, n=10000):
    return f"< {1 / n:.0e}".replace("e-0", "e-") if p <= 1.5 / n else f"= {p:.3f}"


def main():
    sub, vsub, _, _ = C1.main_subset()
    tr = pd.read_csv(PILOT / "test_results.csv")
    tr = tr[(tr.variant == "main (crowd-cut)") & (tr.k == 0)].set_index("metric")
    sel = vsub["main (crowd-cut)"] & sub.wise_usable.values.astype(bool)
    plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "legend.fontsize": 6.5,
                         "xtick.labelsize": 7, "ytick.labelsize": 7})
    fig, axs = plt.subplots(2, 1, figsize=(3.4, 5.0), sharex=True)
    for ax, (metric, ylab), tag in zip(axs, PANELS, "ab"):
        y = np.log(sub[metric].values)
        ids, gid = stats4.compact(sub.tic.values, sel & np.isfinite(y))
        x, yy = sub.disk_excess.values[ids], y[ids]
        n = np.bincount(gid)
        xd = x - (np.bincount(gid, x) / n)[gid]
        yd = yy - (np.bincount(gid, yy) / n)[gid]
        lab = sub.phase_label.values[ids]
        r = tr.loc[metric]
        assert int(r.n_stars) == gid.max() + 1 and int(r.n_sectors) == len(ids), metric
        cnt = {}
        for name, labs, mk, c in PHASES:
            m = np.isin(lab, labs)
            cnt[name] = int(m.sum())
            ax.scatter(xd[m], yd[m], s=3, marker=mk, facecolors="none", edgecolors=c,
                       linewidths=0.4, alpha=0.6, label=f"{name} ({m.sum()})", rasterized=True)
        q = np.quantile(xd, np.linspace(0, 1, 11))
        b = np.clip(np.searchsorted(q, xd, side="right") - 1, 0, 9)
        bx = np.array([xd[b == i].mean() for i in range(10)])
        by = np.array([yd[b == i].mean() for i in range(10)])
        be = np.array([yd[b == i].std(ddof=1) / np.sqrt((b == i).sum()) for i in range(10)])
        g = np.linspace(np.percentile(xd, 0.5), np.percentile(xd, 99.5), 50)
        ax.fill_between(g, r.beta_lo * g, r.beta_hi * g, color="0.75", lw=0, zorder=3)
        ax.plot(g, r.beta * g, "k-", lw=1, zorder=4)
        ax.errorbar(bx, by, be, fmt="s", ms=3.5, c="k", mfc="white", mew=0.9, elinewidth=0.9,
                    capsize=1.5, zorder=5, label="10 quantile bins")
        ax.axhline(0, c="0.5", lw=0.4, ls=":")
        ax.set_ylim(np.percentile(yd, [1, 99]) * 1.1)
        ax.set_xlim(g[0], g[-1])
        ax.set_ylabel(ylab + " $-$ star mean")
        ax.text(0.03, 0.96, f"({tag}) $\\beta$ = {r.beta:.2f} [{r.beta_lo:.2f}, {r.beta_hi:.2f}]"
                f" mag$^{{-1}}$\n$p_{{\\rm perm}}$ {ptxt(r.p_perm)}, $p_{{\\rm swap}}$ {ptxt(r.p_swap)}",
                transform=ax.transAxes, va="top", fontsize=7,
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5), zorder=6)
    hs = [Line2D([], [], ls="none", marker=mk, mfc="none", mec=c, ms=4, mew=0.8,
                 label=f"{name} ({cnt[name]})") for name, _, mk, c in PHASES]
    hs.append(Line2D([], [], ls="none", marker="s", mfc="white", mec="k", ms=4, mew=0.9,
                     label="mean, 10 quantile bins"))
    hs.append(Line2D([], [], c="k", lw=1, label=r"pooled $\beta$, 68% band"))
    fig.legend(handles=hs, loc="upper center", ncol=2, frameon=False, fontsize=6.5,
               handletextpad=0.3, columnspacing=0.8, bbox_to_anchor=(0.55, 1.0))
    axs[1].set_xlabel("W1 disk excess $-$ star mean (mag)")
    n_st, n_sec = int(tr.loc["B1"].n_stars), int(tr.loc["B1"].n_sectors)
    fig.tight_layout(h_pad=0.4, rect=(0, 0, 1, 0.93))
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"rnaas_fig1.{ext}", dpi=300)
    print("binned means B1/flicker written;", n_st, "stars", n_sec, "sectors")


if __name__ == "__main__":
    main()
