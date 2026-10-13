"""
Figure for the EPIC 202126877 Research Note (replaces depth_history_summary.png in the note; review_fable R8).
It shows measured values only, plus simple trends; the rejected constant-rate and quadratic model curves are gone.
  Top: primary and secondary eclipse depths per eclipse (results_depth_events.csv; median flux within +-15 min,
       own TESScut r = 1.5 px), the blind analysis's season depths (independent/results.csv, model depths,
       r = 2 px), the 2004-08 SuperWASP+KELT 3-sigma primary limits (worst 4.49%, typical 1.61%;
       19_ground_combined_limits.py) and the K2 secondary limit (< 0.1%).
  Bottom: inclination per season from our model (base value; bar = range over the systematic variants,
       data/inclination_trend.csv) and from the blind analysis (bar = stat and syst in quadrature), the
       eclipse thresholds (primary 84.4, secondary 86.5 deg), the totality band (89.2-89.9 deg), and linear
       fits to the three TESS seasons of each analysis, extrapolated as dashed lines.
Output: eb_depth_evolution/rnaas/eb_history.png (and a copy in figures/)
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

D = Path(__file__).resolve().parents[1]
BLUE, GREEN, GREY, ORANGE = "#1f6fb4", "#2c9a4b", "#6b7a8c", "#d9682b"


def main():
    ev = pd.read_csv(D / "results_depth_events.csv")
    bl = pd.read_csv(D / "independent" / "results.csv")
    tr = pd.read_csv(D / "data" / "inclination_trend.csv")

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 7.4), sharex=True, gridspec_kw={"height_ratios": [1, 1.05]})

    # --- top: depths
    for kind, col, lab in (("P", BLUE, "primary"), ("S", GREEN, "secondary")):
        e = ev[ev.kind == kind]
        a1.plot(e.year, 100 * e.depth_modelfree, "o", color=col, ms=6, label=f"{lab} (this work)")
    a1.plot(bl.year, 100 * bl.DP_obs, "D", mfc="none", color=BLUE, ms=8, label="primary (blind re-analysis)")
    sec = bl[bl.DS_observed.astype(str) == "True"]
    a1.plot(sec.year, 100 * sec.DS_obs, "D", mfc="none", color=GREEN, ms=8, label="secondary (blind re-analysis)")
    for x, lim in ((2005.0, 4.49), (2007.5, 1.61)):
        a1.errorbar(x, lim, yerr=1.5, uplims=True, fmt="v", color=GREY, ms=7, capsize=0)
    a1.text(2006.2, 6.6, "SuperWASP + KELT\n3$\\sigma$ limits", ha="center", fontsize=8, color=GREY)
    a1.errorbar(2014.36, 0.1, yerr=1.5, uplims=True, fmt="v", color=GREEN, ms=7, capsize=0)
    a1.text(2014.6, 2.3, "K2 secondary < 0.1%", fontsize=8, color=GREEN)
    a1.set_ylabel("eclipse depth (%)")
    a1.set_ylim(-1, 42)
    a1.legend(fontsize=7.5, loc="upper left", frameon=False)

    # --- bottom: inclination
    years = np.array([2014.36, 2021.81, 2023.87, 2026.04])
    cols = ["i2014", "i2021", "i2023", "i2026"]
    base = tr[(tr.variant == "base") & (np.isclose(tr.ph14, 0.3322))].iloc[0]
    ours = np.array([base[c] for c in cols])
    lo = np.array([tr[c].min() for c in cols])
    hi = np.array([tr[c].max() for c in cols])
    a2.errorbar(years, ours, yerr=[ours - lo, hi - ours], fmt="o", color="k", ms=6, capsize=3,
                label="this work (bar: model variants)")
    bi = bl.i_deg.to_numpy()
    be = np.hypot(bl.i_stat, bl.i_syst).to_numpy()
    a2.errorbar(bl.year + 0.25, bi, yerr=be, fmt="D", mfc="none", color=ORANGE, ms=7, capsize=3,
                label="blind re-analysis (bar: stat+syst)")
    t = np.linspace(2004, 2032.5, 200)
    for yy, ii, col, name in ((years[1:], ours[1:], "k", "this work"), (bl.year.to_numpy()[1:], bi[1:], ORANGE, "blind")):
        c = np.polyfit(yy, ii, 1)
        a2.plot(t, np.polyval(c, t), "--", color=col, lw=1, alpha=0.8,
                label=f"TESS-era linear fit, {name}: {c[0]:.2f}$^\\circ$/yr")
    a2.axhline(84.4, color=BLUE, lw=0.8, ls=":")
    a2.text(2004.3, 84.5, "primary eclipses begin (i > 84.4$^\\circ$)", fontsize=7.5, color=BLUE)
    a2.axhline(86.5, color=GREEN, lw=0.8, ls=":")
    a2.text(2004.3, 86.6, "secondary eclipses begin (i > 86.5$^\\circ$)", fontsize=7.5, color=GREEN)
    a2.axhspan(89.2, 89.9, color=GREY, alpha=0.18)
    a2.text(2004.3, 89.35, "primary total (89.2-89.9$^\\circ$)", fontsize=7.5, color=GREY)
    a2.axvspan(2004, 2008.3, color=GREY, alpha=0.07)
    a2.set_ylabel("orbital inclination (deg)")
    a2.set_xlabel("year")
    a2.set_ylim(83.6, 90.3)
    a2.set_xlim(2003.5, 2032.5)
    a2.legend(fontsize=7.2, loc="lower right", frameon=False)
    fig.tight_layout()
    for p in (D / "rnaas" / "eb_history.png", D / "figures" / "eb_history.png"):
        fig.savefig(p, dpi=200)
    print("TESS-era slopes:",
          np.polyfit(years[1:], ours[1:], 1)[0].round(3), np.polyfit(bl.year.to_numpy()[1:], bi[1:], 1)[0].round(3))


if __name__ == "__main__":
    main()
