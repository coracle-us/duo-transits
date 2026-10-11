"""
Diagnostics for one candidate's alias result: what was detected, what excluded what.

Runs the alias pipeline on all TESS data, then
  - lists robust detected dips (the DETECT_* cuts) and every supported epoch, with the
    aliases each one is consistent with;
  - for each alias, the first epoch that excluded it;
  - plots the light curve around Tc1, Tc2 and each detected dip (+ a few excluding windows),
    with the expected box transit (depth, duration from Tc1/Tc2) overplotted.

Usage:  python scripts/diagnose_candidate.py TIC [TIC ...]
Output: results/diagnostics/<TIC>.png, results/diagnostics/<TIC>.txt
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import alias_pipeline as ap

ROOT = Path(__file__).resolve().parent.parent
CAT = ROOT / "data" / "catalog"
OUT = ROOT / "results" / "diagnostics"


def detected_mask(ep, delta):
    side_min = ap.DETECT_SIDE * ap.BASE_W / ep.cadence
    return ((ep.verdict == "supported") & (ep.snr >= ap.DETECT_SNR)
            & (ep.depth <= ap.DEPTH_MAX_FRAC * delta) & (ep.cover >= ap.DETECT_COVER)
            & (ep.n_oot_left >= side_min) & (ep.n_oot_right >= side_min))


def diagnose(tic, cand):
    c = cand.loc[tic]
    lc = ap.load_lightcurves(tic)
    ep, al = ap.run_candidate(c, lc)
    D = c.Dur_h / 24
    delta = 1e-3 * np.nanmean([c.Depth1_ppt, c.Depth2_ppt])
    dT = c.Tc2 - c.Tc1
    lines = [f"TIC {tic}: dT={dT:.4f} d, depth={delta*1e3:.2f} ppt, dur={c.Dur_h:.2f} h, "
             f"N_alias(paper)={c.N_alias}, Tmag={c.Tmag}, b={c.b}",
             "aliases: " + str(al.verdict.value_counts().to_dict())]

    # unique dips (merge epochs tested by several aliases)
    ep["Tobs"] = ep.T_pred + ep.offset.fillna(0)
    sup = ep[ep.verdict == "supported"].copy()
    sup["robust"] = detected_mask(sup, delta)
    sup = sup.sort_values("snr", ascending=False).drop_duplicates(subset=["T_pred"])
    dips = []
    for _, r in sup.sort_values("Tobs").iterrows():
        if dips and abs(r.Tobs - dips[-1]["Tobs"]) < D:
            if r.snr > dips[-1]["snr"]:
                dips[-1] = r.to_dict()
            continue
        dips.append(r.to_dict())
    lines.append(f"\nsupported dips (merged): {len(dips)}")
    for d in dips:
        x = (d["Tobs"] - c.Tc1) / dT          # position in units of dT
        consistent = [n for n in range(1, int(dT / ap.P_FLOOR) + 1)
                      if abs((d["Tobs"] - c.Tc1) / (dT / n) - round((d["Tobs"] - c.Tc1) / (dT / n)))
                      * dT / n < D / 2 + 3 * d["sigma_T"]]
        lines.append(f"  T={d['Tobs']:.3f} s{d['sector']} depth={d['depth']*1e3:.2f}+/-"
                     f"{d['depth_err']*1e3:.2f} ppt snr={d['snr']:.1f} cover={d['cover']:.2f} "
                     f"beta={d['beta']:.2f} robust={d['robust']} (T-Tc1)/dT={x:.4f} "
                     f"consistent n: {consistent[:12]}{' ...' if len(consistent) > 12 else ''}")

    exc = ep[ep.verdict == "excluded"]
    lines.append("\nexcluding windows used most often (time, sector, #aliases killed):")
    cnt = exc.groupby(["T_pred", "sector"]).size().sort_values(ascending=False).head(10)
    for (T, s), k in cnt.items():
        r = exc[exc.T_pred == T].iloc[0]
        lines.append(f"  T={T:.3f} s{s} kills {k} aliases: depth={r.depth*1e3:.2f}+/-"
                     f"{r.depth_err*1e3:.2f} ppt beta={r.beta:.2f} cover={r.cover:.2f}")
    txt = "\n".join(lines)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{tic}.txt").write_text(txt, encoding="utf-8")
    print(txt)

    # ---- plot
    wins = [("Tc1", c.Tc1), ("Tc2", c.Tc2)] + [(f"dip s{d['sector']}{'*' if d['robust'] else ''}",
                                               d["Tobs"]) for d in dips
                                              if min(abs(d["Tobs"] - c.Tc1), abs(d["Tobs"] - c.Tc2)) > D]
    wins += [(f"excl s{s}", T) for (T, s) in list(cnt.index)[:3]]
    wins = wins[:12]
    ncol = 4
    nrow = int(np.ceil(len(wins) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4 * ncol, 2.8 * nrow), squeeze=False)
    t, f = lc.time.to_numpy(), lc.flux.to_numpy()
    for ax, (lab, T) in zip(axs.flat, wins):
        m = np.abs(t - T) < max(1.5 * D, 0.6)
        if m.sum():
            tt, ff = t[m], f[m]
            oot = np.abs(tt - T) > D / 2 + 0.05
            if oot.sum() > 3:
                deg = 1 if (np.any(tt[oot] < T) and np.any(tt[oot] > T)) else 0
                ff = ff / np.polyval(np.polyfit(tt[oot] - T, ff[oot], deg), tt - T)
            ax.plot((tt - T) * 24, (ff - 1) * 1e3, ".", ms=2, color="0.5")
            # 30-min bins
            bins = np.arange(tt.min(), tt.max() + 1 / 48, 1 / 48)
            idx = np.digitize(tt, bins)
            bx = [tt[idx == i].mean() for i in np.unique(idx)]
            by = [ff[idx == i].mean() for i in np.unique(idx)]
            ax.plot((np.array(bx) - T) * 24, (np.array(by) - 1) * 1e3, "o", ms=3, color="C0")
        ax.plot([-D * 12, -D * 12, D * 12, D * 12], [0, -delta * 1e3, -delta * 1e3, 0], "r-", lw=1)
        sec = lc.sector[np.argmin(np.abs(t - T))] if len(t) else "?"
        prod = lc["product"][np.argmin(np.abs(t - T))] if len(t) else "?"
        ax.set_title(f"{lab} T={T:.2f} ({prod})", fontsize=8)
        ax.set_ylim(-2.5 * delta * 1e3, 1.5 * delta * 1e3)
        ax.tick_params(labelsize=7)
    for ax in list(axs.flat)[len(wins):]:
        ax.axis("off")
    fig.suptitle(f"TIC {tic}  (hours from window centre; ppt)  red: expected transit; * = robust dip",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / f"{tic}.png", dpi=90)
    plt.close(fig)


def main():
    cand = pd.read_csv(CAT / "candidates.csv").set_index("TIC", drop=False)
    for tic in map(int, sys.argv[1:]):
        diagnose(tic, cand)
        print()


if __name__ == "__main__":
    main()
