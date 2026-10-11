"""Step 6: vetting diagnostics for every raw candidate.

For each candidate in candidates_raw.csv:
  * background / centroid excursion metrics during the event,
  * cross-star check: other hosts with 2-min data in the same sector AND camera
    (from cache/), flagged if any shows >=2 cadences above 3 sigma within +/-10 min,
  * a vetting plot vetting/<cand_id>_TIC<tic>_s<sec>.png.
Writes candidates_vetmetrics.csv. Human decisions live in manual_vetting.csv and
are merged by summarize.py.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import flarelib as F

WIN_X = 10.0 / 1440  # cross-star window


def robust_z(x, ref):
    med = np.nanmedian(ref)
    sd = 1.4826 * np.nanmedian(np.abs(ref - med))
    return (x - med) / sd if sd > 0 else np.zeros_like(x)


def cross_star(c, sectors):
    same = sectors[(sectors.sector == c.sector) & (sectors.camera == c.camera) & (sectors.tic != c.tic)]
    hits, maxz = [], []
    for tic in same.tic:
        z = np.load(F.CACHE / f"{tic}_s{c.sector:03d}.npz")
        sel = z["good"] & (np.abs(z["t"] - c.t_peak) <= WIN_X + (c.t_stop - c.t_peak))
        sel &= z["t"] >= c.t_start - WIN_X
        if not sel.any():
            maxz.append(np.nan)
            continue
        s = z["resid"][sel] / z["sigma"][sel]
        maxz.append(float(np.nanmax(s)))
        above = s > 3
        if np.any(above[:-1] & above[1:]):
            hits.append(int(tic))
    return len(same), hits, maxz


def main():
    cat = F.load_catalog()
    cand = pd.read_csv(F.SQ / "candidates_raw.csv")
    sectors = pd.read_csv(F.SQ / "sectors.csv")
    out_dir = F.SQ / "vetting"
    out_dir.mkdir(exist_ok=True)
    rows = []
    for _, c in cand.iterrows():
        lc = F.detrend(F.load_sector(int(c.tic), int(c.sector), cat))
        t = lc["t"]
        ev = (t >= c.t_start) & (t <= c.t_stop)
        loc = (np.abs(t - c.t_peak) < 0.5) & ~ev & np.isfinite(lc["bkg"])
        bkg_z = float(np.nanmax(np.abs(robust_z(lc["bkg"][ev], lc["bkg"][loc])))) if ev.any() else np.nan
        c1_z = float(np.nanmax(np.abs(robust_z(lc["c1"][ev], lc["c1"][loc]))))
        c2_z = float(np.nanmax(np.abs(robust_z(lc["c2"][ev], lc["c2"][loc]))))
        sap = lc["sap"] / np.nanmedian(lc["sap"][loc]) - 1
        sap_amp = float(np.nanmax(sap[ev]) - np.nanmedian(sap[loc]))
        n_cmp, hits, maxz = cross_star(c, sectors)
        dt_dump = float(np.min(np.abs(lc["dumps"] - c.t_peak)) * 24) if len(lc["dumps"]) else np.nan
        rows.append(dict(cand_id=c.cand_id, bkg_maxz=round(bkg_z, 2), cen1_maxz=round(c1_z, 2),
                         cen2_maxz=round(c2_z, 2), sap_amp=round(sap_amp, 5), pdc_amp=round(c.amp, 5),
                         hours_to_nearest_dump=round(dt_dump, 2), n_compare_stars=n_cmp,
                         compare_hits=" ".join(map(str, hits)),
                         compare_max_sigma=round(float(np.nanmax(maxz)), 2) if len(maxz) and np.isfinite(maxz).any() else np.nan))

        # ------------------------------------------------------------- plot
        fig, ax = plt.subplots(6, 2, figsize=(13, 13), sharex="col",
                               gridspec_kw=dict(width_ratios=[2, 1]))
        for col, half in enumerate((0.6, 0.06)):
            w = np.abs(t - c.t_peak) < half
            g = w & lc["good"]
            a = ax[0, col]
            a.plot(t[g], lc["f"][g] / lc["med"], ".", ms=2, color="k", label="PDCSAP 2-min")
            a.plot(t[g], lc["trend"][g] / lc["med"], "-", color="C1", lw=1, label="trend")
            if col == 1 and int(c.sector) in F.spoc_files(int(c.tic), 20):
                l20 = F.load_sector(int(c.tic), int(c.sector), cat, exptime=20)
                w20 = (np.abs(l20["t"] - c.t_peak) < half) & l20["good"]
                a.plot(l20["t"][w20], l20["f"][w20] / l20["med"], ".", ms=1.5, color="C0", alpha=0.5,
                       label="PDCSAP 20-s", zorder=0)
            a.axvspan(c.t_start, c.t_stop, color="C3", alpha=0.15)
            a.legend(fontsize=7, loc="upper right")
            a.set_ylabel("norm. flux")
            a = ax[1, col]
            a.plot(t[g], lc["resid"][g] / lc["sigma"][g], ".", ms=2, color="k")
            a.axhline(3, color="C3", lw=0.8, ls="--")
            a.set_ylabel("resid / sigma")
            a = ax[2, col]
            a.plot(t[w], lc["sap"][w] / np.nanmedian(lc["sap"][w]), ".", ms=2, color="C2")
            a.set_ylabel("SAP (norm)")
            a = ax[3, col]
            a.plot(t[w], lc["bkg"][w], ".", ms=2, color="C4")
            a.set_ylabel("SAP_BKG")
            a = ax[4, col]
            a.plot(t[w], lc["c1"][w] - np.nanmedian(lc["c1"][w]), ".", ms=2, label="col")
            a.plot(t[w], lc["c2"][w] - np.nanmedian(lc["c2"][w]), ".", ms=2, label="row")
            a.legend(fontsize=7)
            a.set_ylabel("centroid (pix)")
            a = ax[5, col]
            q = lc["q"]
            for bit, lab, mk in ((F.Q_DESAT, "32 dump", "v"), (F.Q_IMPULSIVE_OUTLIER, "512 impulsive", "o"),
                                 (F.Q_APERTURE_COSMIC | F.Q_COLLATERAL_COSMIC, "CR 64/1024", "x"),
                                 (F.BAD_BITS & ~F.Q_DESAT, "other bad", "s")):
                m = w & ((q & bit) > 0)
                a.plot(t[m], np.full(m.sum(), {32: 3, 512: 2}.get(bit, 1 if lab.startswith("CR") else 0)),
                       mk, ms=4, label=lab)
            a.set_yticks([0, 1, 2, 3], ["bad", "CR", "512", "dump"])
            a.set_ylim(-0.5, 3.5)
            a.set_xlabel("BTJD")
            for r in range(6):
                ax[r, col].axvline(c.t_peak, color="C3", lw=0.5, alpha=0.6)
                for d in lc["dumps"][np.abs(lc["dumps"] - c.t_peak) < half]:
                    ax[r, col].axvline(d, color="C0", lw=0.6, ls=":")
                for (s0, s1) in lc["segs"]:
                    for e in (s0 + F.EDGE_D, s1 - F.EDGE_D):
                        if abs(e - c.t_peak) < half:
                            ax[r, col].axvline(e, color="grey", lw=0.8, ls="-.")
        hit = rows[-1]["compare_hits"] or "none"
        fig.suptitle(f"{c.cand_id}  TIC {c.tic}  S{c.sector} cam{c.camera}  peak {c.t_peak:.4f}  "
                     f"amp {c.amp*1e3:.2f} ppt ({c.amp_sig:.1f} sig)  ED {c.ed_s:.1f} s  dur {c.dur_min:.0f} min\n"
                     f"bkg z={bkg_z:.1f}  cen z={c1_z:.1f}/{c2_z:.1f}  SAP amp {sap_amp*1e3:.2f} ppt  "
                     f"cross-star: {n_cmp} compared, hits: {hit}  nearest dump {dt_dump:.2f} h  "
                     f"flag512 pts={c.n_flag512}", fontsize=9)
        fig.tight_layout()
        fig.savefig(out_dir / f"{c.cand_id}_TIC{c.tic}_s{c.sector:03d}.png", dpi=90)
        plt.close(fig)
    m = pd.DataFrame(rows)
    m.to_csv(F.SQ / "candidates_vetmetrics.csv", index=False)
    print(m.to_string())


if __name__ == "__main__":
    main()
