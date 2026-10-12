"""
Independent re-derivation of the O-C signal for the LTTE candidates (CLAUDE.md rule: a
claim must be re-derived by a method that does not share the first analysis' code).

Shares only the light-curve loader. No template, no per-eclipse fits, no iteration:
  - linear ephemeris = catalogue (Prsa+2022) P and T0, period corrected by straight-line
    fits to the per-sector results below (first with a wide +-0.2 phase window, then the
    narrow one), so only the ephemeris, not the eclipse shape, is iterated;
  - per sector: fold on that ephemeris, take points within +-0.6 eclipse widths
    (catalogue Wp-2g) of the primary minimum, bin the fold (0.002 phase), and fit a parabola
    to the binned medians around the lowest bin; the vertex gives that sector's primary O-C
    (error: 20 bootstrap resamples). Same for the secondary at catalogue phase Phis-2g.
  - compare with the sector-mean O-C from run_etv.py (results/oc/<TIC>_oc.csv).
Usage:  python sidequests/eb_timing/scripts/independent_check.py TIC [TIC ...]
Output: sidequests/eb_timing/results/independent_check.csv and figs/<TIC>_indep.png
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE.parents[1] / "scripts"))
import alias_pipeline as ap  # noqa: E402

ap.LC = HERE / "cache" / "lc"


def vertex(x, y, w):
    """Binned-median fold (bins of 0.002 phase), lowest bin, then a parabola through the
    binned medians within +-0.25 w of it; returns the vertex phase."""
    b = np.floor(x / 0.002).astype(int)
    med = pd.Series(y).groupby(b).median()
    xb = (med.index.values + 0.5) * 0.002
    k = np.argmin(med.values)
    m = np.abs(xb - xb[k]) < max(0.25 * w, 0.006)
    if m.sum() < 4:
        return np.nan
    c = np.polyfit(xb[m] - xb[k], med.values[m], 2)
    if c[0] <= 0:
        return np.nan
    return xb[k] - c[1] / (2 * c[0])


def sector_times(lc, P, T0, ph_ecl, w, win=None):
    rng = np.random.default_rng(0)
    rows = []
    for s, g in lc.groupby("sector"):
        ph = ((g.time.values - T0) / P - ph_ecl + 0.5) % 1 - 0.5
        m = np.abs(ph) < (0.6 * w if win is None else win)
        if m.sum() < 30:
            continue
        x, y = ph[m], g.flux.values[m]
        v = vertex(x, y, w)
        if not np.isfinite(v):
            continue
        boot = [vertex(x[i], y[i], w) for i in (rng.integers(0, len(x), len(x)) for _ in range(20))]
        rows.append(dict(sector=s, t=g.time.mean(), oc=v * P, err=np.nanstd(boot) * P + 1e-6))
    return pd.DataFrame(rows)


def main():
    cat = pd.read_csv(HERE / "data" / "prsa2022_tess_ebs.csv").drop_duplicates("TIC").set_index("TIC")
    out = []
    for tic in [int(x) for x in sys.argv[1:]]:
        r = cat.loc[tic]
        P, T0 = r.Per, r.BJD0
        lc = ap.load_lightcurves(tic, products=["SPOC_120", "TESS-SPOC", "QLP"])
        fig, ax = plt.subplots(figsize=(8, 4))
        mine = pd.read_csv(HERE / "results" / "oc" / f"{tic}_oc.csv")
        mine = mine[mine.flag.isna() | (mine.flag == "")]
        for kind, ph_ecl, w, col in (("p", 0.0, r["Wp-2g"], "C0"), ("s", r["Phis-2g"], r["Ws-2g"], "C3")):
            if not np.isfinite(w) or not np.isfinite(ph_ecl):
                continue
            # pass 1: wide window (+-0.2 phase) to catch a catalogue period error; pass 2-3:
            # narrow window with the period corrected by a line through the previous pass
            P2 = P
            for win in (0.2, None, None):
                st = sector_times(lc, P2, T0, ph_ecl, w, win)
                if len(st) < 5:
                    break
                c = np.polyfit(st.t - T0, st.oc, 1)
                P2 = P2 * (1 + c[0])
            if len(st) < 5:
                continue
            c = np.polyfit(st.t, st.oc, 1, w=1 / st.err)
            st["oc_detr"] = st.oc - np.polyval(c, st.t)
            mk = mine[mine.kind == kind]
            sm = mk.groupby("sector").oc.mean()
            cm = np.polyfit(mk.groupby("sector")["T"].mean().values, sm.values, 1)
            ref = (sm - np.polyval(cm, mk.groupby("sector")["T"].mean())).rename("oc_pipeline")
            st = st.join(ref, on="sector")
            ok = st.oc_pipeline.notna()
            corr = np.corrcoef(st.oc_detr[ok], st.oc_pipeline[ok])[0, 1]
            rms = np.std((st.oc_detr - st.oc_pipeline)[ok]) * 1440
            ax.errorbar(st.t, st.oc_detr * 1440, st.err * 1440, fmt="o", color=col, ms=3, label=f"{kind} indep")
            ax.plot(st.t[ok], st.oc_pipeline[ok] * 1440, "x", color="k", ms=4)
            out.append(dict(tic=tic, kind=kind, n_sectors=int(ok.sum()), corr=round(corr, 3),
                            rms_diff_min=round(rms, 2), ptp_indep_min=round(np.ptp(st.oc_detr) * 1440, 1),
                            ptp_pipeline_min=round(np.ptp(st.oc_pipeline[ok]) * 1440, 1)))
            print(out[-1], flush=True)
        ax.set_title(f"TIC {tic}: independent per-sector fold-vertex timing (o) vs pipeline sector means (x)", fontsize=9)
        ax.set_xlabel("BTJD")
        ax.set_ylabel("O-C, linear removed (min)")
        ax.legend(fontsize=8)
        fig.savefig(HERE / "figs" / f"{tic}_indep.png", dpi=90, bbox_inches="tight")
        plt.close(fig)
    pd.DataFrame(out).to_csv(HERE / "results" / "independent_check.csv", index=False)


if __name__ == "__main__":
    main()
